"""Batched generation of synthetic clusters (and of their log-likelihoods).

Once a :py:class:`asteca.Synthetic` object is created, every random value used by
:py:func:`synth_model.get_model` is drawn in advance, so a model is a deterministic
function of its parameters. The compiled functions in this module loop over a whole
batch of parameter sets in a single call, chaining the numba kernels of
:py:mod:`synth_model` in the same order and with the same operations as
:py:func:`synth_model.get_model`, so that each synthetic cluster is identical to the
one returned by :py:meth:`asteca.Synthetic.generate`. This removes the Python overhead
of each call (which dominates for clusters with few stars), and the batch of
log-likelihoods can be evaluated in parallel threads.

Only available with numba; the :py:class:`asteca.Synthetic` methods fall back to a
loop over :py:meth:`asteca.Synthetic.generate` otherwise.
"""

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ..synthetic import Synthetic

import numpy as np

from . import HAS_NUMBA as _HAS_NUMBA

# Order of the full parameter vector used by the compiled functions
FULL_PARAMS = ("met", "loga", "alpha", "beta", "Av", "DR", "Rv", "dm")


def available(ext_law: str) -> bool:
    """``True`` if the compiled batch functions can be used"""
    return _HAS_NUMBA and ext_law in ("GAIADR3", "CCMO")


def check_n_threads(n_threads) -> None:
    """Raise a ``ValueError`` if ``n_threads`` is not an integer between ``1`` and
    ``NUMBA_NUM_THREADS``. Only called when the threads are actually used (numba
    installed and ``parallel=True``)."""
    import numba

    if (
        isinstance(n_threads, bool)
        or not isinstance(n_threads, (int, np.integer))
        or n_threads < 1
    ):
        raise ValueError(f"'n_threads' must be a positive integer, got {n_threads!r}")
    # NUMBA_NUM_THREADS is set dynamically in numba.config (unknown to type checkers)
    max_threads = numba.config.NUMBA_NUM_THREADS  # pyright: ignore[reportPrivateImportUsage, reportAttributeAccessIssue]
    if n_threads > max_threads:
        raise ValueError(
            f"'n_threads={n_threads}' can not be larger than NUMBA_NUM_THREADS="
            + f"{max_threads}"
        )


def ccmo_array(ext_law: str, ext_coefs: list) -> np.ndarray:
    """CCMO extinction coefficients as an array for the compiled functions: one row
    ``[a, b, 0, 0]`` for the magnitude and one row ``[a1, b1, a2, b2]`` for each
    color (see :py:func:`synth_model.ccmo_coefs`). Empty for other extinction laws.

    :param ext_law: Extinction law
    :type ext_law: str
    :param ext_coefs: List of extinction coefficients, ``[ec_mag, ec_col1, ...]``
    :type ext_coefs: list

    :return: Array of shape ``(N_colors + 1, 4)``, or ``(0, 4)``
    :rtype: np.ndarray
    """
    if ext_law != "CCMO":
        return np.zeros((0, 4))
    rows = [[ext_coefs[0][0], ext_coefs[0][1], 0.0, 0.0]]
    for (a1, b1), (a2, b2) in ext_coefs[1:]:
        rows.append([a1, b1, a2, b2])
    return np.array(rows, dtype=float)


def batch_data(synthcl: "Synthetic") -> dict:
    """Arrays used by the compiled functions that do not change after the
    :py:class:`asteca.Synthetic` object is created. The typed lists hold references
    to the existing IMF samples and mass indexes (no copies).

    :param synthcl: :py:class:`asteca.Synthetic` object
    :type synthcl: Synthetic

    :return: Dictionary with the data
    :rtype: dict
    """
    from numba.typed.typedlist import List

    N_z, N_a = synthcl.theor_tracks.shape[:2]
    st_list = List()
    for ml in range(N_z):
        for al in range(N_a):
            st_list.append(synthcl.st_dist_mass[ml][al])
    # Indexes for the (ml, al) box and its closest node (m0, a0) are stored at
    # position ((ml * N_a + al) * 2 + (m0 - ml)) * 2 + (a0 - al). Nodes beyond the
    # grid's edges are never used (repeated keys keep the positions consistent)
    midx_list = List()
    for ml in range(N_z):
        for al in range(N_a):
            for dz in (0, 1):
                for da in (0, 1):
                    key = (ml, al, min(ml + dz, N_z - 1), min(al + da, N_a - 1))
                    midx_list.append(synthcl.mass_idx[key])

    return {
        "met_grid": np.asarray(synthcl.met_age_dict["met"], dtype=float),
        "age_grid": np.asarray(synthcl.met_age_dict["loga"], dtype=float),
        "N_a": N_a,
        "st_list": st_list,
        "midx_list": midx_list,
        "dr_frac": np.ascontiguousarray(synthcl.rand_floats["dr_frac"]),
        "ext_kind": 0 if synthcl.ext_law == "GAIADR3" else 1,
        "ccm": ccmo_array(synthcl.ext_law, synthcl.ext_coefs),
        "unif1": np.ascontiguousarray(synthcl.rand_floats["unif"][1]),
    }


def full_params(
    params: dict,
    def_params: dict,
    met_age_dict: dict,
    met_tol: float,
    age_tol: float,
) -> np.ndarray:
    """Array with the values of all the parameters for each model, in the order of
    ``FULL_PARAMS``, with the same defaults and checks as
    :py:func:`synth_model.properModel`.

    :param params: Dictionary with an array (or a single value, shared by all the
        models) for each parameter
    :type params: dict
    :param def_params: Dictionary of default parameters
    :type def_params: dict
    :param met_age_dict: Dictionary of metallicity and age values
    :type met_age_dict: dict
    :param met_tol: Tolerance for metallicity values matching
    :type met_tol: float
    :param age_tol: Tolerance for log(age) values matching
    :type age_tol: float

    :raises ValueError: If a parameter is not recognized, if the arrays have
        different lengths, if (z, a) are outside of the isochrones' ranges, if any of
        ``alpha, beta, Av, DR, dm`` is negative, or if ``Rv <= 0``

    :return: Array of shape ``(N_models, 8)``
    :rtype: np.ndarray
    """
    unknown = set(params) - set(FULL_PARAMS)
    if unknown:
        raise ValueError(f"Parameter(s) {sorted(unknown)} not recognized")

    values = {k: np.atleast_1d(np.asarray(v, dtype=float)) for k, v in params.items()}
    for k, v in values.items():
        if v.ndim != 1:
            raise ValueError(f"Parameter {k} must be a single value or a 1D array")
    lengths = {v.size for v in values.values() if v.size != 1}
    if len(lengths) > 1:
        raise ValueError("All the parameter arrays must have the same length")
    N = lengths.pop() if lengths else 1

    P = np.empty((N, len(FULL_PARAMS)))
    for i, k in enumerate(FULL_PARAMS):
        P[:, i] = values[k] if k in values else def_params[k]

    for i, (k, tol) in enumerate((("met", met_tol), ("loga", age_tol))):
        grid, v = met_age_dict[k], P[:, i]
        if len(grid) == 1:
            bad = np.abs(v - grid[0]) > tol
            if bad.any():
                raise ValueError(
                    f"Parameter {k}={v[bad][0]:.5f} != fixed value {grid[0]}"
                )
        else:
            bad = (v < grid[0]) | (v > grid[-1])
            if bad.any():
                raise ValueError(
                    f"Parameter {k}={v[bad][0]:.5f} outside isochrone"
                    + f" range [{grid[0]}, {grid[-1]}]"
                )
    for i, k in enumerate(FULL_PARAMS):
        if k in ("met", "loga", "Rv"):
            continue
        bad = P[:, i] < 0.0
        if bad.any():
            raise ValueError(f"Parameter {k}={P[:, i][bad][0]:.5f} must be >= 0")
    bad = P[:, 6] <= 0.0
    if bad.any():
        raise ValueError(f"Parameter Rv={P[:, 6][bad][0]:.5f} must be > 0")

    return P


if _HAS_NUMBA:
    from numba import njit, prange
    from numba.typed.typedlist import List

    from .likelihood_priv import _tremmel_lkl_nb
    from .synth_model import (
        _binarity_nb,
        _extinction_ccmo_nb,
        _extinction_gaiadr3_nb,
        _mass_interp_nb,
        _scatter_add_nb,
        _zaWAverage_nb,
        per_star_pool,
        per_star_pool_next,
    )

    _per_star_pool_nb = njit(cache=True)(per_star_pool)
    _per_star_pool_next_nb = njit(cache=True)(per_star_pool_next)

    @njit(cache=True)
    def _synth_one(
        p,
        met_grid,
        age_grid,
        theor_tracks,
        m_ini_idx,
        st_list,
        midx_list,
        N_a,
        dr_frac,
        ext_kind,
        ccm,
        unif1,
        max_mag,
        err,
        N_synth,
    ):
        """One synthetic cluster, with the same steps and operations as
        :py:func:`synth_model.get_model` (numba). ``ext_kind`` is ``0`` for the
        GAIADR3 extinction law and ``1`` for CCMO, with the coefficients in ``ccm``
        (see :py:func:`ccmo_array`). Returns an array of shape
        ``(m_ini_idx + 2, N)``, with ``N = 0`` for an empty cluster."""
        # p[6] (Rv) is only used by the CCMO extinction law
        met, loga, alpha, beta, av, dr, dm = p[0], p[1], p[2], p[3], p[4], p[5], p[7]
        N_rows = m_ini_idx + 2

        # properModel(): indexes of the (z, a) box
        ml = mh = 0
        if met_grid.size > 1:
            mh = min(met_grid.size - 1, np.searchsorted(met_grid, met))
            ml = max(0, mh - 1)
        al = ah = 0
        if age_grid.size > 1:
            ah = min(age_grid.size - 1, np.searchsorted(age_grid, loga))
            al = max(0, ah - 1)

        if ml == 0 and al == 0 and mh == 0 and ah == 0:
            isochrone = theor_tracks[0, 0].copy()
            m0, a0 = 0, 0
        else:
            # zaWAverage(): weights and closest grid node
            z1, z2 = met_grid[ml], met_grid[mh]
            a1, a2 = age_grid[al], age_grid[ah]
            dz = z2 - z1
            da = a2 - a1
            t = (met - z1) / dz if dz > 1e-15 else 0.0
            u = (loga - a1) / da if da > 1e-15 else 0.0
            om_t = 1.0 - t
            om_u = 1.0 - u
            m0 = mh if t >= 0.5 else ml
            a0 = ah if u >= 0.5 else al
            isochrone = _zaWAverage_nb(
                theor_tracks,
                ml,
                mh,
                al,
                ah,
                om_t * om_u,
                om_t * u,
                t * om_u,
                t * u,
                m0,
                a0,
                m_ini_idx,
            )

        binar_flag = not (alpha == 0.0 and beta == 0.0)

        # move_isochrone()
        isochrone[0] += dm
        if binar_flag:
            isochrone[m_ini_idx + 2] += dm

        # synth_model.ccmo_coefs()
        ec_mag = ec_col1 = ec_col2 = 0.0
        col2_flag = False
        if ext_kind == 1:
            rv = p[6]
            ec_mag = ccm[0, 0] + ccm[0, 1] / rv
            ec_col1 = (ccm[1, 0] + ccm[1, 1] / rv) - (ccm[1, 2] + ccm[1, 3] / rv)
            col2_flag = ccm.shape[0] > 2
            if col2_flag:
                ec_col2 = (ccm[2, 0] + ccm[2, 1] / rv) - (ccm[2, 2] + ccm[2, 3] / rv)

        # extinction() with the minimum absorption Av, on a copy of the isochrone
        isoch_av = isochrone.copy()
        if ext_kind == 0:
            _extinction_gaiadr3_nb(
                isoch_av, dr_frac, 0, av, av, 0.0, m_ini_idx, binar_flag
            )
        else:
            _extinction_ccmo_nb(
                isoch_av,
                dr_frac,
                0,
                av,
                0.0,
                m_ini_idx,
                binar_flag,
                ec_mag,
                ec_col1,
                ec_col2,
                col2_flag,
            )

        # cut_max_mag()
        msk_cut = isoch_av[0] < max_mag
        idx_cut = np.flatnonzero(msk_cut)
        if idx_cut.size == 0:
            return np.empty((isochrone.shape[0], 0))[:N_rows]

        # Pool of candidate stars, grown until N_synth stars are brighter than the
        # maximum magnitude after their extinction (see synth_model.get_model())
        k = ((ml * N_a + al) * 2 + (m0 - ml)) * 2 + (a0 - al)
        st_mass = st_list[ml * N_a + al]
        dr_kind = 1 if dr > 0.0 else 0
        N_pool = _per_star_pool_nb(N_synth, dr)
        while True:
            # mass_interp()
            isoch_mass = _mass_interp_nb(
                isochrone,
                idx_cut,
                msk_cut,
                m_ini_idx,
                st_mass,
                N_pool,
                midx_list[k],
            )
            if isoch_mass.shape[1] == 0:
                return np.empty((isochrone.shape[0], 0))[:N_rows]
            N_stars = isoch_mass.shape[1]

            # binarity()
            if not binar_flag:
                isoch_mass[m_ini_idx + 1] = np.nan
            else:
                _binarity_nb(alpha, beta, m_ini_idx, unif1, isoch_mass)
            isoch_all = isoch_mass[:N_rows]

            # extinction() of each star
            if ext_kind == 0:
                _extinction_gaiadr3_nb(
                    isoch_all, dr_frac, dr_kind, av, av, dr, m_ini_idx, False
                )
            else:
                _extinction_ccmo_nb(
                    isoch_all,
                    dr_frac,
                    dr_kind,
                    av,
                    dr,
                    m_ini_idx,
                    False,
                    ec_mag,
                    ec_col1,
                    ec_col2,
                    col2_flag,
                )

            keep = np.flatnonzero(isoch_all[0] < max_mag)
            if keep.size >= N_synth or N_stars < N_pool:
                break
            N_new = _per_star_pool_next_nb(
                N_pool, keep.size, N_stars, N_synth, st_mass.size
            )
            if N_new <= N_pool:
                break
            N_pool = N_new

        if keep.size == 0:
            return np.empty((isochrone.shape[0], 0))[:N_rows]
        isoch_binar = isoch_all[:, keep[:N_synth]]

        # add_errors()
        if err.shape[0] > 0:
            mag_sort = np.argsort(-isoch_binar[0])
            for i in range(err.shape[0]):
                _scatter_add_nb(isoch_binar[i], mag_sort, err[i])

        return isoch_binar

    @njit(cache=True)
    def _lnl_one(
        synth,
        n_phot,
        hist_args,
        nz_pos,
        row_off,
        table_flat,
        m_max,
        min_lkl,
        max_lkl,
        lkl_range,
        obs_mag_med,
        obs_col_med,
        w_mag,
        w_col,
        temperature,
    ):
        """Same as :py:meth:`likelihood_priv.Tremmel.log_lkl`, for one synthetic
        cluster"""
        if synth.shape[1] == 0:
            lkl = max_lkl - 10.0 * lkl_range
        else:
            mag = synth[0]
            lkl = _tremmel_lkl_nb(
                mag, synth[1:n_phot], hist_args, nz_pos, row_off, table_flat, m_max
            )
            if lkl <= min_lkl:
                # likelihood_priv.no_overlap_dist()
                d_mag = (obs_mag_med - np.median(mag)) / w_mag
                d_col = (obs_col_med - np.median(synth[1])) / w_col
                lkl = max_lkl - (1 + (d_mag * d_mag + d_col * d_col)) * lkl_range
        return (lkl - max_lkl) / temperature

    @njit(cache=True)
    def _generate_batch(
        P,
        met_grid,
        age_grid,
        theor_tracks,
        m_ini_idx,
        st_list,
        midx_list,
        N_a,
        dr_frac,
        ext_kind,
        ccm,
        unif1,
        max_mag,
        err,
        N_synth,
    ):
        """List with the synthetic clusters for each set of parameters in ``P``"""
        out = List()
        for j in range(P.shape[0]):
            out.append(
                _synth_one(
                    P[j],
                    met_grid,
                    age_grid,
                    theor_tracks,
                    m_ini_idx,
                    st_list,
                    midx_list,
                    N_a,
                    dr_frac,
                    ext_kind,
                    ccm,
                    unif1,
                    max_mag,
                    err,
                    N_synth,
                )
            )
        return out

    @njit(cache=True)
    def _lnl_batch_serial(
        P,
        met_grid,
        age_grid,
        theor_tracks,
        m_ini_idx,
        st_list,
        midx_list,
        N_a,
        dr_frac,
        ext_kind,
        ccm,
        unif1,
        max_mag,
        err,
        N_synth,
        n_phot,
        hist_args,
        nz_pos,
        row_off,
        table_flat,
        m_max,
        min_lkl,
        max_lkl,
        lkl_range,
        obs_mag_med,
        obs_col_med,
        w_mag,
        w_col,
        temperature,
    ):
        """Log-likelihoods of the synthetic clusters for each set of parameters in
        ``P``. The synthetic clusters are not returned"""
        out = np.empty(P.shape[0])
        for j in range(P.shape[0]):
            synth = _synth_one(
                P[j],
                met_grid,
                age_grid,
                theor_tracks,
                m_ini_idx,
                st_list,
                midx_list,
                N_a,
                dr_frac,
                ext_kind,
                ccm,
                unif1,
                max_mag,
                err,
                N_synth,
            )
            out[j] = _lnl_one(
                synth,
                n_phot,
                hist_args,
                nz_pos,
                row_off,
                table_flat,
                m_max,
                min_lkl,
                max_lkl,
                lkl_range,
                obs_mag_med,
                obs_col_med,
                w_mag,
                w_col,
                temperature,
            )
        return out

    @njit(cache=True, parallel=True)
    def _lnl_batch_parallel(
        P,
        met_grid,
        age_grid,
        theor_tracks,
        m_ini_idx,
        st_list,
        midx_list,
        N_a,
        dr_frac,
        ext_kind,
        ccm,
        unif1,
        max_mag,
        err,
        N_synth,
        n_phot,
        hist_args,
        nz_pos,
        row_off,
        table_flat,
        m_max,
        min_lkl,
        max_lkl,
        lkl_range,
        obs_mag_med,
        obs_col_med,
        w_mag,
        w_col,
        temperature,
    ):
        """Same as :py:func:`_lnl_batch_serial`, with the models distributed among
        numba's threads. Each model only reads the shared arrays and writes its own
        output, so the results are identical"""
        out = np.empty(P.shape[0])
        for j in prange(P.shape[0]):
            synth = _synth_one(
                P[j],
                met_grid,
                age_grid,
                theor_tracks,
                m_ini_idx,
                st_list,
                midx_list,
                N_a,
                dr_frac,
                ext_kind,
                ccm,
                unif1,
                max_mag,
                err,
                N_synth,
            )
            out[j] = _lnl_one(
                synth,
                n_phot,
                hist_args,
                nz_pos,
                row_off,
                table_flat,
                m_max,
                min_lkl,
                max_lkl,
                lkl_range,
                obs_mag_med,
                obs_col_med,
                w_mag,
                w_col,
                temperature,
            )
        return out

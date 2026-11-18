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

:py:class:`ModelEngine` is the single entry point used by
:py:meth:`asteca.Synthetic.generate` (with ``batch=True``) and by
:py:class:`asteca.Sampler`: it uses the compiled functions if they are available
(numba installed, see :py:func:`available`), and falls back to a loop over
:py:func:`synth_model.get_model` otherwise.
"""

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .model_data import ModelContext, ObservedData

import numpy as np

from . import HAS_NUMBA as _HAS_NUMBA
from .synth_model import get_model


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


def batch_data(ctx: "ModelContext") -> dict:
    """Arrays used by the compiled functions that do not change after the
    :py:class:`asteca.Synthetic` object is created. The typed lists hold references
    to the existing IMF samples and mass indexes (no copies). Obtained once and
    cached by :py:attr:`ModelContext.batch_data <model_data.ModelContext.batch_data>`.

    :param ctx: Data that defines the synthetic clusters
    :type ctx: ModelContext

    :return: Dictionary with the data
    :rtype: dict
    """
    from numba.typed.typedlist import List

    N_z, N_a = ctx.theor_tracks.shape[:2]
    st_list = List()
    for ml in range(N_z):
        for al in range(N_a):
            st_list.append(ctx.st_dist_mass[ml][al])
    # Indexes for the (ml, al) box and its closest node (m0, a0) are stored at
    # position ((ml * N_a + al) * 2 + (m0 - ml)) * 2 + (a0 - al). Nodes beyond the
    # grid's edges are never used (repeated keys keep the positions consistent)
    midx_list = List()
    for ml in range(N_z):
        for al in range(N_a):
            for dz in (0, 1):
                for da in (0, 1):
                    key = (ml, al, min(ml + dz, N_z - 1), min(al + da, N_a - 1))
                    midx_list.append(ctx.mass_idx[key])

    return {
        "met_grid": np.asarray(ctx.met_age_dict["met"], dtype=float),
        "age_grid": np.asarray(ctx.met_age_dict["loga"], dtype=float),
        "N_a": N_a,
        "st_list": st_list,
        "midx_list": midx_list,
        "dr_frac": {
            k: np.ascontiguousarray(v) for k, v in ctx.rand_floats["dr_frac"].items()
        },
        "ext_kind": 0 if ctx.ext_law == "GAIADR3" else 1,
        "ccm": ccmo_array(ctx.ext_law, ctx.ext_coefs),
        "unif1": np.ascontiguousarray(ctx.rand_floats["unif"][1]),
    }


def batch_args(
    ctx: "ModelContext",
    dr_model: str,
    calib: tuple[float, list[np.ndarray], int],
    P: np.ndarray,
) -> tuple:
    """Arguments of the compiled batch functions :py:func:`_generate_batch`,
    :py:func:`_lnl_batch_serial` and :py:func:`_lnl_batch_parallel` (the
    likelihood arguments are added by :py:meth:`ModelEngine.log_lkl`).

    :param ctx: Data that defines the synthetic clusters
    :type ctx: ModelContext
    :param dr_model: Differential reddening model
    :type dr_model: str
    :param calib: Maximum magnitude, error distribution and number of stars of the
        synthetic clusters, as given by ``Synthetic._calib_values()``
    :type calib: tuple[float, list[np.ndarray], int]
    :param P: Array of parameters, see :py:meth:`ModelContext.params_array
        <model_data.ModelContext.params_array>`
    :type P: np.ndarray

    :return: Tuple of arguments
    :rtype: tuple
    """
    bd = ctx.batch_data

    # Same values as in Synthetic.generate()
    max_mag_syn, err_dist_synth, N_synth_stars = calib
    if len(err_dist_synth):
        err = np.array(err_dist_synth, dtype=float)
    else:
        err = np.empty((0, 0))

    return (
        P,
        bd["met_grid"],
        bd["age_grid"],
        ctx.theor_tracks,
        ctx.m_ini_idx,
        bd["st_list"],
        bd["midx_list"],
        bd["N_a"],
        bd["dr_frac"][dr_model],
        bd["ext_kind"],
        bd["ccm"],
        bd["unif1"],
        float(max_mag_syn),
        err,
        int(N_synth_stars),
    )


class ModelEngine:
    """Generates batches of synthetic clusters and evaluates their log-likelihoods,
    with the compiled functions of this module if they are available (see
    :py:func:`available`), else with a loop over :py:func:`synth_model.get_model`.
    The results are identical either way.

    :param ctx: Data that defines the synthetic clusters
    :type ctx: ModelContext
    :param parallel: Evaluate the log-likelihoods in parallel threads (only used if
        the compiled functions are available)
    :type parallel: bool
    :param n_threads: Number of threads used if ``parallel=True``. If ``None``,
        numba's current setting is used
    :type n_threads: int | None

    :raises ValueError: If the threads are used and ``n_threads`` is not an integer
        between ``1`` and ``NUMBA_NUM_THREADS``
    """

    def __init__(  # noqa: DOC502, RUF100
        self,
        ctx: "ModelContext",
        parallel: bool = False,
        n_threads: int | None = None,
    ) -> None:
        self.ctx = ctx
        self.compiled = available(ctx.ext_law)
        self.parallel = parallel
        self.n_threads = n_threads
        # 'n_threads' is only validated if it is used
        if self.compiled and parallel and n_threads is not None:
            check_n_threads(n_threads)

    def generate(
        self,
        P: np.ndarray,
        dr_model: str,
        calib: tuple[float, list[np.ndarray], int],
    ) -> list[np.ndarray]:
        """Synthetic clusters of a batch of models, identical to those generated
        one at a time by :py:func:`synth_model.get_model`

        :param P: Array of parameters, see :py:meth:`ModelContext.params_array
            <model_data.ModelContext.params_array>`
        :type P: np.ndarray
        :param dr_model: Differential reddening model
        :type dr_model: str
        :param calib: Maximum magnitude, error distribution and number of stars of
            the synthetic clusters, as given by ``Synthetic._calib_values()``
        :type calib: tuple[float, list[np.ndarray], int]

        :return: One synthetic cluster per model (an empty 1D array if the cluster
            is empty)
        :rtype: list[np.ndarray]
        """
        if not self.compiled:
            return [get_model(self.ctx, p, dr_model, calib)[1] for p in P]
        out = _generate_batch(*batch_args(self.ctx, dr_model, calib, P))
        # An empty synthetic cluster is returned as an empty 1D array, as with
        # get_model()
        return [a if a.shape[1] else np.array([]) for a in out]

    def log_lkl(
        self,
        P: np.ndarray,
        dr_model: str,
        obs: "ObservedData",
        temperature: float,
    ) -> np.ndarray:
        """Log-likelihoods of a batch of models, identical to
        ``obs.lkl.log_lkl(get_model(ctx, p, dr_model, obs.calib_values())[1],
        temperature)`` for each row ``p`` of ``P``. With the compiled functions,
        each synthetic cluster is generated and evaluated in the same call (the
        synthetic clusters never go back to Python, which is where the gain comes
        from), optionally distributing the models among parallel threads.

        :param P: Array of parameters, see :py:meth:`ModelContext.params_array
            <model_data.ModelContext.params_array>`
        :type P: np.ndarray
        :param dr_model: Differential reddening model
        :type dr_model: str
        :param obs: Data of the observed cluster
        :type obs: ObservedData
        :param temperature: Likelihood temperature, assumed positive (validated by
            :py:class:`asteca.Sampler`)
        :type temperature: float

        :return: Array with the log-likelihood of each model, ``0`` is the best value
        :rtype: np.ndarray
        """
        lk = obs.lkl
        calib = obs.calib_values()
        if not self.compiled:
            return np.array(
                [
                    lk.log_lkl(get_model(self.ctx, p, dr_model, calib)[1], temperature)
                    for p in P
                ]
            )

        args: tuple = (
            *batch_args(self.ctx, dr_model, calib, P),
            lk.m_ini_idx,
            lk.hist_args,
            lk.nz_pos,
            lk.row_off,
            lk.table_flat,
            lk.m_max,
            lk.min_lkl,
            lk.max_lkl,
            lk.lkl_range,
            float(lk.obs_mag_median),
            float(lk.obs_col_median),
            float(lk.cmd_widths[0]),
            float(lk.cmd_widths[1]),
            float(temperature),
        )
        if not self.parallel:
            return _lnl_batch_serial(*args)
        if self.n_threads is None:
            return _lnl_batch_parallel(*args)

        import numba

        # Set the number of threads for this call only
        prev_threads = numba.get_num_threads()
        numba.set_num_threads(int(self.n_threads))
        try:
            return _lnl_batch_parallel(*args)
        finally:
            numba.set_num_threads(prev_threads)


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

        # synth_model.grid_box(): indexes of the (z, a) box
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

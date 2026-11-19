"""Compiled generation of synthetic clusters (and of their log-likelihoods).

Once a :py:class:`asteca.Synthetic` object is created, every random value used by
:py:func:`synth_model.get_model` is drawn in advance, so a model is a deterministic
function of its parameters. :py:func:`_synth_one` generates a synthetic cluster in a
single compiled call, chaining the numba kernels of :py:mod:`synth_model` in the same
order and with the same operations as :py:func:`synth_model.get_model` (the Python
implementation, used without numba), so that the results are identical. The batch
functions loop over many parameter sets in a single call, which removes the Python
overhead of each call (which dominates for clusters with few stars), and the batch of
log-likelihoods can be evaluated in parallel threads.

:py:class:`model_data.ModelEngine` is the single entry point used to generate
synthetic clusters (one or many) by :py:class:`asteca.Synthetic` and
:py:class:`asteca.Sampler`: it uses the compiled functions of this module if numba is
installed, and :py:func:`synth_model.get_model` otherwise.
"""

from collections import namedtuple
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from .model_data import ModelEngine

import numpy as np

from . import HAS_NUMBA as _HAS_NUMBA

# Output of _synth_one(): the synthetic cluster and the data used by
# ModelEngine.get_model() to obtain the isochrone, the mass information and the reason
# for an empty cluster
SynthOut = namedtuple(
    "SynthOut",
    ("synth", "status", "isoch_moved", "isoch_av", "idx_cut", "N_pool", "keep"),
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


def batch_data(ctx: "ModelEngine") -> dict:
    """Arrays used by the compiled functions that do not change after the
    :py:class:`asteca.Synthetic` object is created. The typed lists hold references
    to the existing IMF samples and mass indexes (no copies). Obtained once and
    cached by :py:attr:`ModelEngine.batch_data <model_data.ModelEngine.batch_data>`.

    :param ctx: Data that defines the synthetic clusters
    :type ctx: ModelEngine

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


def model_args(
    ctx: "ModelEngine",
    dr_model: str,
    calib: tuple[float, list[np.ndarray], int],
) -> tuple[tuple, Any, Any]:
    """Data of the synthetic clusters passed to the compiled functions
    :py:func:`_synth_one`, :py:func:`_generate_batch`, :py:func:`_lnl_batch_serial`
    and :py:func:`_lnl_batch_parallel`: the ``model`` tuple, and the typed lists of
    the IMF samples and of their mass indexes. The typed lists are separate
    arguments because numba types a tuple that contains them on every call (~0.1 ms)

    :param ctx: Data that defines the synthetic clusters
    :type ctx: ModelEngine
    :param dr_model: Differential reddening model
    :type dr_model: str
    :param calib: Maximum magnitude, error distribution and number of stars of the
        synthetic clusters, as given by ``Synthetic._calib_values()``
    :type calib: tuple[float, list[np.ndarray], int]

    :return: ``(model, st_list, midx_list)``, with ``model = (met_grid, age_grid,
        theor_tracks, m_ini_idx, N_a, dr_frac, ext_kind, ccm, unif1, max_mag, err,
        N_synth)``
    :rtype: tuple[tuple, Any, Any]
    """
    bd = ctx.batch_data

    # Same values as in Synthetic.generate()
    max_mag_syn, err_dist_synth, N_synth_stars = calib
    if len(err_dist_synth):
        err = np.array(err_dist_synth, dtype=float)
    else:
        err = np.empty((0, 0))

    model = (
        bd["met_grid"],
        bd["age_grid"],
        ctx.theor_tracks,
        ctx.m_ini_idx,
        bd["N_a"],
        bd["dr_frac"][dr_model],
        bd["ext_kind"],
        bd["ccm"],
        bd["unif1"],
        float(max_mag_syn),
        err,
        int(N_synth_stars),
    )
    return model, bd["st_list"], bd["midx_list"]


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
    def _synth_one(p, model, st_list, midx_list):
        """One synthetic cluster, with the same steps and operations as
        :py:func:`synth_model.get_model` (numba). ``model, st_list, midx_list`` are
        returned by :py:func:`model_args`; the ``ext_kind`` of ``model`` is ``0`` for
        the GAIADR3 extinction law and ``1`` for CCMO, with the coefficients in
        ``ccm`` (see :py:func:`ccmo_array`).

        Returns a ``SynthOut`` named tuple ``(synth, status, isoch_moved, isoch_av,
        idx_cut, N_pool, keep)``: the synthetic cluster, an array of shape ``(m_ini_idx + 2, N)`` with
        ``N = 0`` for an empty cluster; the reason for an empty cluster (see
        ``synth_model.EMPTY_REASONS``, ``0`` if not empty); and the data used by
        :py:meth:`model_data.ModelEngine.get_model` to obtain the isochrone and the mass
        information: the isochrone moved by the distance modulus, the isochrone
        with the minimum absorption ``Av``, the indexes of its points brighter than
        the maximum magnitude, the final size of the pool of candidate stars, and
        the indexes of the candidates brighter than the maximum magnitude."""
        (
            met_grid,
            age_grid,
            theor_tracks,
            m_ini_idx,
            N_a,
            dr_frac,
            ext_kind,
            ccm,
            unif1,
            max_mag,
            err,
            N_synth,
        ) = model
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
        empty = np.empty((isochrone.shape[0], 0))[:N_rows]
        keep = np.empty(0, dtype=np.intp)
        if idx_cut.size == 0:
            return SynthOut(empty, 1, isochrone, isoch_av, idx_cut, 0, keep)

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
                return SynthOut(empty, 2, isochrone, isoch_av, idx_cut, N_pool, keep)
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
            return SynthOut(empty, 3, isochrone, isoch_av, idx_cut, N_pool, keep)
        isoch_binar = isoch_all[:, keep[:N_synth]]

        # add_errors()
        if err.shape[0] > 0:
            mag_sort = np.argsort(-isoch_binar[0])
            for i in range(err.shape[0]):
                _scatter_add_nb(isoch_binar[i], mag_sort, err[i])

        return SynthOut(isoch_binar, 0, isochrone, isoch_av, idx_cut, N_pool, keep)

    @njit(cache=True)
    def _lnl_one(synth, lkl, temperature):
        """Same as :py:meth:`likelihood_priv.Tremmel.log_lkl`, for one synthetic
        cluster. ``lkl`` is :py:attr:`Tremmel.kernel_args
        <likelihood_priv.Tremmel.kernel_args>`"""
        (
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
        ) = lkl
        if synth.shape[1] == 0:
            lnl = max_lkl - 10.0 * lkl_range
        else:
            mag = synth[0]
            lnl = _tremmel_lkl_nb(
                mag, synth[1:n_phot], hist_args, nz_pos, row_off, table_flat, m_max
            )
            if lnl <= min_lkl:
                # likelihood_priv.no_overlap_dist()
                d_mag = (obs_mag_med - np.median(mag)) / w_mag
                d_col = (obs_col_med - np.median(synth[1])) / w_col
                lnl = max_lkl - (1 + (d_mag * d_mag + d_col * d_col)) * lkl_range
        return (lnl - max_lkl) / temperature

    @njit(cache=True)
    def _generate_batch(P, model, st_list, midx_list):
        """List with the synthetic clusters for each set of parameters in ``P``"""
        out = List()
        for j in range(P.shape[0]):
            out.append(_synth_one(P[j], model, st_list, midx_list).synth)
        return out

    @njit(cache=True)
    def _lnl_batch_serial(P, model, st_list, midx_list, lkl, temperature):
        """Log-likelihoods of the synthetic clusters for each set of parameters in
        ``P``. The synthetic clusters are not returned"""
        out = np.empty(P.shape[0])
        for j in range(P.shape[0]):
            synth = _synth_one(P[j], model, st_list, midx_list).synth
            out[j] = _lnl_one(synth, lkl, temperature)
        return out

    # A separate function, not _lnl_batch_serial() compiled with parallel=True:
    # numba's cache does not tell the two versions of the same function apart
    @njit(cache=True, parallel=True)
    def _lnl_batch_parallel(P, model, st_list, midx_list, lkl, temperature):
        """Same as :py:func:`_lnl_batch_serial`, with the models distributed among
        numba's threads. Each model only reads the shared arrays and writes its own
        output, so the results are identical"""
        out = np.empty(P.shape[0])
        for j in prange(P.shape[0]):
            synth = _synth_one(P[j], model, st_list, midx_list).synth
            out[j] = _lnl_one(synth, lkl, temperature)
        return out

"""Compiled generation of synthetic clusters (and of their log-likelihoods).

Once a :py:class:`asteca.Synthetic` object is created, every random value used by
:py:func:`generate.synth_one` is drawn in advance, so a model is a deterministic
function of its parameters. :py:func:`synth_one_nb` generates a synthetic cluster in a
single compiled call, chaining the numba kernels of :py:mod:`generate` in the same
order and with the same operations as :py:func:`generate.synth_one` (the Python
implementation, used without numba), so that the results are identical. The batch
functions loop over many parameter sets in a single call, which removes the Python
overhead of each call (which dominates for clusters with few stars), and the batch of
log-likelihoods can be evaluated in parallel threads.

:py:class:`engine.SynthEngine` is the single entry point used to generate
synthetic clusters (one or many) by :py:class:`asteca.Synthetic` and
:py:class:`asteca.Sampler`: it uses the compiled functions of this module if numba is
installed, and :py:func:`generate.synth_one` otherwise.
"""

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .engine import SynthEngine

import numpy as np

from .. import HAS_NUMBA as _HAS_NUMBA
from .generate import SynthOut


def kernel_data(ctx: "SynthEngine") -> dict:
    """Arrays used by the compiled functions, obtained once for each engine and
    cached by :py:attr:`SynthEngine.kernel_data <engine.SynthEngine.kernel_data>`.

    :param ctx: Data that defines the synthetic clusters
    :type ctx: SynthEngine

    :return: Dictionary with the data
    :rtype: dict
    """
    return {
        "met_grid": np.asarray(ctx.met_age_dict["met"], dtype=float),
        "age_grid": np.asarray(ctx.met_age_dict["loga"], dtype=float),
        "N_a": ctx.theor_tracks.shape[1],
        "imf": (
            ctx.imf_cdf,
            ctx.imf_lut,
            ctx.rand_floats["imf_u"],
            np.ascontiguousarray(ctx.rand_floats["imf_offsets"], dtype=np.int64),
        ),
        "dr_frac": {
            k: np.ascontiguousarray(v) for k, v in ctx.rand_floats["dr_frac"].items()
        },
        "ext_kind": 0 if ctx.ext_law == "GAIADR3" else 1,
        # Rows [a, b, 0, 0] (magnitude) and [a1, b1, a2, b2] (each color), see
        # precompute.ccmo_filter_coefs(); no rows for GAIADR3
        "ccm": np.ascontiguousarray(ctx.ext_coefs, dtype=float),
        "unif1": np.ascontiguousarray(ctx.rand_floats["unif"][1]),
        "err": (
            np.ascontiguousarray(ctx.err_dist, dtype=float)
            if len(ctx.err_dist)
            else np.empty((0, 0))
        ),
    }


def synth_args(ctx: "SynthEngine", dr_model: str) -> tuple[tuple, tuple]:
    """Data of the synthetic clusters passed to the compiled functions
    :py:func:`synth_one_nb`, :py:func:`synth_batch_nb`, :py:func:`log_lkl_batch_nb`
    and :py:func:`log_lkl_batch_parallel_nb`: the ``model`` tuple, and the ``imf`` tuple
    with the data used to draw the masses of the candidate stars: the IMF CDF of each
    grid node and their lookup tables (see :py:func:`precompute.imf_cdf_tables`), and
    the uniform values and the offset of the window of each IMF node (see
    :py:func:`precompute.draw_rand_floats`)

    :param ctx: Data that defines the synthetic clusters
    :type ctx: SynthEngine
    :param dr_model: Differential reddening model
    :type dr_model: str

    :return: ``(model, imf)``, with ``model = (met_grid, age_grid, theor_tracks,
        m_ini_idx, N_a, dr_frac, ext_kind, ccm, unif1, max_mag, err, N_synth)`` and
        ``imf = (imf_cdf, imf_lut, imf_u, imf_offsets)``
    :rtype: tuple[tuple, tuple]
    """
    bd = ctx.kernel_data
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
        float(ctx.max_mag),
        bd["err"],
        int(ctx.N_synth),
    )
    return model, bd["imf"]


if _HAS_NUMBA:
    from numba import njit, prange
    from numba.typed.typedlist import List

    from .generate import (
        _apply_extinction_ccmo_nb,
        _apply_extinction_gaiadr3_nb,
        _assign_binaries_nb,
        _draw_candidates_nb,
        _isochrone_at_met_age_nb,
        _scatter_add_nb,
        initial_pool_size,
        next_pool_size,
    )
    from .likelihood import _tremmel_sum_nb

    _initial_pool_size_nb = njit(cache=True)(initial_pool_size)
    _next_pool_size_nb = njit(cache=True)(next_pool_size)

    @njit(cache=True)
    def synth_one_nb(p, model, imf):
        """One synthetic cluster, with the same steps and operations as
        :py:func:`generate.synth_one` (numba). ``model, imf`` are
        returned by :py:func:`synth_args`; the ``ext_kind`` of ``model`` is ``0`` for
        the GAIADR3 extinction law and ``1`` for CCMO, with the coefficients in
        ``ccm`` (see :py:func:`precompute.ccmo_filter_coefs`).

        Returns a ``generate.SynthOut`` named tuple: the synthetic cluster (an array
        of shape ``(m_ini_idx + 2, N)`` with ``N = 0`` for an empty cluster), the
        reason for an empty cluster (see ``generate.EMPTY_REASONS``, ``0`` if not
        empty), and the data used by :py:meth:`engine.SynthEngine.synth_one` to
        obtain the isochrone and the mass information."""
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

        # generate.grid_box(): indexes of the (z, a) box
        ml = mh = 0
        if met_grid.size > 1:
            mh = min(met_grid.size - 1, np.searchsorted(met_grid, met))
            ml = max(0, mh - 1)
        al = ah = 0
        if age_grid.size > 1:
            ah = min(age_grid.size - 1, np.searchsorted(age_grid, loga))
            al = max(0, ah - 1)

        # IMF node of the model (see generate.imf_uniforms())
        node = ml * N_a + al

        if ml == 0 and al == 0 and mh == 0 and ah == 0:
            isochrone = theor_tracks[0, 0].copy()
            m0, a0 = 0, 0
        else:
            # isochrone_at_met_age(): weights and closest grid node
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
            isochrone = _isochrone_at_met_age_nb(
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

        # apply_distance()
        isochrone[0] += dm
        if binar_flag:
            isochrone[m_ini_idx + 2] += dm

        # generate.ccmo_ratios()
        ec_mag = ec_col1 = ec_col2 = 0.0
        col2_flag = False
        if ext_kind == 1:
            rv = p[6]
            ec_mag = ccm[0, 0] + ccm[0, 1] / rv
            ec_col1 = (ccm[1, 0] + ccm[1, 1] / rv) - (ccm[1, 2] + ccm[1, 3] / rv)
            col2_flag = ccm.shape[0] > 2
            if col2_flag:
                ec_col2 = (ccm[2, 0] + ccm[2, 1] / rv) - (ccm[2, 2] + ccm[2, 3] / rv)

        # apply_extinction() with the minimum absorption Av, on a copy of the isochrone
        isoch_av = isochrone.copy()
        if ext_kind == 0:
            _apply_extinction_gaiadr3_nb(
                isoch_av, dr_frac, 0, av, av, 0.0, m_ini_idx, binar_flag
            )
        else:
            _apply_extinction_ccmo_nb(
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

        # bright_points()
        msk_cut = isoch_av[0] < max_mag
        idx_cut = np.flatnonzero(msk_cut)
        empty = np.empty((isochrone.shape[0], 0))[:N_rows]
        keep = np.empty(0, dtype=np.intp)
        if idx_cut.size == 0:
            return SynthOut(
                empty, 1, isochrone, isoch_av, idx_cut, 0, keep, m0, a0, node
            )

        # Pool of candidate stars, grown until N_synth stars are brighter than the
        # maximum magnitude after their extinction (see generate.synth_one())
        # Uniform values of the candidates of the IMF node (see
        # generate.imf_uniforms()), and the IMF CDF of the grid node (m0, a0)
        imf_cdf, imf_lut, imf_u_all, imf_offsets = imf
        off = imf_offsets[node]
        imf_u = imf_u_all[off : off + imf_u_all.size // 2]
        dr_kind = 1 if dr > 0.0 else 0
        N_pool = _initial_pool_size_nb(N_synth, dr)
        while True:
            # draw_candidates()
            isoch_mass = _draw_candidates_nb(
                isochrone,
                idx_cut,
                msk_cut,
                m_ini_idx,
                imf_u,
                N_pool,
                imf_cdf[m0, a0],
                imf_lut[m0, a0],
            )
            if isoch_mass.shape[1] == 0:
                return SynthOut(
                    empty, 2, isochrone, isoch_av, idx_cut, N_pool, keep, m0, a0, node
                )
            N_stars = isoch_mass.shape[1]

            # assign_binaries()
            if not binar_flag:
                isoch_mass[m_ini_idx + 1] = np.nan
            else:
                _assign_binaries_nb(alpha, beta, m_ini_idx, unif1, isoch_mass)
            isoch_all = isoch_mass[:N_rows]

            # apply_extinction() of each star
            if ext_kind == 0:
                _apply_extinction_gaiadr3_nb(
                    isoch_all, dr_frac, dr_kind, av, av, dr, m_ini_idx, False
                )
            else:
                _apply_extinction_ccmo_nb(
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
            N_new = _next_pool_size_nb(N_pool, keep.size, N_stars, N_synth, imf_u.size)
            if N_new <= N_pool:
                break
            N_pool = N_new

        if keep.size == 0:
            return SynthOut(
                empty, 3, isochrone, isoch_av, idx_cut, N_pool, keep, m0, a0, node
            )
        isoch_binar = isoch_all[:, keep[:N_synth]]

        # add_errors()
        if err.shape[0] > 0:
            mag_sort = np.argsort(-isoch_binar[0])
            for i in range(err.shape[0]):
                _scatter_add_nb(isoch_binar[i], mag_sort, err[i])

        return SynthOut(
            isoch_binar, 0, isochrone, isoch_av, idx_cut, N_pool, keep, m0, a0, node
        )

    @njit(cache=True)
    def log_lkl_one_nb(synth, lkl, temperature):
        """Same as :py:meth:`likelihood.Tremmel.log_lkl`, for one synthetic
        cluster. ``lkl`` is :py:attr:`Tremmel.kernel_args
        <likelihood.Tremmel.kernel_args>`"""
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
            lnl = _tremmel_sum_nb(
                mag, synth[1:n_phot], hist_args, nz_pos, row_off, table_flat, m_max
            )
            if lnl <= min_lkl:
                # likelihood.no_overlap_dist()
                d_mag = (obs_mag_med - np.median(mag)) / w_mag
                d_col = (obs_col_med - np.median(synth[1])) / w_col
                lnl = max_lkl - (1 + (d_mag * d_mag + d_col * d_col)) * lkl_range
        return (lnl - max_lkl) / temperature

    @njit(cache=True)
    def synth_batch_nb(P, model, imf):
        """List with the synthetic clusters for each set of parameters in ``P``"""
        out = List()
        for j in range(P.shape[0]):
            out.append(synth_one_nb(P[j], model, imf).synth)
        return out

    @njit(cache=True)
    def log_lkl_batch_nb(P, model, imf, lkl, temperature):
        """Log-likelihoods of the synthetic clusters for each set of parameters in
        ``P``. The synthetic clusters are not returned"""
        out = np.empty(P.shape[0])
        for j in range(P.shape[0]):
            synth = synth_one_nb(P[j], model, imf).synth
            out[j] = log_lkl_one_nb(synth, lkl, temperature)
        return out

    # A separate function, not log_lkl_batch_nb() compiled with parallel=True:
    # numba's cache does not tell the two versions of the same function apart
    @njit(cache=True, parallel=True)
    def log_lkl_batch_parallel_nb(P, model, imf, lkl, temperature):
        """Same as :py:func:`log_lkl_batch_nb`, with the models distributed among
        numba's threads. Each model only reads the shared arrays and writes its own
        output, so the results are identical"""
        out = np.empty(P.shape[0])
        for j in prange(P.shape[0]):
            synth = synth_one_nb(P[j], model, imf).synth
            out[j] = log_lkl_one_nb(synth, lkl, temperature)
        return out

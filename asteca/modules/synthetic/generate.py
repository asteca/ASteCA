from bisect import bisect_left
from collections import namedtuple
from typing import TYPE_CHECKING

import numpy as np

from .. import HAS_NUMBA as _HAS_NUMBA
from .. import njit

if TYPE_CHECKING:
    from .engine import SynthEngine

# The 'coeffs' values are the main sequence values taken from:
# https://www.cosmos.esa.int/web/gaia/edr3-extinction-law
_G_COEFFS = (
    0.995969721536602,
    -0.159726460302015,
    0.0122380738156057,
    0.00090726555099859,
    -0.0377160263914123,
    0.00151347495244888,
    -2.52364537395142e-05,
    0.0114522658102451,
    -0.000936914989014318,
    -0.000260296774134201,
)
# "BP" minus "RP" coefficients
_BPRP_COEFFS = (
    0.49042409542357,
    -0.0634165826723407,
    -0.036506793426665456,
    0.02189430261385411,
    -0.01588332677592224,
    0.0008075446901998552,
    -1.468074307777017e-05,
    0.006681041701047307,
    -0.001137927314099498,
    -0.0001222634137884233,
)
# Range of intrinsic BP-RP colors for which the coefficients above were fitted. The
# coefficients should not be extrapolated beyond it (the polynomials diverge), so
# the colors used to evaluate them are clamped to this range
_X_MIN, _X_MAX = -0.06, 2.5


# Reasons for an empty synthetic cluster, indexed by the status code returned by
# compiled.synth_one_nb() (0: not empty)
EMPTY_REASONS = ("", "max magnitude cut", "IMF interpolation", "the extinction")

# Output of compiled.synth_one_nb() (also built by synth_one() for fill_mass_info()): the
# synthetic cluster, the reason for an empty cluster (see EMPTY_REASONS), and the data
# used to obtain the isochrone and the mass information: the isochrone moved by the
# distance modulus, the isochrone with the minimum absorption Av, the indexes of its
# points brighter than the maximum magnitude, the final size of the pool of candidate
# stars, the indexes of the candidates brighter than the maximum magnitude, the grid
# node (m0, a0) the masses are taken from, and the IMF node (see imf_uniforms())
SynthOut = namedtuple(
    "SynthOut",
    (
        "synth",
        "status",
        "isoch_moved",
        "isoch_av",
        "idx_cut",
        "N_pool",
        "keep",
        "m0",
        "a0",
        "node",
    ),
)


def synth_one(
    ctx: "SynthEngine",
    p: np.ndarray,
    dr_model: str,
    show_errors: bool = False,
    isoch_flag: bool = False,
    isoch_dr_frac: float = 0.0,
    mass_info: dict | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Generate a synthetic cluster.

    The synthetic cluster is generated according to the parameters in ``p`` and the
    data of the :py:class:`Synthetic` object. This is the Python implementation of
    the synthetic cluster recipe, used when numba is not installed; with numba, the
    compiled ``compiled.synth_one_nb()`` is used instead (through
    ``engine.SynthEngine.synth_one()``), with identical results.

    :param ctx: Data that defines the synthetic clusters, including their maximum
        magnitude, error distribution and number of stars
    :type ctx: SynthEngine
    :param p: Values of the fundamental parameters, in the order of ``FULL_PARAMS``:
        a row of the array returned by :py:meth:`SynthEngine.validate_params
        <engine.SynthEngine.validate_params>` (the values are not validated here)
    :type p: np.ndarray
    :param dr_model: Differential reddening model, it selects the fractions of the
        differential reddening applied to each star (see
        :py:func:`precompute.dr_fractions`)
    :type dr_model: str
    :param show_errors: If ``True``, print the reason when an empty synthetic cluster
        is returned
    :type show_errors: bool
    :param isoch_flag: ``True`` if the isochrone should be returned along with the
        synthetic cluster, ``False`` otherwise.
    :type isoch_flag: bool
    :param isoch_dr_frac: Fraction ``f`` of the differential reddening applied to
        the returned isochrone, i.e.: it is reddened with ``Av + f * DR`` instead of
        ``Av``. It does not affect the synthetic cluster
    :type isoch_dr_frac: float
    :param mass_info: If a dictionary is given, it is filled with the data required
        to estimate the total mass of the simulated cluster, see
        :py:func:`fill_mass_info`. It is not modified if the synthetic cluster is
        empty
    :type mass_info: dict | None

    :return: Synthetic isochrone plus synthetic array; if ``isoch_flag=False``,
        an empty array is returned instead of the isochrone.
        The synthetic array contains the data ``[mag, c1, (c2), mass, mass_b]``,
        where ``mag`` is the magnitude, ``c1`` is the color, ``c2`` is the
        optional second color, and ``mass, mass_b`` are the masses of the single
        and secondary components of the binary systems, respectively
        (if generated). If the system is a single star, then ``mass_b==np.nan``.
    :rtype: tuple[np.ndarray, np.ndarray]
    """

    met_age_dict = ctx.met_age_dict
    m_ini_idx, theor_tracks = ctx.m_ini_idx, ctx.theor_tracks
    ext_law, ext_coefs, rand_floats = ctx.ext_law, ctx.ext_coefs, ctx.rand_floats
    dr_frac = rand_floats["dr_frac"][dr_model]
    max_mag_syn, err_dist_synth, N_synth_stars = ctx.max_mag, ctx.err_dist, ctx.N_synth

    met, loga, alpha, beta, av, dr, rv, dm = p.tolist()
    # Indexes of the (z, log(age)) box used for the isochrone averaging
    ml, mh, al, ah = grid_box(met_age_dict, met, loga)
    # If (z, a) are both fixed, use the single processed isochrone
    if ml == al == mh == ah == 0:
        # The np.array() is important to avoid overwriting 'theor_tracks'
        isochrone = np.array(theor_tracks[0][0])
        # Grid node that the masses are taken from
        m0, a0 = 0, 0
    else:
        # Generate a weighted average isochrone from the (z, log(age)) values in
        # the 'model'. (m0, a0) is the closest grid node, from which the masses
        # are taken
        isochrone, m0, a0 = isochrone_at_met_age(
            theor_tracks,
            met_age_dict,
            m_ini_idx,
            met,
            loga,
            ml,
            mh,
            al,
            ah,
        )

    binar_flag = not (alpha == 0.0 and beta == 0.0)

    # Move theoretical isochrone using the distance modulus
    isoch_moved = apply_distance(isochrone, dm)

    # Isochrone with the minimum absorption Av (no differential reddening). It is
    # only used to identify the isochrone points that can produce stars brighter
    # than the maximum magnitude (a larger absorption only makes a star fainter),
    # and returned as the isochrone if requested. The copy avoids reddening
    # 'isoch_moved', which is interpolated below
    isoch_av = apply_extinction(
        ext_law,
        ext_coefs,
        dr_frac,
        av,
        0.0,
        rv,
        isoch_moved.copy(),
    )
    isoch_extin_dr_0 = isochrone_at_dr_frac(
        ctx, p, dr_model, isoch_flag, isoch_dr_frac, isoch_moved, isoch_av
    )

    # Identify the isochrone points below the maximum magnitude. The isochrone
    # itself is not cut (which would copy it): the points are selected through
    # their indexes in draw_candidates()
    idx_cut, msk_cut = bright_points(isoch_av, max_mag_syn)
    if idx_cut.size == 0:
        if show_errors:
            print(empty_message(1, p))
        return isoch_extin_dr_0, np.array([])

    # Uniform values of the candidate stars of this IMF node, and the IMF CDF at the
    # masses of the isochrone (those of the grid node (m0, a0))
    node = ml * theor_tracks.shape[1] + al
    imf_u = imf_uniforms(rand_floats, node)
    imf_F = ctx.imf_cdf[m0, a0]

    # The extinction is applied to each synthetic star (after the interpolation and
    # the binarity), so that each star gets its own differential reddening draw
    # (see initial_pool_size()). As some stars end up below the maximum magnitude, a
    # pool of candidates is generated, with masses drawn from the IMF within the
    # isochrone's mass range (see draw_candidates()), sized with the expected
    # fraction of lost stars (see lost_fraction(), for clusters of at least
    # POOL_ESTIMATE_MIN_STARS stars). The pool is grown until 'N_synth_stars' stars are brighter than the maximum magnitude, or all the
    # candidates are used. The candidates (and their random values) are always taken
    # in the same order, so the result does not depend on the pool sizes
    q_lost = 0.0
    if dr > 0.0 and N_synth_stars < POOL_ESTIMATE_MIN_STARS:
        q_lost = Q_LOST_DEFAULT
    elif dr > 0.0:
        k_av0 = k_max = GAIADR3_K_MAX
        if ext_law == "CCMO":
            k_av0 = ccmo_ratios(ext_coefs, rv)[0]
            k_max = 1.01 * k_av0
        q_lost = lost_fraction(
            isoch_av[0],
            isoch_moved[0],
            av,
            dr,
            k_av0,
            k_max,
            imf_F,
            int(idx_cut[0]),
            int(idx_cut[-1]),
            max_mag_syn,
            rand_floats["dr_frac_q"][dr_model],
        )
    N_pool = initial_pool_size(N_synth_stars, q_lost)
    while True:
        # Draw the masses of the candidates and interpolate them into the isochrone
        isoch_mass = draw_candidates(
            isoch_moved,
            idx_cut,
            msk_cut,
            m_ini_idx,
            imf_u,
            N_pool,
            imf_F,
        )
        if isoch_mass.size == 0:
            if show_errors:
                print(empty_message(2, p))
            return isoch_extin_dr_0, np.array([])
        N_stars = isoch_mass.shape[1]

        # Assignment of binarity. The secondary of each binary system is
        # interpolated in the (full) isochrone moved by the distance modulus
        isoch_binar = assign_binaries(
            alpha,
            beta,
            binar_flag,
            m_ini_idx,
            rand_floats["unif"][1],
            ctx.binar_tables,
            isoch_moved,
            isoch_mass,
        )

        # Extinction of each star, from the intrinsic color of the system
        isoch_binar = apply_extinction(
            ext_law,
            ext_coefs,
            dr_frac,
            av,
            dr,
            rv,
            isoch_binar,
        )

        keep = np.flatnonzero(isoch_binar[0] < max_mag_syn)
        # Enough stars, or all the candidates used
        if keep.size >= N_synth_stars or N_stars < N_pool:
            break
        N_new = next_pool_size(N_pool, keep.size, N_stars, N_synth_stars, imf_u.size)
        if N_new <= N_pool:
            break
        N_pool = N_new

    if keep.size == 0:
        if show_errors:
            print(empty_message(3, p))
        return isoch_extin_dr_0, np.array([])

    if mass_info is not None:
        out = SynthOut(
            isoch_binar, 0, isoch_moved, isoch_av, idx_cut, N_pool, keep, m0, a0, node
        )
        fill_mass_info(mass_info, ctx, p, out)

    isoch_binar = isoch_binar[:, keep[:N_synth_stars]]

    # Assign errors according to errors distribution.
    synth_clust = add_errors(isoch_binar, err_dist_synth)

    return isoch_extin_dr_0, synth_clust


def empty_message(status: int, p: np.ndarray) -> str:
    """Message printed when an empty synthetic cluster is returned and
    ``show_errors=True``

    :param status: Reason for the empty cluster, see ``EMPTY_REASONS``
    :type status: int
    :param p: Values of the fundamental parameters, see :py:func:`synth_one`
    :type p: np.ndarray

    :return: Message
    :rtype: str
    """
    met, loga, _, _, av, dr, _, dm = p.tolist()
    return (
        f"No stars left after {EMPTY_REASONS[status]} (met={met:.5f},"
        + f"log(age)={loga:.3f}, Av={av:.3f}, DR={dr:.3f}, dm={dm:.3f})"
    )


def isochrone_at_dr_frac(
    ctx: "SynthEngine",
    p: np.ndarray,
    dr_model: str,
    isoch_flag: bool,
    isoch_dr_frac: float,
    isoch_moved: np.ndarray,
    isoch_av: np.ndarray,
) -> np.ndarray:
    """Isochrone returned by :py:func:`synth_one`: the isochrone with the minimum
    absorption ``Av``, or reddened with ``Av + isoch_dr_frac * DR`` if
    ``isoch_dr_frac > 0``; an empty array if ``isoch_flag=False``

    :param ctx: Data that defines the synthetic clusters
    :type ctx: SynthEngine
    :param p: Values of the fundamental parameters, see :py:func:`synth_one`
    :type p: np.ndarray
    :param dr_model: Differential reddening model
    :type dr_model: str
    :param isoch_flag: Return the isochrone
    :type isoch_flag: bool
    :param isoch_dr_frac: Fraction of the differential reddening applied
    :type isoch_dr_frac: float
    :param isoch_moved: Isochrone moved by the distance modulus (not modified)
    :type isoch_moved: np.ndarray
    :param isoch_av: Isochrone with the minimum absorption ``Av``
    :type isoch_av: np.ndarray

    :return: Isochrone
    :rtype: np.ndarray
    """
    if not isoch_flag:
        return np.array([])
    _, _, _, _, av, dr, rv, _ = p.tolist()
    if isoch_dr_frac > 0 and dr > 0:
        return apply_extinction(
            ctx.ext_law,
            ctx.ext_coefs,
            ctx.rand_floats["dr_frac"][dr_model],
            av + isoch_dr_frac * dr,
            0.0,
            rv,
            isoch_moved.copy(),
        )
    return isoch_av


def fill_mass_info(
    mass_info: dict, ctx: "SynthEngine", p: np.ndarray, out: SynthOut
) -> None:
    """Fill ``mass_info`` with the data required to estimate the total mass of a
    (non empty) simulated cluster (see :py:func:`cluster_mass.actual_mass`).

    The candidates of the synthetic cluster are drawn from the IMF restricted to the
    mass range of the isochrone points brighter than the maximum magnitude (see
    :py:func:`draw_candidates`). The simulated cluster is made of the candidates up to
    the last one kept (or all of them, if fewer than ``ctx.N_synth`` were kept),
    plus the stars of the IMF outside of that range that they represent.

    Keys: the masses of those candidates (``pool_masses``), the indexes of the
    stars kept (``kept_idx``), the number of stars requested (``N_target``), the IMF
    CDF at the limits of the mass range (``F_lo``, ``F_hi``) and at the largest
    mass of the isochrone (``F_alive``), the knots of the inverse CDF
    (``imf_knots``), the initial masses of the isochrone (``mass_row``), the
    binary fraction parameters (``alpha``, ``beta``), and the tables of the
    mass-ratio distribution (``q_dist``, see :py:func:`precompute.mass_ratio_dist`)

    :param mass_info: Dictionary to fill
    :type mass_info: dict
    :param ctx: Data that defines the synthetic clusters
    :type ctx: SynthEngine
    :param p: Values of the fundamental parameters, see :py:func:`synth_one`
    :type p: np.ndarray
    :param out: Data of the synthetic cluster, see ``SynthOut``
    :type out: SynthOut
    """
    m_ini_idx, N_synth_stars = ctx.m_ini_idx, ctx.N_synth
    alpha, beta = p[2:4].tolist()
    # IMF node and grid node of the model (see synth_one())
    imf_u = imf_uniforms(ctx.rand_floats, out.node)
    imf_F = ctx.imf_cdf[out.m0, out.a0]
    mass_row = out.isoch_moved[m_ini_idx]
    idx_cut, keep = out.idx_cut, out.keep

    kept_idx = keep[:N_synth_stars]
    N_end = int(kept_idx[-1]) + 1 if keep.size >= N_synth_stars else out.N_pool
    N_end = min(N_end, imf_u.size)
    i_lo, i_hi = idx_cut[0], idx_cut[-1]
    F_lo = float(imf_F[i_lo])
    _, pool_masses = imf_masses(
        mass_row,
        imf_F,
        F_lo,
        float(imf_F[i_hi]) - F_lo,
        mass_row[i_lo],
        mass_row[i_hi],
        imf_u[:N_end],
    )
    mass_info.update(
        {
            "pool_masses": pool_masses,
            "kept_idx": kept_idx,
            "N_target": N_synth_stars,
            "F_lo": float(F_lo),
            "F_hi": float(imf_F[i_hi]),
            "F_alive": float(imf_F[-1]),
            "imf_knots": ctx.imf_knots,
            "mass_row": mass_row.copy(),
            "alpha": alpha,
            "beta": beta,
            "q_dist": ctx.q_dist,
        }
    )


def imf_uniforms(rand_floats: dict, node: int) -> np.ndarray:
    """Uniform values of the candidate stars of the IMF node ``node`` (the index of
    the lower corner ``(ml, al)`` of the model's (z, log(age)) box, ``ml * N_a +
    al``): a window of ``imf_u`` starting at its offset, see
    :py:func:`precompute.draw_rand_floats`

    :param rand_floats: Random values, see :py:func:`precompute.draw_rand_floats`
    :type rand_floats: dict
    :param node: Index of the IMF node
    :type node: int

    :return: Uniform values (a view)
    :rtype: np.ndarray
    """
    imf_u = rand_floats["imf_u"]
    off = int(rand_floats["imf_offsets"][node])
    return imf_u[off : off + imf_u.size // 2]


# Largest ratio of the magnitude absorption to Av of the GAIADR3 law (G band, over the
# colors of the PARSEC isochrones: 0.63 - 1.005, decreasing with the absorption). Used
# by lost_fraction() as an upper limit when Av=0 (the ratio cannot be measured)
GAIADR3_K_MAX = 1.01
# Largest expected fraction of lost candidates used to size the initial pool (at most
# 10x the stars required). Larger losses are handled by growing the pool
Q_LOST_MAX = 0.9
# The fraction of lost candidates is only estimated (see lost_fraction()) for
# synthetic clusters of at least this many stars: for smaller ones the cost of the
# estimate is larger than that of the candidates it saves (without numba). Below it,
# the fraction Q_LOST_DEFAULT is used (an initial pool of 1.5x the stars required, +10)
POOL_ESTIMATE_MIN_STARS = 1000
Q_LOST_DEFAULT = 0.3
# Largest number of values in one bucket of the bucket sort used to order the
# synthetic stars by magnitude (see _faint_first_order_nb()); above it a stable
# np.argsort() is used. Synthetic clusters have ~10-30 at most
BUCKET_SORT_MAX = 64


def initial_pool_size(N_synth_stars: int, q_lost: float) -> int:
    """Initial number of candidate stars generated by :py:func:`synth_one`: the
    number required if a fraction ``q_lost`` of them are lost (see
    :py:func:`lost_fraction`), plus a 5% margin and 10 stars. Without differential
    reddening (``q_lost=0``) only stars interpolated close to the magnitude limit can
    be lost, which the margin covers. If the pool is still too small it is grown
    (see :py:func:`next_pool_size`); the result does not depend on its size.

    :param N_synth_stars: Number of synthetic stars required
    :type N_synth_stars: int
    :param q_lost: Expected fraction of candidates fainter than the maximum magnitude
    :type q_lost: float

    :return: Number of candidates
    :rtype: int
    """
    N_eff = N_synth_stars
    if q_lost > 0.0:
        N_eff = int(np.ceil(N_synth_stars / (1.0 - min(q_lost, Q_LOST_MAX))))
    return N_eff + N_eff // 20 + 10


def lost_fraction(
    mag_av: np.ndarray,
    mag_0: np.ndarray,
    av: float,
    dr: float,
    k_av0: float,
    k_max: float,
    imf_F: np.ndarray,
    i_lo: int,
    i_hi: int,
    max_mag: float,
    f_q: np.ndarray,
) -> float:
    """Expected fraction of the candidate stars of :py:func:`synth_one` that end up
    fainter than the maximum magnitude because of their differential reddening.

    The candidates are drawn from the IMF between the isochrone points ``i_lo`` and
    ``i_hi`` (see :py:func:`draw_candidates`), and each one gets the absorption
    ``Av + f * DR``. A star at an isochrone point with magnitude ``mag_av`` (at the
    absorption ``Av``) is lost if ``f * DR * k > max_mag - mag_av``, where ``k`` is
    the ratio of the magnitude absorption to the absorption measured at ``Av`` (or
    ``k_av0`` if ``Av=0``), and its probability is obtained from the quantiles ``f_q``
    of the fractions ``f``. Each point is weighted with the IMF probability of half of
    each of its adjacent intervals. Only the points within ``k_max * DR`` of the
    maximum magnitude (``k_max >= k``) can be lost, so the others are skipped. The
    estimate is conservative: the binary systems (brighter) are ignored, and for the
    GAIADR3 law ``k`` decreases with the absorption. Also compiled for
    :py:func:`compiled.synth_one_nb` (same operations).

    :param mag_av: Magnitudes of the isochrone at the absorption ``Av``
    :type mag_av: np.ndarray
    :param mag_0: Magnitudes of the isochrone with no absorption
    :type mag_0: np.ndarray
    :param av: Minimum absorption
    :type av: float
    :param dr: Differential reddening (``> 0``)
    :type dr: float
    :param k_av0: Ratio of the magnitude absorption to the absorption, used if
        ``Av=0``
    :type k_av0: float
    :param k_max: Upper limit of the ratio of the magnitude absorption to the
        absorption
    :type k_max: float
    :param imf_F: IMF CDF at the initial masses of the isochrone
    :type imf_F: np.ndarray
    :param i_lo: Index of the first isochrone point brighter than ``max_mag``
    :type i_lo: int
    :param i_hi: Index of the last isochrone point brighter than ``max_mag``
    :type i_hi: int
    :param max_mag: Maximum magnitude
    :type max_mag: float
    :param f_q: Quantiles of the fractions of the differential reddening, from ``0``
        to ``1`` (see :py:func:`precompute.draw_rand_floats`)
    :type f_q: np.ndarray

    :return: Expected fraction of lost candidates
    :rtype: float
    """
    dF = imf_F[i_hi] - imf_F[i_lo]
    if i_hi <= i_lo or dF <= 0.0:
        return 0.0
    # Points that can be lost
    j = i_lo + np.flatnonzero(mag_av[i_lo : i_hi + 1] > max_mag - k_max * dr)
    if j.size == 0:
        return 0.0
    mags = mag_av[j]
    if av > 0.0:
        k = (mags - mag_0[j]) / av
    else:
        k = np.full(j.size, k_av0)
    # Magnitude shift of each point for f=1 (bounded away from zero)
    shift = np.maximum(k, 1e-6) * dr
    # Probability that f > (max_mag - mag) / shift
    P = 1.0 - np.searchsorted(f_q, (max_mag - mags) / shift, side="right") / f_q.size
    # IMF probability of half of each adjacent interval (twice)
    w = imf_F[np.minimum(j + 1, i_hi)] - imf_F[np.maximum(j - 1, i_lo)]
    return float(np.sum(P * w) * 0.5 / dF)


def next_pool_size(
    N_pool: int, N_kept: int, N_stars: int, N_synth_stars: int, N_max: int
) -> int:
    """Size of the next pool of candidates in :py:func:`synth_one`, when the
    current one gave too few stars above the maximum magnitude: the size required
    by the observed fraction of surviving stars, plus a 20% margin (at least 1.5x
    the current size), and at most the number of candidates available.

    :param N_pool: Current number of candidates
    :type N_pool: int
    :param N_kept: Number of candidates that survived the magnitude cut
    :type N_kept: int
    :param N_stars: Number of candidates generated (``<= N_pool``)
    :type N_stars: int
    :param N_synth_stars: Number of synthetic stars required
    :type N_synth_stars: int
    :param N_max: Number of candidates available
    :type N_max: int

    :return: Number of candidates
    :rtype: int
    """
    # 1.2 * N_synth_stars / (N_kept / N_stars), in integer arithmetic
    N_new = (12 * N_synth_stars * N_stars) // (10 * max(N_kept, 1)) + 10
    N_new = max(N_new, N_pool + N_pool // 2)
    return min(N_new, N_max)


def grid_box(met_age_dict: dict, met: float, loga: float) -> tuple[int, int, int, int]:
    """Indexes of the (z, a) values in the grid that define the box that encloses
    the ``(met, loga)`` values (all ``0`` for a grid with a single value).

    :param met_age_dict: Dictionary of metallicity and age values.
    :type met_age_dict: dict
    :param met: Metallicity
    :type met: float
    :param loga: Logarithmic age
    :type loga: float

    :returns: Indexes ``(ml, mh, al, ah)``
    :rtype: tuple[int, int, int, int]
    """
    # bisect_left() returns the same index as np.searchsorted() (side="left"), and
    # is much faster on grids this small
    met_grid, age_grid = met_age_dict["met"], met_age_dict["loga"]
    ml = mh = 0
    if len(met_grid) > 1:
        mh = min(len(met_grid) - 1, bisect_left(met_grid, met))
        # prevent negative index
        ml = max(0, mh - 1)

    al = ah = 0
    if len(age_grid) > 1:
        ah = min(len(age_grid) - 1, bisect_left(age_grid, loga))
        al = max(0, ah - 1)

    return int(ml), int(mh), int(al), int(ah)


def isochrone_at_met_age(
    theor_tracks: np.ndarray,
    met_age_dict: dict,
    m_ini_idx: int,
    z_model: float,
    a_model: float,
    ml: int,
    mh: int,
    al: int,
    ah: int,
) -> tuple[np.ndarray, int, int]:
    """Generate a new weighted isochrone using bilinear interpolation over the
    orthogonal (z, loga) grid space.

    theor_tracks = [m1, m2, .., mN]
    mX = [age1, age2, ..., ageM]
    ageX = [mag, c1, (c2), M_ini, M_b, magb, c1b, (c2b)]
    where:
    mag:  magnitude
    cX:  colors
    M_ini: initial mass
    M_b:  binary masses
    magb: binary magnitude
    cXb: colors with the binary data

    It is important that the returned isochrone is a *copy* of the
    theor_tracks[mx][ax] array if no averaging is done. Otherwise the
    'theor_tracks' are modified.

    :param theor_tracks: Array of theoretical isochrones.
    :type theor_tracks: np.ndarray
    :param met_age_dict: Dictionary of metallicity and age values.
    :type met_age_dict: dict
    :param m_ini_idx: Index of the initial mass.
    :type m_ini_idx: int
    :param z_model: Metallicity value.
    :type z_model: float
    :param a_model: Age value.
    :type a_model: float
    :param ml: Index of the lower metallicity.
    :type ml: int
    :param mh: Index of the higher metallicity.
    :type mh: int
    :param al: Index of the lower age.
    :type al: int
    :param ah: Index of the higher age.
    :type ah: int

    :returns: Weighted isochrone, and the (metallicity, age) indexes of the closest
     grid node (the masses of the weighted isochrone are taken from this node).
    :rtype: tuple[np.ndarray, int, int]
    """

    # The four points in the (z, age) grid that define the box that contains
    # the model value (z_model, a_model)
    z1 = met_age_dict["met"][ml]
    z2 = met_age_dict["met"][mh]
    a1 = met_age_dict["loga"][al]
    a2 = met_age_dict["loga"][ah]

    # Compute normalized coordinates (t, u) within the local cell bounds
    # Safe division handling for edge cases where grid width might be zero
    dz = z2 - z1
    da = a2 - a1

    t = (z_model - z1) / dz if dz > 1e-15 else 0.0
    u = (a_model - a1) / da if da > 1e-15 else 0.0

    # Compute bilinear basis weights
    # Let om_t = (1 - t) and om_u = (1 - u)
    om_t = 1.0 - t
    om_u = 1.0 - u

    # Weights
    w11 = om_t * om_u  # Weight for (ml, al)
    w12 = om_t * u  # Weight for (ml, ah)
    w21 = t * om_u  # Weight for (mh, al)
    w22 = t * u  # Weight for (mh, ah)

    # Identify the closest grid point using a geometric boundary shortcut.
    # Since t and u map strictly to [0, 1], crossing 0.5 determines the closest node.
    m0 = mh if t >= 0.5 else ml
    a0 = ah if u >= 0.5 else al

    # Core interpolation via in-place accumulation
    isochrone = theor_tracks[ml][al] * w11
    isochrone += theor_tracks[ml][ah] * w12
    isochrone += theor_tracks[mh][al] * w21
    isochrone += theor_tracks[mh][ah] * w22

    # DO NOT average the masses or their distribution will be lost.
    # Overwrite mutated mass values with the exact discrete distribution from the closest node
    isochrone[m_ini_idx] = theor_tracks[m0][a0][m_ini_idx]

    return isochrone, m0, a0


if _HAS_NUMBA:

    @njit(cache=True)
    def _isochrone_at_met_age_nb(
        theor_tracks: np.ndarray,
        ml: int,
        mh: int,
        al: int,
        ah: int,
        w11: float,
        w12: float,
        w21: float,
        w22: float,
        m0: int,
        a0: int,
        m_ini_idx: int,
    ) -> np.ndarray:
        """Compiled (numba) version of the weighted average in
        :py:func:`isochrone_at_met_age`. A single pass reads the four grid isochrones and
        writes the weighted sum, with the same operations in the same order as
        the NumPy version (identical results), without temporary arrays. The
        mass row is copied from the closest node (m0, a0) instead of being
        averaged and then overwritten.
        """
        t11 = theor_tracks[ml, al]
        t12 = theor_tracks[ml, ah]
        t21 = theor_tracks[mh, al]
        t22 = theor_tracks[mh, ah]
        t00 = theor_tracks[m0, a0]
        N_rows, N_pts = t11.shape
        isochrone = np.empty((N_rows, N_pts), dtype=t11.dtype)
        for r in range(N_rows):
            if r == m_ini_idx:
                # DO NOT average the masses (see isochrone_at_met_age())
                for k in range(N_pts):
                    isochrone[r, k] = t00[r, k]
            else:
                for k in range(N_pts):
                    v = t11[r, k] * w11
                    v += t12[r, k] * w12
                    v += t21[r, k] * w21
                    v += t22[r, k] * w22
                    isochrone[r, k] = v
        return isochrone


def apply_distance(isochrone: np.ndarray, dm: float) -> np.ndarray:
    """Receives an isochrone of a given age and metallicity and modifies
    its magnitude values according to a given distance modulus.

    The filters stored after the initial mass (see
    :py:func:`precompute.add_color_filters`) are not moved: they are only used to
    combine the photometry of the binary systems, which depends on the differences
    between the magnitudes of their components.

    :param isochrone: Isochrone array.
    :type isochrone: np.ndarray
    :param dm: Distance modulus.
    :type dm: float

    :returns: Modified isochrone array.
    :rtype: np.ndarray
    """
    isochrone = np.ascontiguousarray(isochrone)
    # Move magnitude
    isochrone[0] += dm

    return isochrone


def apply_extinction(
    ext_law: str,
    ext_coefs: np.ndarray,
    dr_frac: np.ndarray,
    Av: float,
    dr: float,
    Rv: float,
    isochrone: np.ndarray,
) -> np.ndarray:
    """Modifies magnitude and color(s) according to given values for the
    total absorption Av. Using this parameter instead of the E(B-V) extinction
    reduces the correlation with Rv.

    The distance modulus was applied before this function.

    Only the magnitude and color(s) (the first rows of ``isochrone``) are modified:
    the rows after them (initial mass, filters, secondary mass) are not.

    For the CCMO model ``ext_coefs`` has a row ``[a, b, 0, 0]`` for the magnitude
    and a row ``[a1, b1, a2, b2]`` for each color (see
    :py:func:`precompute.ccmo_filter_coefs`), where the coefficient of each filter is:

    ccm_coef = a + b / Rv

    With differential reddening (dr>0) the total absorption of each column (a
    synthetic star, see :py:func:`synth_one`) is:

    A_i = Av + f_i * dr

    where 'f_i' (in [0, 1]) are the fractions in 'dr_frac', drawn from the
    distribution selected with 'DR_distribution' (see
    :py:func:`precompute.dr_fractions`). Hence Av is the minimum and
    Av + dr the maximum absorption.

    Ax = ef * Av
    m_obs = M_int + Ax + dist_mod
          = M_int + ef * R_V * E(B-V) + dist_mod
          = M_int + (a + b / Rv) * R_V * E(B-V) + dist_mod

    E(m1 - m2) = A_m1 - A_m2
               = (ef_m1 - ef_m2) * Av
               = [(a1 + b1/Rv) - (a2 + b2/Rv)] * Av
               = [(a1 - a2) + (b1 - b2)/Rv] * Av
               = (a12 + b12/Rv) * Av
               = (a12 + b12/Rv) * R_V * E(B-V)
    (m1 - m2)_obs = (m1 - m2)_int + E(m1 - m2)
    (m1 - m2)_obs = (m1 - m2)_int + (a12 + b12/Rv) * R_V * E(B-V)

    :param ext_law: Extinction law to be used.
    :type ext_law: str
    :param ext_coefs: Extinction coefficients (empty for ``GAIADR3``).
    :type ext_coefs: np.ndarray
    :param dr_frac: Fractions of the differential reddening applied to each column.
    :type dr_frac: np.ndarray
    :param Av: Total absorption.
    :type Av: float
    :param dr: Differential reddening.
    :type dr: float
    :param Rv: Total-to-selective extinction ratio.
    :type Rv: float
    :param isochrone: Isochrone array.
    :type isochrone: np.ndarray

    :raises ValueError: If the extinction law is not recognized.

    :returns: Modified isochrone array.
    :rtype: np.ndarray
    """
    isochrone = np.ascontiguousarray(isochrone)

    if _HAS_NUMBA and ext_law == "GAIADR3":
        # Compiled version: identical results, see _apply_extinction_gaiadr3_nb()
        # dr_kind: 0 -> no differential reddening (scalar A), 1 -> A = Av + f*dr
        dr_kind = 1 if dr > 0.0 else 0
        _apply_extinction_gaiadr3_nb(isochrone, dr_frac, dr_kind, Av, Av, dr)
        return isochrone

    if _HAS_NUMBA and ext_law == "CCMO":
        # Compiled version: identical results, see _apply_extinction_ccmo_nb()
        ec_mag, ec_col1, ec_col2 = ccmo_ratios(ext_coefs, Rv)
        _apply_extinction_ccmo_nb(
            isochrone,
            dr_frac,
            1 if dr > 0.0 else 0,
            Av,
            dr,
            ec_mag,
            ec_col1,
            ec_col2,
            len(ext_coefs) > 2,
        )
        return isochrone

    Av_dr = Av
    if dr > 0.0:
        Ns = isochrone.shape[-1]
        # A_i = Av + f_i * dr, in [Av, Av + dr]
        Av_dr = dr_frac[:Ns] * dr
        Av_dr += Av

    ec_col2 = 0.0
    if ext_law == "CCMO":
        ec_mag, ec_col1, ec_col2 = ccmo_ratios(ext_coefs, Rv)
        Ax = ec_mag * Av_dr
        Ex1 = ec_col1 * Av_dr
    elif ext_law == "GAIADR3":
        # If this model is used the first color is always expected to be BP-RP
        # BP_RP = isochrone[1]
        Ax, Ex1 = gaiadr3_extinction(isochrone[1], Av_dr)
    else:
        raise ValueError(f"Unknown extinction law: {ext_law}")

    isochrone[0] += Ax  # Magnitude
    isochrone[1] += Ex1  # First color

    # Second color (CCMO only)
    if len(ext_coefs) > 2:
        Ex2 = ec_col2 * Av_dr
        isochrone[2] += Ex2

    return isochrone


if _HAS_NUMBA:

    @njit(cache=True)
    def _apply_extinction_gaiadr3_nb(
        isochrone: np.ndarray,
        dr_frac: np.ndarray,
        dr_kind: int,
        A_scalar: float,
        Av: float,
        dr: float,
    ) -> None:
        """Compiled (numba) version of the GAIADR3 branch of
        :py:func:`apply_extinction`, modifying 'isochrone' in place.

        A single loop over the isochrone points computes the (differentially
        reddened) absorption, A_G and E(BP-RP) (same nested form as
        :py:func:`_gaiadr3_poly` and :py:func:`_gaiadr3_poly_fixed_A`), and adds them to
        the magnitude and color rows. The floating point operations are the same
        and in the same order as in the NumPy version, so the results are
        identical.

        dr_kind: 0 -> scalar absorption 'A_scalar' (no differential reddening),
        1 -> differential reddening A = Av + f * dr, with f from 'dr_frac'.
        """
        g = _G_COEFFS
        b = _BPRP_COEFFS

        # Scalar absorption: coefficients of the cubic polynomial in X
        A = A_scalar
        gk0 = A * (g[0] + A * (g[4] + A * (g[5] + A * g[6])))
        gk1 = A * (g[1] + A * (g[7] + A * g[9]))
        gk2 = A * (g[2] + A * g[8])
        gk3 = A * g[3]
        bk0 = A * (b[0] + A * (b[4] + A * (b[5] + A * b[6])))
        bk1 = A * (b[1] + A * (b[7] + A * b[9]))
        bk2 = A * (b[2] + A * b[8])
        bk3 = A * b[3]

        for i in range(isochrone.shape[1]):
            # Clamped to the range of validity of the coefficients, as in
            # gaiadr3_extinction()
            X = isochrone[1, i]
            if X < _X_MIN:
                X = _X_MIN
            elif X > _X_MAX:
                X = _X_MAX

            if dr_kind == 0:
                Ax = gk3 * X
                Ax += gk2
                Ax *= X
                Ax += gk1
                Ax *= X
                Ax += gk0
                Ex = bk3 * X
                Ex += bk2
                Ex *= X
                Ex += bk1
                Ex *= X
                Ex += bk0
            else:
                A = dr_frac[i] * dr
                A += Av

                # G: A * [P0(X) + A * (P1(X) + A * (P2(X) + c6 * A))]
                R = g[9] * X
                R += g[5]
                R += g[6] * A
                P = g[8] * X
                P += g[7]
                P *= X
                P += g[4]
                R *= A
                R += P
                P = X * g[3]
                P += g[2]
                P *= X
                P += g[1]
                P *= X
                P += g[0]
                R *= A
                R += P
                R *= A
                Ax = R
                # BP-RP
                R = b[9] * X
                R += b[5]
                R += b[6] * A
                P = b[8] * X
                P += b[7]
                P *= X
                P += b[4]
                R *= A
                R += P
                P = X * b[3]
                P += b[2]
                P *= X
                P += b[1]
                P *= X
                P += b[0]
                R *= A
                R += P
                R *= A
                Ex = R

            isochrone[0, i] += Ax  # Magnitude
            isochrone[1, i] += Ex  # First color

    @njit(cache=True)
    def _apply_extinction_ccmo_nb(
        isochrone: np.ndarray,
        dr_frac: np.ndarray,
        dr_kind: int,
        Av: float,
        dr: float,
        ec_mag: float,
        ec_col1: float,
        ec_col2: float,
        col2_flag: bool,
    ) -> None:
        """Compiled (numba) version of the CCMO branch of :py:func:`apply_extinction`,
        modifying 'isochrone' in place. The coefficients ``ec_*`` (``a + b / Rv``
        values, see :py:func:`ccmo_ratios`) multiply the (differentially reddened)
        absorption of each point. The floating point operations are the same and in
        the same order as in the NumPy version, so the results are identical.

        dr_kind: 0 -> scalar absorption Av (no differential reddening),
        1 -> differential reddening A = Av + f * dr, with f from 'dr_frac'.
        """
        A = Av
        for i in range(isochrone.shape[1]):
            if dr_kind == 1:
                A = dr_frac[i] * dr
                A += Av
            Ax = ec_mag * A
            Ex1 = ec_col1 * A
            isochrone[0, i] += Ax  # Magnitude
            isochrone[1, i] += Ex1  # First color
            if col2_flag:
                Ex2 = ec_col2 * A
                isochrone[2, i] += Ex2  # Second color


def ccmo_ratios(ext_coefs: np.ndarray, Rv: float) -> tuple[float, float, float]:
    """Coefficients of the CCMO extinction law for a given ``Rv``: the ratios
    ``A_mag / Av`` and ``E(color) / Av`` for the magnitude and the color(s), see
    :py:func:`apply_extinction`. Same operations as in ``compiled.synth_one_nb()``.

    :param ext_coefs: Extinction coefficients, see
        :py:func:`precompute.ccmo_filter_coefs`
    :type ext_coefs: np.ndarray
    :param Rv: Total-to-selective extinction ratio.
    :type Rv: float

    :returns: Coefficients for the magnitude, the first color and the second color
        (``0`` if there is no second color).
    :rtype: tuple[float, float, float]
    """
    c = ext_coefs
    ec_mag = c[0, 0] + c[0, 1] / Rv
    ec_col1 = (c[1, 0] + c[1, 1] / Rv) - (c[1, 2] + c[1, 3] / Rv)
    ec_col2 = 0.0
    if len(c) > 2:
        ec_col2 = (c[2, 0] + c[2, 1] / Rv) - (c[2, 2] + c[2, 3] / Rv)
    return float(ec_mag), float(ec_col1), float(ec_col2)


def gaiadr3_extinction(
    X: np.ndarray, A: float | np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Gaia DR3 extinction law for main sequence stars, returning the extinction
    in G and the color excess in BP-RP directly:

    A_G = ec_G(X, A) * A ; E(BP-RP) = ec_BPRP(X, A) * A

    where ec_G, ec_BPRP are the coefficients of the law: polynomials in X and A with
    the coefficients in ``_G_COEFFS`` and ``_BPRP_COEFFS``,

    ec = c0 + c1*X + c2*X^2 + c3*X^3 + c4*A + c5*A^2 + c6*A^3 + c7*A*X + c8*A*X^2
         + c9*X*A^2

    The result is equal (to floating point rounding, ~1e-15 mag) to evaluating these
    polynomials and multiplying them by ``A``, but it is evaluated faster by
    writing the polynomial in nested (Horner) form, both in X and A, with
    in-place operations to avoid allocating temporary arrays. If ``A`` is a scalar
    (no differential reddening) the powers of ``A`` are absorbed into the
    coefficients and only a cubic polynomial in X is evaluated.

    The coefficients are evaluated with the colors clamped to the range for which
    they were fitted, ``[_X_MIN, _X_MAX]``: stars outside of it are assigned the
    coefficients of the closest valid color.

    :param X: Array of (intrinsic) BP-RP colors.
    :type X: np.ndarray
    :param A: Total absorption (eventually containing differential reddening).
    :type A: float | np.ndarray

    :returns: Extinction in G and color excess in BP-RP.
    :rtype: tuple[np.ndarray, np.ndarray]
    """
    X = np.clip(X, _X_MIN, _X_MAX)
    if isinstance(A, np.ndarray):
        return _gaiadr3_poly(_G_COEFFS, X, A), _gaiadr3_poly(_BPRP_COEFFS, X, A)
    return _gaiadr3_poly_fixed_A(_G_COEFFS, X, A), _gaiadr3_poly_fixed_A(
        _BPRP_COEFFS, X, A
    )


def _gaiadr3_poly(c: tuple, X: np.ndarray, A: np.ndarray) -> np.ndarray:
    """Evaluate ``A * ec(X, A)`` for array ``A``, written as:

    A * [P0(X) + A * (P1(X) + A * (P2(X) + c6 * A))]

    with P0 = c0 + c1*X + c2*X^2 + c3*X^3 ; P1 = c4 + c7*X + c8*X^2 ;
    P2 = c5 + c9*X

    :param c: The 10 coefficients of the extinction law polynomial
    :type c: tuple
    :param X: Array of (intrinsic) BP-RP colors
    :type X: np.ndarray
    :param A: Array of absorptions, one per color
    :type A: np.ndarray

    :return: Extinction (or color excess) for each color
    :rtype: np.ndarray
    """
    # P2(X) + c6 * A
    R = c[9] * X
    R += c[5]
    R += c[6] * A
    # P1(X)
    P = c[8] * X
    P += c[7]
    P *= X
    P += c[4]
    R *= A
    R += P
    # P0(X), reusing the 'P' buffer
    np.multiply(X, c[3], out=P)
    P += c[2]
    P *= X
    P += c[1]
    P *= X
    P += c[0]
    R *= A
    R += P
    R *= A
    return R


def _gaiadr3_poly_fixed_A(c: tuple, X: np.ndarray, A: float) -> np.ndarray:
    """Evaluate ``A * ec(X, A)`` for scalar ``A`` as a cubic polynomial in X:

    k0 + k1*X + k2*X^2 + k3*X^3

    where the 'k' coefficients contain all the powers of A.

    :param c: The 10 coefficients of the extinction law polynomial
    :type c: tuple
    :param X: Array of (intrinsic) BP-RP colors
    :type X: np.ndarray
    :param A: Absorption, the same for all colors
    :type A: float

    :return: Extinction (or color excess) for each color
    :rtype: np.ndarray
    """
    k0 = A * (c[0] + A * (c[4] + A * (c[5] + A * c[6])))
    k1 = A * (c[1] + A * (c[7] + A * c[9]))
    k2 = A * (c[2] + A * c[8])
    k3 = A * c[3]
    R = k3 * X
    R += k2
    R *= X
    R += k1
    R *= X
    R += k0
    return R


def bright_points(
    isoch_moved: np.ndarray, max_mag_syn: float
) -> tuple[np.ndarray, np.ndarray]:
    """Identify the isochrone points with magnitude values smaller than the maximum
    observed value.

    The isochrone is not cut here: ``isoch_moved[:, idx]`` would copy the entire
    array, which is not needed since :py:func:`draw_candidates` only reads the
    points it interpolates. Instead, the indexes of the kept points are returned.

    :param isoch_moved: Isochrone array.
    :type isoch_moved: np.ndarray
    :param max_mag_syn: Maximum magnitude value.
    :type max_mag_syn: float

    :returns: Indexes of the points kept (sorted), and the boolean mask of the
     points kept.
    :rtype: tuple[np.ndarray, np.ndarray]
    """
    # Discard stars in isochrone beyond max_mag_syn limit.
    msk = isoch_moved[0] < max_mag_syn
    return np.flatnonzero(msk), msk


def imf_masses(
    mass_ini: np.ndarray,
    imf_F: np.ndarray,
    F_lo: float,
    dF: float,
    v_min: float,
    v_max: float,
    u: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Masses of the candidate stars drawn from the IMF restricted to ``[v_min,
    v_max]``: each uniform value ``u`` gives the CDF value ``F = F_lo + u * dF`` (with
    ``F_lo, dF`` the CDF at ``v_min`` and its range), whose mass is interpolated
    linearly between the two grid points of the isochrone that bracket it.

    Same operations, in the same order, as the compiled version in
    :py:func:`_draw_candidates_nb`, so the results are identical.

    :param mass_ini: Initial masses of the isochrone (sorted)
    :type mass_ini: np.ndarray
    :param imf_F: IMF CDF at ``mass_ini``, see :py:func:`precompute.imf_cdf_tables`
    :type imf_F: np.ndarray
    :param F_lo: IMF CDF at ``v_min``
    :type F_lo: float
    :param dF: IMF CDF at ``v_max`` minus ``F_lo``
    :type dF: float
    :param v_min: Smallest mass
    :type v_min: float
    :param v_max: Largest mass
    :type v_max: float
    :param u: Uniform values
    :type u: np.ndarray

    :return: Indexes of the first grid point with a CDF value not smaller than that
        of each candidate (``np.searchsorted(imf_F, F)``), and the masses
    :rtype: tuple[np.ndarray, np.ndarray]
    """
    F = F_lo + u * dF
    s = np.searchsorted(imf_F, F)
    N_g = imf_F.size
    s0 = np.clip(s - 1, 0, N_g - 1)
    s1 = np.clip(s, 0, N_g - 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        m = mass_ini[s0] + (F - imf_F[s0]) * (mass_ini[s1] - mass_ini[s0]) / (
            imf_F[s1] - imf_F[s0]
        )
    m[s == 0] = mass_ini[0]
    m[s == N_g] = mass_ini[N_g - 1]
    m = np.minimum(np.maximum(m, v_min), v_max)
    return s, m


def draw_candidates(
    isochrone: np.ndarray,
    idx_cut: np.ndarray,
    msk_cut: np.ndarray,
    m_ini_idx: int,
    imf_u: np.ndarray,
    N_pool: int,
    imf_F: np.ndarray,
) -> np.ndarray:
    """Draw the masses of (up to) ``N_pool`` candidate stars and interpolate them
    (and all the sub-arrays of the isochrone) into the isochrone.

    The masses are drawn from the IMF restricted to the mass range ``[v_min,
    v_max]`` of the isochrone points given by ``idx_cut`` (those brighter than the
    maximum magnitude, see :py:func:`bright_points`): no star outside of that range
    can be brighter than the maximum magnitude, so every candidate can end up in the
    synthetic cluster (see :py:func:`imf_masses`). Only those isochrone points are
    used for the interpolation, without having to copy the cut isochrone.

    :param isochrone: Isochrone array (not cut).
    :type isochrone: np.ndarray
    :param idx_cut: Sorted indexes of the isochrone points kept by
        :py:func:`bright_points`.
    :type idx_cut: np.ndarray
    :param msk_cut: Boolean mask of the isochrone points kept by
        :py:func:`bright_points`.
    :type msk_cut: np.ndarray
    :param m_ini_idx: Index of the initial mass.
    :type m_ini_idx: int
    :param imf_u: Uniform values of the candidates, see :py:func:`imf_uniforms`
    :type imf_u: np.ndarray
    :param N_pool: Number of candidates.
    :type N_pool: int
    :param imf_F: IMF CDF at the initial masses of the isochrone, see
        :py:func:`precompute.imf_cdf_tables`
    :type imf_F: np.ndarray

    :returns: Interpolated isochrone array (empty if the mass range is empty)
    :rtype: np.ndarray
    """
    # ==IMPORTANT==
    # Assumes `mass_ini=isochrone[m_ini_idx]` is sorted min to max. This is checked
    # when loading the isochrones in isochrone_loader.interp_isochrones()
    # ==IMPORTANT==
    mass_ini = isochrone[m_ini_idx]
    i_lo, i_hi = idx_cut[0], idx_cut[-1]
    F_lo = float(imf_F[i_lo])
    dF = float(imf_F[i_hi]) - F_lo
    # A single isochrone point (or an empty CDF range): no mass in the range
    if idx_cut.size < 2 or dF <= 0.0:
        return np.array([])

    N = min(N_pool, imf_u.size)
    s, mass_dist = imf_masses(
        mass_ini, imf_F, F_lo, dF, mass_ini[i_lo], mass_ini[i_hi], imf_u[:N]
    )

    # Indexes into the *cut* mass grid: the number of grid points before 's' that
    # survived the cut
    N_kept = np.empty(msk_cut.size + 1, dtype=np.intp)
    N_kept[0] = 0
    np.cumsum(msk_cut, out=N_kept[1:])
    hi_cut = N_kept[s]

    # Map back to indexes into the full isochrone. Note that hi_cut=0 gives
    # lo=idx_cut[-1], exactly as indexing the cut isochrone with -1 would
    x_new_indices = idx_cut[hi_cut]
    lo = idx_cut[hi_cut - 1]

    # Interpolate the candidates (masses) into the isochrone
    isoch_mass = interp_at_masses(isochrone, mass_ini, mass_dist, x_new_indices, lo)

    return isoch_mass


if _HAS_NUMBA:

    @njit(cache=True)
    def _draw_candidates_nb(
        isochrone: np.ndarray,
        idx_cut: np.ndarray,
        msk_cut: np.ndarray,
        m_ini_idx: int,
        imf_u: np.ndarray,
        N_pool: int,
        imf_F: np.ndarray,
        imf_lut: np.ndarray,
    ) -> np.ndarray:
        """Compiled (numba) version of :py:func:`draw_candidates`. It performs the same
        floating point operations in the same order, so the results are identical
        to the NumPy version, but it avoids most temporary arrays and Python calls:

        1. Bracketing points: for each point 'k' of the full mass grid, the next
           kept point at or after 'k' and the last kept point before 'k' are
           tabulated once.
        2. For each candidate, the grid index ``np.searchsorted(imf_F, F)`` is
           found by a short walk from the bucket lookup table 'imf_lut' (see
           :py:func:`precompute.imf_cdf_tables`), and its mass as in
           :py:func:`imf_masses`.
        3. Interpolate the isochrone linearly, one row at a time (rows are
           contiguous in memory).
        """
        mass_ini = isochrone[m_ini_idx]
        N_cut = idx_cut.size
        i_lo, i_hi = idx_cut[0], idx_cut[N_cut - 1]
        v_min, v_max = mass_ini[i_lo], mass_ini[i_hi]
        F_lo = imf_F[i_lo]
        dF = imf_F[i_hi] - F_lo
        N_rows = isochrone.shape[0]
        if N_cut < 2 or dF <= 0.0:
            return np.empty((N_rows, 0), dtype=isochrone.dtype)
        N = min(N_pool, imf_u.size)
        N_g = imf_F.size
        N_b = imf_lut.size - 1

        # 1. Next kept point at or after k (nxt), and last kept point before k
        # (prv). If there is no kept point before k, prv=idx_cut[-1], as indexing
        # the cut isochrone with -1 does in the NumPy version.
        N_grid = msk_cut.size
        nxt = np.empty(N_grid + 1, dtype=np.intp)
        prv = np.empty(N_grid + 1, dtype=np.intp)
        last = idx_cut[N_cut - 1]
        for k in range(N_grid + 1):
            prv[k] = last
            if k < N_grid and msk_cut[k]:
                last = k
        nxt_k = N_grid
        for k in range(N_grid, -1, -1):
            if k < N_grid and msk_cut[k]:
                nxt_k = k
            nxt[k] = nxt_k

        # 2. Grid index and mass of each candidate
        hi = np.empty(N, dtype=np.intp)
        lo = np.empty(N, dtype=np.intp)
        w = np.empty(N, dtype=isochrone.dtype)
        for j in range(N):
            F = F_lo + imf_u[j] * dF
            s = np.intp(imf_lut[min(int(F * N_b), N_b)])
            while s < N_g and imf_F[s] < F:
                s += 1
            while s > 0 and imf_F[s - 1] >= F:
                s -= 1
            if s == 0:
                m = mass_ini[0]
            elif s == N_g:
                m = mass_ini[N_g - 1]
            else:
                m = mass_ini[s - 1] + (F - imf_F[s - 1]) * (
                    mass_ini[s] - mass_ini[s - 1]
                ) / (imf_F[s] - imf_F[s - 1])
            m = min(max(m, v_min), v_max)
            hi[j] = nxt[s]
            lo[j] = prv[s]
            x_lo = mass_ini[lo[j]]
            w_j = m - x_lo
            w_j /= mass_ini[hi[j]] - x_lo
            w[j] = w_j

        # 3. Linear interpolation, y_lo + (y_hi - y_lo) * w
        isoch_mass = np.empty((N_rows, N), dtype=isochrone.dtype)
        for r in range(N_rows):
            row = isochrone[r]
            out = isoch_mass[r]
            for j in range(N):
                y_lo = row[lo[j]]
                out[j] = (row[hi[j]] - y_lo) * w[j] + y_lo

        return isoch_mass


def interp_at_masses(
    isoch_cut: np.ndarray,
    mass_ini: np.ndarray,
    mass_dist: np.ndarray,
    x_new_indices: np.ndarray | None = None,
    lo: np.ndarray | None = None,
) -> np.ndarray:
    """Find where in the original data, the values to interpolate would be inserted.

    NOTE: I already tried to speed this block up using numba (@jit(nopython=True))
    but it does not help. The code runs even slower.

    :param isoch_cut: Isochrone array.
    :type isoch_cut: np.ndarray
    :param mass_ini: Array of initial masses.
    :type mass_ini: np.ndarray
    :param mass_dist: Array of sampled masses.
    :type mass_dist: np.ndarray
    :param x_new_indices: Result of ``np.searchsorted(mass_ini, mass_dist)`` if
        already known; computed here otherwise.
    :type x_new_indices: np.ndarray | None
    :param lo: Indexes of the lower bracketing points, ``x_new_indices - 1`` if
        not given. Used when the isochrone points are selected through indexes
        (see :py:func:`draw_candidates`).
    :type lo: np.ndarray | None

    :returns: Interpolated isochrone array.
    :rtype: np.ndarray
    """
    if x_new_indices is None:
        x_new_indices = np.searchsorted(mass_ini, mass_dist)
    #
    # Calculate the slope of regions that each x_new value falls in.
    if lo is None:
        lo = x_new_indices - 1
    x_lo = mass_ini[lo]
    # Interpolation weights: (x_new - x_lo) / (x_hi - x_lo)
    w = mass_dist - x_lo
    w /= mass_ini[x_new_indices] - x_lo

    # Calculate the actual value for each entry in x_new as:
    # y_lo + (y_hi - y_lo) * w
    # using in-place operations to avoid allocating temporary (D, N) arrays
    y_lo = np.take(isoch_cut, lo, axis=1)
    isoch_mass = np.take(isoch_cut, x_new_indices, axis=1)
    isoch_mass -= y_lo
    isoch_mass *= w
    isoch_mass += y_lo

    return isoch_mass


def assign_binaries(
    alpha: float,
    beta: float,
    binar_flag: bool,
    m_ini_idx: int,
    rand_unif_vals: np.ndarray,
    binar_tables: tuple,
    isochrone: np.ndarray,
    isoch_mass: np.ndarray,
) -> np.ndarray:
    """Select a fraction of stars to be binaries, and obtain the photometry of the
    binary systems.

    Each star is a binary system with probability ``alpha + beta / (1 + 1.4 / m1)``
    (Offner et al. 2022; Fig 1 (left), Table 1), with ``m1`` its (primary) mass. The
    mass ratio ``q`` of each binary system is obtained from its own uniform random
    value and its primary mass (see :py:func:`precompute.binary_tables`), and the
    secondary mass ``m2 = q * m1`` is interpolated linearly in the full isochrone
    (moved by the distance modulus, before the extinction is applied). The flux of
    the secondary is added to that of the primary in the magnitude and in both
    filters of each color, the second one being the first minus the color. A
    secondary with a mass below the minimum mass of the isochrone adds no flux.

    With numba installed the compiled version is used (the same function called by
    the compiled ``compiled.synth_one_nb()``), else a NumPy version with the same
    operations.

    :param alpha: Alpha parameter for the binary fraction.
    :type alpha: float
    :param beta: Beta parameter for the binary fraction.
    :type beta: float
    :param binar_flag: Flag to indicate if binarity is being used.
    :type binar_flag: bool
    :param m_ini_idx: Index of the initial mass.
    :type m_ini_idx: int
    :param rand_unif_vals: Uniform random values of the candidate stars (binarity)
    :type rand_unif_vals: np.ndarray
    :param binar_tables: Data used to obtain the mass ratios and the photometry of
        the binary systems, see :py:func:`precompute.binary_tables`
    :type binar_tables: tuple
    :param isochrone: Full isochrone moved by the distance modulus, with rows
        ``[mag, c1, (c2), mass, f1, (f2)]`` (``f*``: first filter of each color)
    :type isochrone: np.ndarray
    :param isoch_mass: Candidate stars interpolated in the isochrone, same rows as
        ``isochrone``. Modified in place
    :type isoch_mass: np.ndarray

    :returns: Array of the stars, with rows ``[mag, c1, (c2), mass, mass_b]``, where
        ``mass_b`` is the secondary mass (``nan`` for single systems)
    :rtype: np.ndarray
    """
    # No binarity process defined
    if binar_flag is False:
        # Update the binary systems' masses so that the secondary masses for
        # SINGLE systems are identified with a 'nan' value.
        isoch_mass[m_ini_idx + 1] = np.nan
        return isoch_mass[: m_ini_idx + 2]

    q_pre, log_u, e_tab, e_m_lo, e_inv_step, flux_tab, flux_inv_step = binar_tables

    if _HAS_NUMBA:
        _assign_binaries_nb(
            alpha,
            beta,
            m_ini_idx,
            rand_unif_vals,
            q_pre,
            log_u,
            e_tab,
            e_m_lo,
            e_inv_step,
            flux_tab,
            flux_inv_step,
            isochrone,
            isoch_mass,
        )
        return isoch_mass[: m_ini_idx + 2]

    Ns = isoch_mass.shape[-1]
    m1 = isoch_mass[m_ini_idx]

    # Distribution of probability of binarity (multiplicity fraction) versus
    # primary masses: alpha + beta / (1.0 + 1.4 / m1). Clipping to [0, 1] is not
    # needed: the probabilities are only compared with uniform random values in
    # [0, 1), so a value < 0 (> 1) gives the same result as 0 (1)
    b_p = 1.4 / m1
    b_p += 1.0
    b_p = beta / b_p
    b_p += alpha
    b_idx = np.flatnonzero(b_p > rand_unif_vals[:Ns])

    # Secondary masses
    m1_b = m1[b_idx]
    if e_tab.size > 0:
        q = np.exp(log_u[b_idx] * _table_interp(m1_b, e_tab, e_m_lo, e_inv_step))
    else:
        q = q_pre[b_idx]
    m2 = q * m1_b

    # Secondaries with a mass in the range of the isochrone, interpolated in it
    mass_iso = isochrone[m_ini_idx]
    flux = m2 >= mass_iso[0]
    j_f, m2_f = b_idx[flux], m2[flux]
    s = np.minimum(np.maximum(np.searchsorted(mass_iso, m2_f), 1), mass_iso.size - 1)
    lo = s - 1
    x_lo = mass_iso[lo]
    w = (m2_f - x_lo) / (mass_iso[s] - x_lo)

    def interp_row(r: int) -> np.ndarray:
        row = isochrone[r]
        return (row[s] - row[lo]) * w + row[lo]

    def flux_add(dmag: np.ndarray) -> np.ndarray:
        """Change of a magnitude when the flux of a star 'dmag' fainter is added"""
        neg = dmag < 0.0
        out = _table_interp(np.abs(dmag), flux_tab, 0.0, flux_inv_step)
        out[neg] += dmag[neg]
        return out

    for k in range(m_ini_idx - 1):
        c1, f1 = isoch_mass[1 + k, j_f], isoch_mass[m_ini_idx + 1 + k, j_f]
        c2, f2 = interp_row(1 + k), interp_row(m_ini_idx + 1 + k)
        # First filter, and second filter (first filter minus color)
        d_f1 = f2 - f1
        d_f2 = (f2 - c2) - (f1 - c1)
        isoch_mass[1 + k, j_f] = c1 + (flux_add(d_f1) - flux_add(d_f2))
    mag1 = isoch_mass[0, j_f]
    isoch_mass[0, j_f] = mag1 + flux_add(interp_row(0) - mag1)

    # Secondary masses, 'nan' for single systems
    isoch_mass[m_ini_idx + 1] = np.nan
    isoch_mass[m_ini_idx + 1, b_idx] = m2

    # Return [mag, c1, (c2), mass, mass_b]
    return isoch_mass[: m_ini_idx + 2]


def _table_interp(
    v: np.ndarray, tab: np.ndarray, v_lo: float, inv_step: float
) -> np.ndarray:
    """Linear interpolation in the table ``tab`` of a function tabulated at the
    values ``v_lo + k / inv_step`` (``k = 0, 1, ...``), with the first (last) value
    of the table for ``v`` below (above) its range. NumPy version of
    ``_table_interp_nb()``, with the same operations"""
    x = (v - v_lo) * inv_step
    out = np.empty(x.size)
    below = x <= 0.0
    above = ~below & (x >= tab.size - 1)
    mid = ~(below | above)
    k = x[mid].astype(np.intp)
    y0 = tab[k]
    out[mid] = (tab[k + 1] - y0) * (x[mid] - k) + y0
    out[below] = tab[0]
    out[above] = tab[-1]
    return out


if _HAS_NUMBA:

    @njit(cache=True)
    def _scatter_add_nb(row: np.ndarray, order: np.ndarray, sigma: np.ndarray) -> None:
        """Compiled equivalent of ``row[order] += sigma[:order.size]``, used by
        :py:func:`add_errors`"""
        for j in range(order.size):
            row[order[j]] += sigma[j]

    @njit(cache=True)
    def _faint_first_order_nb(mag: np.ndarray) -> np.ndarray:
        """Compiled equivalent of ``np.argsort(-mag, kind="stable")``, used by
        :py:func:`add_errors`: indexes that sort the stars from the largest to the
        smallest magnitude (stars with equal magnitudes keep their order).

        A bucket sort, several times faster than ``np.argsort()`` for the sizes of
        the synthetic clusters: the values are distributed in ``N`` buckets of equal
        width (counting sort, which keeps their order) and each bucket is sorted by
        insertion (also stable). If too many values fall in the same bucket (values
        concentrated in a small part of their range), the quadratic cost of the
        insertion sort is avoided by using a stable ``np.argsort()`` instead.
        """
        N = mag.size
        x = -mag
        if N < 2:
            return np.arange(N)
        lo, hi = x.min(), x.max()
        scale = N / (hi - lo) if hi > lo else 0.0
        if not (scale > 0.0 and np.isfinite(scale)):
            return np.argsort(x, kind="mergesort")

        # Bucket of each value (non-decreasing with the value), and the number of
        # values in each bucket
        bucket = np.empty(N, dtype=np.intp)
        start = np.zeros(N + 1, dtype=np.intp)
        for i in range(N):
            b = min(int((x[i] - lo) * scale), N - 1)
            bucket[i] = b
            start[b + 1] += 1
        largest = 0
        for b in range(N):
            largest = max(largest, start[b + 1])
            start[b + 1] += start[b]
        if largest > BUCKET_SORT_MAX:
            return np.argsort(x, kind="mergesort")

        # Counting sort of the buckets
        order = np.empty(N, dtype=np.intp)
        pos = start[:N].copy()
        for i in range(N):
            b = bucket[i]
            order[pos[b]] = i
            pos[b] += 1
        # Insertion sort of each bucket (strict comparison: stable)
        for b in range(N):
            s, e = start[b], start[b + 1]
            for j in range(s + 1, e):
                v = order[j]
                xv = x[v]
                k = j - 1
                while k >= s and x[order[k]] > xv:
                    order[k + 1] = order[k]
                    k -= 1
                order[k + 1] = v
        return order


if _HAS_NUMBA:

    @njit(cache=True)
    def _table_interp_nb(v: float, tab: np.ndarray, v_lo: float, inv_step: float):
        """Linear interpolation in the table ``tab`` of a function tabulated at the
        values ``v_lo + k / inv_step`` (``k = 0, 1, ...``), for a single value
        ``v``, with the first (last) value of the table below (above) its range"""
        x = (v - v_lo) * inv_step
        if x <= 0.0:
            return tab[0]
        if x >= tab.size - 1:
            return tab[tab.size - 1]
        k = int(x)
        y0 = tab[k]
        return (tab[k + 1] - y0) * (x - k) + y0

    @njit(cache=True)
    def _flux_add_nb(dmag: float, flux_tab: np.ndarray, inv_step: float) -> float:
        """Change of a magnitude when the flux of a star 'dmag' fainter is added (see
        :py:func:`precompute.binary_tables`)"""
        if dmag < 0.0:
            return _table_interp_nb(-dmag, flux_tab, 0.0, inv_step) + dmag
        return _table_interp_nb(dmag, flux_tab, 0.0, inv_step)

    @njit(cache=True)
    def _assign_binaries_nb(
        alpha: float,
        beta: float,
        m_ini_idx: int,
        rand_unif_vals: np.ndarray,
        q_pre: np.ndarray,
        log_u: np.ndarray,
        e_tab: np.ndarray,
        e_m_lo: float,
        e_inv_step: float,
        flux_tab: np.ndarray,
        flux_inv_step: float,
        isochrone: np.ndarray,
        isoch_mass: np.ndarray,
    ) -> None:
        """Compiled (numba) version of :py:func:`assign_binaries` (with binaries),
        modifying 'isoch_mass' in place. A single loop over the stars computes the
        binary probability and compares it with the uniform random value. A single
        star gets a 'nan' secondary mass. For a binary system the mass ratio is
        obtained (see :py:func:`precompute.binary_tables`), the secondary is
        interpolated in the full isochrone (binary search in its masses) and its
        flux is added to the magnitude and to both filters of each color; the
        secondary mass replaces the first filter (read before) in row
        ``m_ini_idx + 1``.
        """
        mass = isoch_mass[m_ini_idx]
        mass_iso = isochrone[m_ini_idx]
        N_iso = mass_iso.size
        m_min = mass_iso[0]
        N_colors = m_ini_idx - 1
        mass_dep = e_tab.size > 0
        for j in range(isoch_mass.shape[1]):
            m1 = mass[j]
            # alpha + beta / (1.0 + 1.4 / mass)
            b_p = 1.4 / m1
            b_p += 1.0
            b_p = beta / b_p
            b_p += alpha
            if not b_p > rand_unif_vals[j]:
                isoch_mass[m_ini_idx + 1, j] = np.nan
                continue

            # Mass ratio and secondary mass
            if mass_dep:
                q = np.exp(log_u[j] * _table_interp_nb(m1, e_tab, e_m_lo, e_inv_step))
            else:
                q = q_pre[j]
            m2 = q * m1

            # Photometry of the system. A secondary below the minimum mass of the
            # isochrone adds no flux
            if m2 >= m_min:
                s = np.searchsorted(mass_iso, m2)
                s = min(max(s, 1), N_iso - 1)
                lo = s - 1
                x_lo = mass_iso[lo]
                w = (m2 - x_lo) / (mass_iso[s] - x_lo)
                for k in range(N_colors):
                    c_r = 1 + k
                    f_r = m_ini_idx + 1 + k
                    row_c, row_f = isochrone[c_r], isochrone[f_r]
                    c2 = (row_c[s] - row_c[lo]) * w + row_c[lo]
                    f2 = (row_f[s] - row_f[lo]) * w + row_f[lo]
                    c1 = isoch_mass[c_r, j]
                    f1 = isoch_mass[f_r, j]
                    # First filter, and second filter (first filter minus color)
                    d_f1 = f2 - f1
                    d_f2 = (f2 - c2) - (f1 - c1)
                    isoch_mass[c_r, j] = c1 + (
                        _flux_add_nb(d_f1, flux_tab, flux_inv_step)
                        - _flux_add_nb(d_f2, flux_tab, flux_inv_step)
                    )
                row_m = isochrone[0]
                mag2 = (row_m[s] - row_m[lo]) * w + row_m[lo]
                mag1 = isoch_mass[0, j]
                isoch_mass[0, j] = mag1 + _flux_add_nb(
                    mag2 - mag1, flux_tab, flux_inv_step
                )

            isoch_mass[m_ini_idx + 1, j] = m2


def add_errors(isoch_binar: np.ndarray, err_dist: np.ndarray) -> np.ndarray:
    """Add random synthetic uncertainties to the magnitude and color(s)

    :param isoch_binar: Isochrone array.
    :type isoch_binar: np.ndarray
    :param err_dist: Error distributions, one row per magnitude and color (no rows
        if the synthetic clusters have no uncertainties)
    :type err_dist: np.ndarray

    :returns: Isochrone array with added errors.
    :rtype: np.ndarray
    """
    if len(err_dist) == 0:
        return isoch_binar

    N = isoch_binar.shape[1]
    # Indexes that sort the stars from the largest to the smallest magnitude, the
    # same order used by the error distributions (see error_distribution()). Stars
    # with equal magnitudes keep their order (stable sort), as in the compiled
    # version (_faint_first_order_nb()), so that both give identical results
    mag_sort = np.argsort(-isoch_binar[0], kind="stable")

    for i, sigma in enumerate(err_dist):
        # Index the row first: fancy indexing a 1D view is noticeably faster
        # than the mixed (row, indexes) indexing of the 2D array
        row = isoch_binar[i]
        row[mag_sort] += sigma[:N]

    return isoch_binar

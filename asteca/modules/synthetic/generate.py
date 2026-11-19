from bisect import bisect_left
from typing import TYPE_CHECKING

import numpy as np

from .. import HAS_NUMBA as _HAS_NUMBA
from .. import njit

if TYPE_CHECKING:
    from .engine import ModelEngine

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
# batch._synth_one() (0: not empty)
EMPTY_REASONS = ("", "max magnitude cut", "IMF interpolation", "the extinction")


def get_model(
    ctx: "ModelEngine",
    p: np.ndarray,
    dr_model: str,
    calib: tuple[float, list[np.ndarray], int],
    show_errors: bool = False,
    isoch_flag: bool = False,
    isoch_dr_frac: float = 0.0,
    mass_info: dict | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Generate a synthetic cluster.

    The synthetic cluster is generated according to the parameters in ``p`` and the
    data of the :py:class:`Synthetic` object. This is the Python implementation of
    the synthetic cluster recipe, used when numba is not installed; with numba, the
    compiled ``batch._synth_one()`` is used instead (through
    ``engine.ModelEngine.get_model()``), with identical results.

    :param ctx: Data that defines the synthetic clusters
    :type ctx: ModelEngine
    :param p: Values of the fundamental parameters, in the order of ``FULL_PARAMS``:
        a row of the array returned by :py:meth:`ModelEngine.params_array
        <engine.ModelEngine.params_array>` (the values are not validated here)
    :type p: np.ndarray
    :param dr_model: Differential reddening model, it selects the fractions of the
        differential reddening applied to each star (see
        :py:func:`precompute.dr_fractions`)
    :type dr_model: str
    :param calib: Maximum synthetic magnitude allowed, error distributions used to
        perturb the synthetic photometry, and number of synthetic stars to generate
        (see ``Synthetic._calib_values()``)
    :type calib: tuple[float, list[np.ndarray], int]
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
        to estimate the total mass of the simulated cluster (see
        :py:func:`cluster_mass.get_M_actual`): the IMF sampled masses used to
        generate it (``imf_masses``: all the masses up to the last star kept in the
        synthetic cluster, or the full sample if it was exhausted), the indexes of
        the stars kept (``kept_idx``), the number of stars requested
        (``N_target``), the initial and secondary masses of the
        isochrone (``mass_row``, ``m2_row``), and the binary fraction parameters
        (``alpha``, ``beta``). It is not modified if the synthetic cluster is empty
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
    max_mag_syn, err_dist_synth, N_synth_stars = calib

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
        isochrone, m0, a0 = zaWAverage(
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

    binar_flag = True
    if alpha == 0.0 and beta == 0.0:
        binar_flag = False

        # # TODO: this was not tested thoroughly (April 2025)
        # # Remove binary photometry
        # isochrone = isochrone[: self.m_ini_idx + 2]
        # # TODO: this was not tested thoroughly

    # Move theoretical isochrone using the distance modulus
    isoch_moved = move_isochrone(isochrone, binar_flag, m_ini_idx, dm)

    # Isochrone with the minimum absorption Av (no differential reddening). It is
    # only used to identify the isochrone points that can produce stars brighter
    # than the maximum magnitude (a larger absorption only makes a star fainter),
    # and returned as the isochrone if requested. The copy avoids reddening
    # 'isoch_moved', which is interpolated below
    isoch_av = extinction(
        ext_law,
        ext_coefs,
        dr_frac,
        m_ini_idx,
        binar_flag,
        av,
        0.0,
        rv,
        isoch_moved.copy(),
    )
    isoch_extin_dr_0 = output_isochrone(
        ctx, p, dr_model, isoch_flag, isoch_dr_frac, isoch_moved, isoch_av
    )

    # Identify the isochrone points below the maximum magnitude. The isochrone
    # itself is not cut (which would copy it): the points are selected through
    # their indexes in mass_interp()
    idx_cut, msk_cut = cut_max_mag(isoch_av, max_mag_syn)
    if idx_cut.size == 0:
        if show_errors:
            print(empty_message(1, p))
        return isoch_extin_dr_0, np.array([])

    # Pre-computed search indexes for this IMF sample and grid node
    mass_idx = ctx.mass_idx[(ml, al, m0, a0)]
    st_mass = ctx.st_dist_mass[ml][al]

    # The extinction is applied to each synthetic star (after the interpolation and
    # the binarity), so that each star gets its own differential reddening draw
    # (see per_star_pool()). As some stars end up below the maximum magnitude, a
    # pool of candidates is generated: the first stars of the IMF sample within the
    # isochrone's mass range. The pool is grown until 'N_synth_stars' stars are
    # brighter than the maximum magnitude, or the IMF sample is exhausted. The
    # stars (and their random values) are always taken in the same order, so the
    # result does not depend on the pool sizes
    N_pool = per_star_pool(N_synth_stars, dr)
    while True:
        # Interpolate IMF's sampled masses into the isochrone.
        isoch_mass = mass_interp(
            isoch_moved,
            idx_cut,
            msk_cut,
            m_ini_idx,
            st_mass,
            N_pool,
            mass_idx,
        )
        if isoch_mass.size == 0:
            if show_errors:
                print(empty_message(2, p))
            return isoch_extin_dr_0, np.array([])
        N_stars = isoch_mass.shape[1]

        # Assignment of binarity.
        isoch_binar = binarity(
            alpha,
            beta,
            binar_flag,
            m_ini_idx,
            rand_floats["unif"][1],
            isoch_mass,
        )

        # Extinction of each star, from the intrinsic color of the system. The
        # binary photometry was already combined, hence 'binar_flag=False'
        isoch_binar = extinction(
            ext_law,
            ext_coefs,
            dr_frac,
            m_ini_idx,
            False,
            av,
            dr,
            rv,
            isoch_binar,
        )

        keep = np.flatnonzero(isoch_binar[0] < max_mag_syn)
        # Enough stars, or IMF sample exhausted
        if keep.size >= N_synth_stars or N_stars < N_pool:
            break
        N_new = per_star_pool_next(
            N_pool, keep.size, N_stars, N_synth_stars, st_mass.size
        )
        if N_new <= N_pool:
            break
        N_pool = N_new

    if keep.size == 0:
        if show_errors:
            print(empty_message(3, p))
        return isoch_extin_dr_0, np.array([])

    if mass_info is not None:
        fill_mass_info(
            mass_info, ctx, p, N_synth_stars, isoch_moved, idx_cut, N_pool, keep
        )

    isoch_binar = isoch_binar[:, keep[:N_synth_stars]]

    # Assign errors according to errors distribution.
    synth_clust = add_errors(isoch_binar, err_dist_synth)

    return isoch_extin_dr_0, synth_clust


def empty_message(status: int, p: np.ndarray) -> str:
    """Message printed when an empty synthetic cluster is returned and
    ``show_errors=True``

    :param status: Reason for the empty cluster, see ``EMPTY_REASONS``
    :type status: int
    :param p: Values of the fundamental parameters, see :py:func:`get_model`
    :type p: np.ndarray

    :return: Message
    :rtype: str
    """
    met, loga, _, _, av, dr, _, dm = p.tolist()
    return (
        f"No stars left after {EMPTY_REASONS[status]} (met={met:.5f},"
        + f"log(age)={loga:.3f}, Av={av:.3f}, DR={dr:.3f}, dm={dm:.3f})"
    )


def output_isochrone(
    ctx: "ModelEngine",
    p: np.ndarray,
    dr_model: str,
    isoch_flag: bool,
    isoch_dr_frac: float,
    isoch_moved: np.ndarray,
    isoch_av: np.ndarray,
) -> np.ndarray:
    """Isochrone returned by :py:func:`get_model`: the isochrone with the minimum
    absorption ``Av``, or reddened with ``Av + isoch_dr_frac * DR`` if
    ``isoch_dr_frac > 0``; an empty array if ``isoch_flag=False``

    :param ctx: Data that defines the synthetic clusters
    :type ctx: ModelEngine
    :param p: Values of the fundamental parameters, see :py:func:`get_model`
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
    _, _, alpha, beta, av, dr, rv, _ = p.tolist()
    if isoch_dr_frac > 0 and dr > 0:
        return extinction(
            ctx.ext_law,
            ctx.ext_coefs,
            ctx.rand_floats["dr_frac"][dr_model],
            ctx.m_ini_idx,
            not (alpha == 0.0 and beta == 0.0),
            av + isoch_dr_frac * dr,
            0.0,
            rv,
            isoch_moved.copy(),
        )
    return isoch_av


def fill_mass_info(
    mass_info: dict,
    ctx: "ModelEngine",
    p: np.ndarray,
    N_synth_stars: int,
    isoch_moved: np.ndarray,
    idx_cut: np.ndarray,
    N_pool: int,
    keep: np.ndarray,
) -> None:
    """Fill ``mass_info`` with the data required to estimate the total mass of a
    (non empty) simulated cluster, see :py:func:`get_model`

    :param mass_info: Dictionary to fill
    :type mass_info: dict
    :param ctx: Data that defines the synthetic clusters
    :type ctx: ModelEngine
    :param p: Values of the fundamental parameters, see :py:func:`get_model`
    :type p: np.ndarray
    :param N_synth_stars: Number of synthetic stars requested
    :type N_synth_stars: int
    :param isoch_moved: Isochrone moved by the distance modulus
    :type isoch_moved: np.ndarray
    :param idx_cut: Indexes of the isochrone points brighter than the maximum
        magnitude
    :type idx_cut: np.ndarray
    :param N_pool: Final size of the pool of candidate stars
    :type N_pool: int
    :param keep: Indexes of the candidates brighter than the maximum magnitude
    :type keep: np.ndarray
    """
    m_ini_idx = ctx.m_ini_idx
    met, loga, alpha, beta = p[:4].tolist()
    # Sampled IMF masses used (see get_model())
    ml, _, al, _ = grid_box(ctx.met_age_dict, met, loga)
    st_mass = ctx.st_dist_mass[ml][al]
    binar_flag = not (alpha == 0.0 and beta == 0.0)
    # Indexes in the IMF sample of the stars in the pool (the same selection done by
    # mass_interp()), and of those kept in the synthetic cluster
    mass_row = isoch_moved[m_ini_idx]
    pool_idx = first_in_range(
        st_mass, mass_row[idx_cut[0]], mass_row[idx_cut[-1]], N_pool
    )
    kept_idx = pool_idx[keep[:N_synth_stars]]
    # The simulated cluster is made of the IMF masses up to the last kept star. If
    # fewer than 'N_synth_stars' stars were kept, the whole IMF sample was used (see
    # get_model())
    N_end = kept_idx[-1] + 1 if keep.size >= N_synth_stars else st_mass.size
    mass_info.update(
        {
            "imf_masses": st_mass[:N_end],
            "kept_idx": kept_idx,
            "N_target": N_synth_stars,
            "mass_row": mass_row.copy(),
            "m2_row": isoch_moved[m_ini_idx + 1].copy() if binar_flag else None,
            "alpha": alpha,
            "beta": beta,
        }
    )


def per_star_pool(N_synth_stars: int, dr: float) -> int:
    """Initial number of candidate stars generated by :py:func:`get_model`.

    Each synthetic star gets its own absorption ``A = Av + f * DR``, applied after
    the interpolation of the IMF masses into the isochrone. The stars are selected
    with the magnitude of the isochrone at the minimum absorption ``Av``, so with
    ``DR > 0`` some of them end up below the maximum magnitude and are discarded.
    50% more candidates (plus 10) are generated in that case, and 5% more without
    differential reddening (where only stars interpolated close to the magnitude
    limit can be lost).

    :param N_synth_stars: Number of synthetic stars required
    :type N_synth_stars: int
    :param dr: Differential reddening
    :type dr: float

    :return: Number of candidates
    :rtype: int
    """
    if dr > 0.0:
        return N_synth_stars + N_synth_stars // 2 + 10
    return N_synth_stars + N_synth_stars // 20 + 10


def per_star_pool_next(
    N_pool: int, N_kept: int, N_stars: int, N_synth_stars: int, N_max: int
) -> int:
    """Size of the next pool of candidates in :py:func:`get_model`, when the
    current one gave too few stars above the maximum magnitude: the size required
    by the observed fraction of surviving stars, plus a 20% margin (at least 1.5x
    the current size), and at most the size of the IMF sample.

    :param N_pool: Current number of candidates
    :type N_pool: int
    :param N_kept: Number of candidates that survived the magnitude cut
    :type N_kept: int
    :param N_stars: Number of candidates generated (``<= N_pool``)
    :type N_stars: int
    :param N_synth_stars: Number of synthetic stars required
    :type N_synth_stars: int
    :param N_max: Size of the IMF sample
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


def zaWAverage(
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
    closest_iso = theor_tracks[m0][a0]
    isochrone[m_ini_idx] = closest_iso[m_ini_idx]
    # Now for the secondary masses
    isochrone[m_ini_idx + 1] = closest_iso[m_ini_idx + 1]

    return isochrone, m0, a0


if _HAS_NUMBA:

    @njit(cache=True)
    def _zaWAverage_nb(
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
        :py:func:`zaWAverage`. A single pass reads the four grid isochrones and
        writes the weighted sum, with the same operations in the same order as
        the NumPy version (identical results), without temporary arrays. The
        two mass rows are copied from the closest node (m0, a0) instead of being
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
            if r == m_ini_idx or r == m_ini_idx + 1:
                # DO NOT average the masses (see zaWAverage())
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


def move_isochrone(
    isochrone: np.ndarray, binar_flag: bool, m_ini_idx: int, dm: float
) -> np.ndarray:
    """Receives an isochrone of a given age and metallicity and modifies
    its magnitude values according to a given distance modulus.

    :param isochrone: Isochrone array.
    :type isochrone: np.ndarray
    :param binar_flag: Binary flag
    :type binar_flag: bool
    :param m_ini_idx: Index of the initial mass.
    :type m_ini_idx: int
    :param dm: Distance modulus.
    :type dm: float

    :returns: Modified isochrone array.
    :rtype: np.ndarray
    """
    isochrone = np.ascontiguousarray(isochrone)
    # Move magnitude
    isochrone[0] += dm
    # Move binary magnitude if required
    if binar_flag:
        # m_ini_idx + 1 --> Secondary binary mass
        # m_ini_idx + 2 --> Binary magnitude
        isochrone[m_ini_idx + 2] += dm

    return isochrone


def extinction(
    ext_law: str,
    ext_coefs: list | np.ndarray,
    dr_frac: np.ndarray,
    m_ini_idx: int,
    binar_flag: bool,
    Av: float,
    dr: float,
    Rv: float,
    isochrone: np.ndarray,
) -> np.ndarray:
    """Modifies magnitude and color(s) according to given values for the
    total absorption Av. Using this parameter instead of the E(B-V) extinction
    reduces the correlation with Rv.

    The distance modulus was applied before this function.

    isochrone = [mag, c1, (c2), Mini, Mini_b, mag_b, c1b, (c2b)]

    For the CCMO model:

    ext_coefs = [mag_ec, c1_ec, ...]

    where:
    mag_ec = [a, b]  ; cX_ec = [[a1, b1], [a2, b2]]

    and:
    ccm_coef = a + b / Rv = ext_coefs[i][0] + ext_coefs[i][1] / Rv

    With differential reddening (dr>0) the total absorption of each column (a
    synthetic star, see :py:func:`get_model`) is:

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
    :param ext_coefs: List of extinction coefficients.
    :type ext_coefs: list | np.ndarray
    :param dr_frac: Fractions of the differential reddening applied to each column.
    :type dr_frac: np.ndarray
    :param m_ini_idx: Index of the initial mass.
    :type m_ini_idx: int
    :param binar_flag: Flag to indicate if binarity is being used.
    :type binar_flag: bool
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
        # Compiled version: identical results, see _extinction_gaiadr3_nb()
        # dr_kind: 0 -> no differential reddening (scalar A), 1 -> A = Av + f*dr
        dr_kind = 1 if dr > 0.0 else 0
        _extinction_gaiadr3_nb(
            isochrone, dr_frac, dr_kind, Av, Av, dr, m_ini_idx, binar_flag
        )
        return isochrone

    if _HAS_NUMBA and ext_law == "CCMO":
        # Compiled version: identical results, see _extinction_ccmo_nb()
        ec_mag, ec_col1, ec_col2 = ccmo_coefs(ext_coefs, Rv)
        _extinction_ccmo_nb(
            isochrone,
            dr_frac,
            1 if dr > 0.0 else 0,
            Av,
            dr,
            m_ini_idx,
            binar_flag,
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
        ec_mag, ec_col1, ec_col2 = ccmo_coefs(ext_coefs, Rv)
        Ax = ec_mag * Av_dr
        Ex1 = ec_col1 * Av_dr
    elif ext_law == "GAIADR3":
        # If this model is used the first color is always expected to be BP-RP
        # BP_RP = isochrone[1]
        Ax, Ex1 = gaiadr3_ext(isochrone[1], Av_dr)
    else:
        raise ValueError(f"Unknown extinction law: {ext_law}")

    isochrone[0] += Ax  # Magnitude
    isochrone[1] += Ex1  # First color

    # Move binary data.
    if binar_flag:
        isochrone[m_ini_idx + 2] += Ax  # Magnitude (binary system)
        isochrone[m_ini_idx + 3] += Ex1  # First color (binary system)

    # Second color (CCMO only)
    if len(ext_coefs) > 2:
        Ex2 = ec_col2 * Av_dr
        isochrone[2] += Ex2
        # Move color with binary data.
        if binar_flag:
            isochrone[m_ini_idx + 4] += Ex2  # Second color (binary system)

    return isochrone


if _HAS_NUMBA:

    @njit(cache=True)
    def _extinction_gaiadr3_nb(
        isochrone: np.ndarray,
        dr_frac: np.ndarray,
        dr_kind: int,
        A_scalar: float,
        Av: float,
        dr: float,
        m_ini_idx: int,
        binar_flag: bool,
    ) -> None:
        """Compiled (numba) version of the GAIADR3 branch of
        :py:func:`extinction`, modifying 'isochrone' in place.

        A single loop over the isochrone points computes the (differentially
        reddened) absorption, A_G and E(BP-RP) (same nested form as
        :py:func:`_ext_horner_XA` and :py:func:`_ext_horner_X`), and adds them to
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

        i_mb = m_ini_idx + 2
        i_cb = m_ini_idx + 3
        for i in range(isochrone.shape[1]):
            # Clamped to the range of validity of the coefficients, as in
            # gaiadr3_ext()
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
            if binar_flag:
                isochrone[i_mb, i] += Ax  # Magnitude (binary system)
                isochrone[i_cb, i] += Ex  # First color (binary system)

    @njit(cache=True)
    def _extinction_ccmo_nb(
        isochrone: np.ndarray,
        dr_frac: np.ndarray,
        dr_kind: int,
        Av: float,
        dr: float,
        m_ini_idx: int,
        binar_flag: bool,
        ec_mag: float,
        ec_col1: float,
        ec_col2: float,
        col2_flag: bool,
    ) -> None:
        """Compiled (numba) version of the CCMO branch of :py:func:`extinction`,
        modifying 'isochrone' in place. The coefficients ``ec_*`` (``a + b / Rv``
        values, see :py:func:`ccmo_coefs`) multiply the (differentially reddened)
        absorption of each point. The floating point operations are the same and in
        the same order as in the NumPy version, so the results are identical.

        dr_kind: 0 -> scalar absorption Av (no differential reddening),
        1 -> differential reddening A = Av + f * dr, with f from 'dr_frac'.
        """
        i_mb = m_ini_idx + 2
        i_cb = m_ini_idx + 3
        i_c2b = m_ini_idx + 4
        A = Av
        for i in range(isochrone.shape[1]):
            if dr_kind == 1:
                A = dr_frac[i] * dr
                A += Av
            Ax = ec_mag * A
            Ex1 = ec_col1 * A
            isochrone[0, i] += Ax  # Magnitude
            isochrone[1, i] += Ex1  # First color
            if binar_flag:
                isochrone[i_mb, i] += Ax  # Magnitude (binary system)
                isochrone[i_cb, i] += Ex1  # First color (binary system)
            if col2_flag:
                Ex2 = ec_col2 * A
                isochrone[2, i] += Ex2  # Second color
                if binar_flag:
                    isochrone[i_c2b, i] += Ex2  # Second color (binary system)


def ccmo_coefs(ext_coefs: list | np.ndarray, Rv: float) -> tuple[float, float, float]:
    """Coefficients of the CCMO extinction law for a given ``Rv``: the ratios
    ``A_mag / Av`` and ``E(color) / Av`` for the magnitude and the color(s), see
    :py:func:`extinction`.

    :param ext_coefs: List of extinction coefficients, ``[ec_mag, ec_col1, ...]``
    :type ext_coefs: list | np.ndarray
    :param Rv: Total-to-selective extinction ratio.
    :type Rv: float

    :returns: Coefficients for the magnitude, the first color and the second color
        (``0`` if there is no second color).
    :rtype: tuple[float, float, float]
    """
    ec_mag = ext_coefs[0][0] + ext_coefs[0][1] / Rv
    ec_col1 = (ext_coefs[1][0][0] + ext_coefs[1][0][1] / Rv) - (
        ext_coefs[1][1][0] + ext_coefs[1][1][1] / Rv
    )
    ec_col2 = 0.0
    if len(ext_coefs) > 2:
        ec_col2 = (ext_coefs[2][0][0] + ext_coefs[2][0][1] / Rv) - (
            ext_coefs[2][1][0] + ext_coefs[2][1][1] / Rv
        )
    return float(ec_mag), float(ec_col1), float(ec_col2)


def gaiadr3_ext(X: np.ndarray, A: float | np.ndarray) -> tuple[np.ndarray, np.ndarray]:
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
        return _ext_horner_XA(_G_COEFFS, X, A), _ext_horner_XA(_BPRP_COEFFS, X, A)
    return _ext_horner_X(_G_COEFFS, X, A), _ext_horner_X(_BPRP_COEFFS, X, A)


def _ext_horner_XA(c: tuple, X: np.ndarray, A: np.ndarray) -> np.ndarray:
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


def _ext_horner_X(c: tuple, X: np.ndarray, A: float) -> np.ndarray:
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


def cut_max_mag(
    isoch_moved: np.ndarray, max_mag_syn: float
) -> tuple[np.ndarray, np.ndarray]:
    """Identify the isochrone points with magnitude values smaller than the maximum
    observed value.

    The isochrone is not cut here: ``isoch_moved[:, idx]`` would copy the entire
    array, which is not needed since :py:func:`mass_interp` only reads the
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


def mass_interp(
    isochrone: np.ndarray,
    idx_cut: np.ndarray,
    msk_cut: np.ndarray,
    m_ini_idx: int,
    st_dist_mass: np.ndarray,
    N_synth_stars: int,
    mass_idx: np.ndarray,
    # binar_flag: bool,
) -> np.ndarray:
    """For each mass in the sampled IMF mass distribution, interpolate its value
    (and those of all the sub-arrays in the isochrone) into the isochrone.

    Only the isochrone points given by 'idx_cut' (those below the maximum
    magnitude, see :py:func:`cut_max_mag`) are used. The result is identical to
    interpolating into the cut isochrone ``isochrone[:, idx_cut]``, without
    having to copy it.

    Masses that fall outside of the (cut) isochrone's mass range are rejected.

    :param isochrone: Isochrone array (not cut).
    :type isochrone: np.ndarray
    :param idx_cut: Sorted indexes of the isochrone points kept by
        :py:func:`cut_max_mag`.
    :type idx_cut: np.ndarray
    :param msk_cut: Boolean mask of the isochrone points kept by
        :py:func:`cut_max_mag`.
    :type msk_cut: np.ndarray
    :param m_ini_idx: Index of the initial mass.
    :type m_ini_idx: int
    :param st_dist_mass: Array of sampled masses.
    :type st_dist_mass: np.ndarray
    :param N_synth_stars: Number of observed stars.
    :type N_synth_stars: int
    :param mass_idx: Pre-computed search indexes of 'st_dist_mass' into the full
        (before the maximum magnitude cut) mass grid of the isochrone, see
        :py:func:`precompute.mass_idx_lookup`.
    :type mass_idx: np.ndarray

    :returns: Interpolated isochrone array.
    :rtype: np.ndarray
    """
    # ==IMPORTANT==
    # Assumes `mass_ini=isochrone[m_ini_idx]` is sorted min to max. This is checked
    # when loading the isochrones in isochrone_loader.interp_isochrones()
    # ==IMPORTANT==
    mass_ini = isochrone[m_ini_idx]

    # Extract up to 'N_synth_stars' masses sampled from an IMF, within the mass range
    # of the cut isochrone (the first ones, as np.nonzero(msk)[0][:N_synth_stars]
    # with 'msk' the masses in that range)
    idx = first_in_range(
        st_dist_mass, mass_ini[idx_cut[0]], mass_ini[idx_cut[-1]], N_synth_stars
    )

    mass_dist = st_dist_mass[idx]

    if mass_dist.size == 0:
        return np.array([])

    mass_dist = np.ascontiguousarray(mass_dist, dtype=mass_ini.dtype)

    # Indexes into the *cut* mass grid, i.e.:
    # np.searchsorted(mass_ini[idx_cut], mass_dist)
    # Obtained by converting the pre-computed indexes into the full mass grid, into
    # indexes into the cut grid: the number of grid points with a mass smaller than
    # the sampled mass that survived the cut
    N_kept = np.empty(msk_cut.size + 1, dtype=np.intp)
    N_kept[0] = 0
    np.cumsum(msk_cut, out=N_kept[1:])
    hi_cut = N_kept[mass_idx[idx]]

    # Map back to indexes into the full isochrone. Note that hi_cut=0 gives
    # lo=idx_cut[-1], exactly as indexing the cut isochrone with -1 would
    x_new_indices = idx_cut[hi_cut]
    lo = idx_cut[hi_cut - 1]

    # Interpolate the sampled stars (masses) into the isochrone
    isoch_mass = interp_mass_isoch(isochrone, mass_ini, mass_dist, x_new_indices, lo)

    return isoch_mass


if _HAS_NUMBA:

    @njit(cache=True)
    def _mass_interp_nb(
        isochrone: np.ndarray,
        idx_cut: np.ndarray,
        msk_cut: np.ndarray,
        m_ini_idx: int,
        st_dist_mass: np.ndarray,
        N_synth_stars: int,
        mass_idx: np.ndarray,
    ) -> np.ndarray:
        """Compiled (numba) version of :py:func:`mass_interp` (with the
        pre-computed 'mass_idx' search indexes). It performs the same floating
        point operations in the same order, so the results are identical to the
        NumPy version, but it avoids most temporary arrays and Python calls:

        1. Select the first 'N_synth_stars' sampled masses within the mass range
           of the cut isochrone, stopping as soon as they are found. The loop is
           written without branches, which is faster.
        2. Bracketing points: for each point 'k' of the full mass grid, the next
           kept point at or after 'k' and the last kept point before 'k' are
           tabulated once. With 'mass_idx' this replaces the cumulative count
           used in the NumPy version.
        3. Interpolate the isochrone linearly, one row at a time (rows are
           contiguous in memory).
        """
        mass_ini = isochrone[m_ini_idx]
        N_cut = idx_cut.size
        v_min = mass_ini[idx_cut[0]]
        v_max = mass_ini[idx_cut[N_cut - 1]]

        # 1. First 'N_synth_stars' masses within [v_min, v_max]. The index is
        # always written and the counter only advances for masses in range.
        sel = np.empty(st_dist_mass.size, dtype=np.intp)
        N = 0
        for i in range(st_dist_mass.size):
            m = st_dist_mass[i]
            sel[N] = i
            N += (m >= v_min) & (m <= v_max)
            if N >= N_synth_stars:
                break
        N = min(N, N_synth_stars)

        # 2. Next kept point at or after k (nxt), and last kept point before k
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

        hi = np.empty(N, dtype=np.intp)
        lo = np.empty(N, dtype=np.intp)
        w = np.empty(N, dtype=isochrone.dtype)
        for j in range(N):
            i = sel[j]
            k = mass_idx[i]
            hi[j] = nxt[k]
            lo[j] = prv[k]
            x_lo = mass_ini[lo[j]]
            w_j = st_dist_mass[i] - x_lo
            w_j /= mass_ini[hi[j]] - x_lo
            w[j] = w_j

        # 3. Linear interpolation, y_lo + (y_hi - y_lo) * w
        N_rows = isochrone.shape[0]
        isoch_mass = np.empty((N_rows, N), dtype=isochrone.dtype)
        for r in range(N_rows):
            row = isochrone[r]
            out = isoch_mass[r]
            for j in range(N):
                y_lo = row[lo[j]]
                out[j] = (row[hi[j]] - y_lo) * w[j] + y_lo

        return isoch_mass


def first_in_range(
    arr: np.ndarray, v_min: float, v_max: float, N_max: int
) -> np.ndarray:
    """Return the indexes of the first (at most) 'N_max' elements of 'arr' that
    fall within the closed range [v_min, v_max].

    Equivalent to ``np.nonzero((arr >= v_min) & (arr <= v_max))[0][:N_max]`` but
    the mask is evaluated over progressively larger chunks of 'arr', stopping as
    soon as 'N_max' elements are found. Since the IMF sampled array is usually
    much larger than the number of required stars (e.g.: ~55000 vs ~1500), this
    avoids masking the entire array in most cases.

    :param arr: Array of values.
    :type arr: np.ndarray
    :param v_min: Minimum value of the range.
    :type v_min: float
    :param v_max: Maximum value of the range.
    :type v_max: float
    :param N_max: Maximum number of indexes to return.
    :type N_max: int

    :returns: Array of indexes.
    :rtype: np.ndarray
    """
    N_tot = arr.size
    # The initial chunk size of 4*N_max was found to be a good compromise for
    # typical in-range fractions; the chunk size is doubled in each iteration
    i_end = min(N_tot, 4 * N_max)
    i_start, N_found, idx_lst = 0, 0, []
    while True:
        chunk = arr[i_start:i_end]
        idx = np.flatnonzero((chunk >= v_min) & (chunk <= v_max))
        if i_start > 0:
            idx += i_start
        idx_lst.append(idx)
        N_found += idx.size
        if N_found >= N_max or i_end == N_tot:
            break
        i_start, i_end = i_end, min(N_tot, 2 * i_end)

    if len(idx_lst) == 1:
        return idx[:N_max]
    return np.concatenate(idx_lst)[:N_max]


def interp_mass_isoch(
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
        (see :py:func:`mass_interp`).
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


def binarity(
    alpha: float,
    beta: float,
    binar_flag: bool,
    m_ini_idx: int,
    rand_unif_vals: np.ndarray,
    isoch_mass: np.ndarray,
) -> np.ndarray:
    """Select a fraction of stars to be binaries, given a chosen method.

    :param alpha: Alpha parameter for the binary fraction.
    :type alpha: float
    :param beta: Beta parameter for the binary fraction.
    :type beta: float
    :param binar_flag: Flag to indicate if binarity is being used.
    :type binar_flag: bool
    :param m_ini_idx: Index of the initial mass.
    :type m_ini_idx: int
    :param rand_unif_vals: Array of random uniform values.
    :type rand_unif_vals: np.ndarray
    :param isoch_mass: Isochrone array.
    :type isoch_mass: np.ndarray

    :returns: Isochrone array with binary stars.
    :rtype: np.ndarray
    """
    # No binarity process defined
    if binar_flag is False:
        # Update the binary systems' masses so that the secondary masses for
        # SINGLE systems are identified with a 'nan' value.
        isoch_mass[m_ini_idx + 1] = np.nan
        return isoch_mass[: m_ini_idx + 2]

    Ns = isoch_mass.shape[-1]

    # Distribution of probability of binarity (multiplicity fraction) versus
    # primary masses.

    # Offner et al. (2022); Fig 1 (left), Table 1
    # b_p = np.clip(alpha + beta / (1.0 + 1.4 / isoch_mass[m_ini_idx]), 0.0, 1.0)
    # Same operations (hence identical values) as the line above, but in place.
    # Clipping to [0, 1] is not needed: the probabilities are only compared below
    # with uniform random values in [0, 1), so a value < 0 (> 1) gives the same
    # result as 0 (1)
    b_p = 1.4 / isoch_mass[m_ini_idx]
    b_p += 1.0
    np.divide(beta, b_p, out=b_p)
    b_p += alpha

    # Stars (masses) with the largest binary probabilities are selected
    # proportional to their probability
    bin_indxs = b_p > rand_unif_vals[:Ns]

    # Index of the binary magnitude: mag_binar = m_ini_idx + 1
    # Update array with new values of magnitudes and colors for the binary
    # systems. This does not change the primary mass values, just the magnitude and
    # color(s).
    for i in range(m_ini_idx):
        # isoch_mass[i][bin_indxs] = isoch_mass[m_ini_idx + 2 + i][bin_indxs]
        # Use np.copyto with where parameter (faster than boolean indexing)
        np.copyto(isoch_mass[i], isoch_mass[m_ini_idx + 2 + i], where=bin_indxs)

    # Update the binary systems' masses so that the secondary masses for
    # SINGLE systems are identified with a 'nan' value. Using np.copyto() on the
    # row is faster than the 2D boolean indexing isoch_mass[m_ini_idx + 1, ~bin_indxs]
    np.copyto(isoch_mass[m_ini_idx + 1], np.nan, where=~bin_indxs)

    # Return [mag, c1, (c2), mass, mass_b]
    return isoch_mass[: m_ini_idx + 2]


if _HAS_NUMBA:

    @njit(cache=True)
    def _scatter_add_nb(row: np.ndarray, order: np.ndarray, sigma: np.ndarray) -> None:
        """Compiled equivalent of ``row[order] += sigma[:order.size]``, used by
        :py:func:`add_errors`"""
        for j in range(order.size):
            row[order[j]] += sigma[j]


if _HAS_NUMBA:

    @njit(cache=True)
    def _binarity_nb(
        alpha: float,
        beta: float,
        m_ini_idx: int,
        rand_unif_vals: np.ndarray,
        isoch_mass: np.ndarray,
    ) -> None:
        """Compiled (numba) version of :py:func:`binarity` (with binaries),
        modifying 'isoch_mass' in place. A single loop over the stars computes
        the binary probability (same operations, in the same order, as the
        NumPy version: identical results), compares it with the uniform random
        value, and either replaces the magnitude and color(s) with those of the
        binary system, or flags the secondary mass of a single star with 'nan'.
        """
        mass = isoch_mass[m_ini_idx]
        mass_b = isoch_mass[m_ini_idx + 1]
        for j in range(isoch_mass.shape[1]):
            # alpha + beta / (1.0 + 1.4 / mass)
            b_p = 1.4 / mass[j]
            b_p += 1.0
            b_p = beta / b_p
            b_p += alpha
            if b_p > rand_unif_vals[j]:
                for i in range(m_ini_idx):
                    isoch_mass[i, j] = isoch_mass[m_ini_idx + 2 + i, j]
            else:
                mass_b[j] = np.nan


def add_errors(isoch_binar: np.ndarray, err_dist: list[np.ndarray]) -> np.ndarray:
    """Add random synthetic uncertainties to the magnitude and color(s)

    :param isoch_binar: Isochrone array.
    :type isoch_binar: np.ndarray
    :param err_dist: List of arrays of error distributions.
    :type err_dist: list[np.ndarray]

    :returns: Isochrone array with added errors.
    :rtype: np.ndarray
    """
    if not err_dist:
        return isoch_binar

    N = isoch_binar.shape[1]
    # Indexes that sort the stars from the largest to the smallest magnitude, the
    # same order used by the error distributions (see error_distribution())
    mag_sort = np.argsort(-isoch_binar[0])

    for i, sigma in enumerate(err_dist):
        # Index the row first: fancy indexing a 1D view is noticeably faster
        # than the mixed (row, indexes) indexing of the 2D array
        row = isoch_binar[i]
        row[mag_sort] += sigma[:N]

    return isoch_binar

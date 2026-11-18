import warnings
from collections.abc import Sequence
from typing import TYPE_CHECKING

import numpy as np
from scipy.spatial import KDTree

from . import HAS_NUMBA as _HAS_NUMBA
from . import njit

if TYPE_CHECKING:
    from ..isochrones import Isochrones


def turn_off_points(
    isoch_tracks: list[np.ndarray],
    max_masses: list[float],
    color_idx: int,
    color_rise: float = 0.3,
) -> dict:
    """Turn-off points of the isochrones of the sampled synthetic clusters, used to
    define the BSS region.

    Each isochrone is cut at the maximum mass of its synthetic cluster and its points
    are ordered by initial mass. Starting from the lowest mass, the main sequence
    gets bluer with increasing mass until the turn-off, after which the color
    increases (subgiant and giant branches). The turn-off point is the bluest point
    before the color increases by more than ``color_rise``. If it never does (young
    isochrones), it is the bluest point of the whole isochrone.

    :param isoch_tracks: Magnitude, color and initial mass of each isochrone, arrays
        of shape ``(3, N_points)``
    :type isoch_tracks: list[np.ndarray]
    :param max_masses: Maximum initial mass of each synthetic cluster
    :type max_masses: list[float]
    :param color_idx: Index of the color used, stored in the returned dictionary
    :type color_idx: int
    :param color_rise: Color increase (mag) along the isochrone, after its bluest
        point, that marks the turn-off point (avoids stopping at small wiggles of the
        main sequence, such as the knee of the overshooting hook)
    :type color_rise: float

    :return: Dictionary with the keys ``color_idx, to_col, to_mag, manual``
        (``manual=False``). The turn-off point is ``(nan, nan)`` for an isochrone
        with no finite points below the maximum mass
    :rtype: dict
    """
    turn_off = {"color_idx": color_idx, "to_col": [], "to_mag": [], "manual": False}
    for track, max_mass in zip(isoch_tracks, max_masses):
        # Isochrone up to the maximum mass of the synthetic cluster, with finite
        # photometry, ordered by initial mass
        msk = (track[2] < max_mass) & np.isfinite(track[0]) & np.isfinite(track[1])
        mag, col, mass = track[:, msk]
        order = np.argsort(mass, kind="stable")
        mag, col = mag[order], col[order]

        to_col, to_mag = np.nan, np.nan
        if col.size > 0:
            # First point redder than the bluest previous point by more than
            # 'color_rise'
            rise = np.flatnonzero(col > np.minimum.accumulate(col) + color_rise)
            end = rise[0] if rise.size else col.size
            i_to = int(np.argmin(col[:end]))
            to_col, to_mag = float(col[i_to]), float(mag[i_to])

        turn_off["to_col"].append(to_col)
        turn_off["to_mag"].append(to_mag)

    return turn_off


def turn_off_points_manual(
    rng: np.random.Generator,
    turn_off: Sequence[float],
    N_points: int,
    color_idx: int,
) -> dict:
    """Turn-off points sampled from a normal distribution around a turn-off point
    given by the user, used to define the BSS region instead of the points estimated
    by :py:func:`turn_off_points`

    :param rng: Random number generator
    :type rng: np.random.Generator
    :param turn_off: Turn-off point and its STDDEVs, ``(color, magnitude,
        stddev_color, stddev_magnitude)``
    :type turn_off: Sequence[float]
    :param N_points: Number of points to sample (one per synthetic cluster)
    :type N_points: int
    :param color_idx: Index of the color used, stored in the returned dictionary
    :type color_idx: int

    :return: Dictionary with the keys ``color_idx, to_col, to_mag, manual``
        (``manual=True``)
    :rtype: dict
    """
    col, mag, col_std, mag_std = turn_off
    return {
        "color_idx": color_idx,
        "to_col": list(rng.normal(col, col_std, N_points)),
        "to_mag": list(rng.normal(mag, mag_std, N_points)),
        "manual": True,
    }


def stellar_parameters(
    cluster_mag: np.ndarray,
    cluster_e_mag: np.ndarray,
    cluster_colors: list,
    cluster_e_colors: list,
    m_ini_idx: int,
    synthcls: list[np.ndarray],
    turn_off_points: dict,
    isochs_young: "Isochrones | None",
    extrinsic_pars: dict,
    ext_law: str,
    ext_coefs: list | np.ndarray,
    rng: np.random.Generator,
) -> tuple[dict, np.ndarray, np.ndarray]:
    """Per-star masses and probabilities of the observed stars, used by
    :py:meth:`asteca.Synthetic.cluster_properties`

    :param cluster_mag: Magnitudes of the observed stars
    :type cluster_mag: np.ndarray
    :param cluster_e_mag: Magnitude uncertainties of the observed stars
    :type cluster_e_mag: np.ndarray
    :param cluster_colors: Colors of the observed stars, ``[color, color2]``
    :type cluster_colors: list
    :param cluster_e_colors: Color uncertainties of the observed stars
    :type cluster_e_colors: list
    :param m_ini_idx: Index of the initial mass in the synthetic clusters
    :type m_ini_idx: int
    :param synthcls: Synthetic clusters of the sampled models
    :type synthcls: list[np.ndarray]
    :param turn_off_points: Turn-off points of the synthetic clusters, see
        :py:func:`turn_off_points` and :py:func:`turn_off_points_manual`
    :type turn_off_points: dict
    :param isochs_young: Isochrones with ``log(age)<=7``, or ``None``
    :type isochs_young: Isochrones | None
    :param extrinsic_pars: Values and STDDEVs of the extrinsic parameters, and the
        valid range of each parameter (``values, stddevs, limits``)
    :type extrinsic_pars: dict
    :param ext_law: Extinction law
    :type ext_law: str
    :param ext_coefs: Extinction coefficients
    :type ext_coefs: list | np.ndarray
    :param rng: Random number generator
    :type rng: np.random.Generator

    :return: Dictionary of per-star arrays, the secondary masses assigned to the
        observed stars by each synthetic cluster (``nan`` for single systems), and
        the mask of the stars with valid photometry
    :rtype: tuple[dict, np.ndarray, np.ndarray]
    """
    # Stellar masses
    m1_med, m1_std, m2_med, m2_std, nan_msk, m2_obs = get_stellar_masses(
        cluster_mag,
        cluster_e_mag,
        cluster_colors,
        cluster_e_colors,
        m_ini_idx,
        synthcls,
    )

    # Binary probabilities
    binar_probs = get_binar_probabilities(m2_obs)

    # BSS probabilities
    bss_probs = get_bss_probabilities(
        cluster_mag, cluster_colors, turn_off_points, len(synthcls)
    )

    binar_probs = (1 - bss_probs) * binar_probs
    single_probs = 1 - binar_probs - bss_probs

    N_cluster = len(cluster_mag)
    bss_mass_median, bss_mass_stddev = (
        np.full(N_cluster, np.nan),
        np.full(N_cluster, np.nan),
    )
    if isochs_young is not None:
        bss_mass_median, bss_mass_stddev = get_bss_masses(
            cluster_mag,
            cluster_colors[0],
            extrinsic_pars,
            ext_law,
            ext_coefs,
            rng,
            isochs_young.met_age_dict,
            isochs_young.theor_tracks,
        )
    else:
        warnings.warn(
            "No young isochrones were provided, BSSs masses can not be estimated"
        )

    data = {
        "mass_single": m1_med,
        "mass_single_std": m1_std,
        "mass_binar": m2_med,
        "mass_binar_std": m2_std,
        "mass_bss": bss_mass_median,
        "mass_bss_std": bss_mass_stddev,
        "single_prob": single_probs,
        "binar_prob": binar_probs,
        "bss_prob": bss_probs,
    }

    # Assign all nans to stars with a photometric nan in any dimension
    for value in data.values():
        value[nan_msk] = np.nan

    return data, m2_obs, ~nan_msk


def get_stellar_masses(
    cluster_mag: np.ndarray,
    cluster_e_mag: np.ndarray,
    cluster_colors: list,
    cluster_e_colors: list,
    m_ini_idx: int,
    sampled_synthcls: list,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Assign primary and secondary synthetic masses to each observed star,
    for each sampled synthetic model.

    Photometric dimensions are normalized by representative observational
    uncertainties before performing the nearest-neighbour search.

    :param cluster_mag: Observed magnitude.
    :type cluster_mag: np.ndarray
    :param cluster_e_mag: Magnitude uncertainties.
    :type cluster_e_mag: np.ndarray
    :param cluster_colors: Observed colors.
    :type cluster_colors: list
    :param cluster_e_colors: Color uncertainties.
    :type cluster_e_colors: list
    :param m_ini_idx: Index of the initial mass column.
    :type m_ini_idx: int
    :param sampled_synthcls: Sampled synthetic cluster data.
    :type sampled_synthcls: list

    :return: Primary and secondary mass values (median + stddev), a mask for
        NaN values, and secondary masses assigned across all sampled models.
    :rtype: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]
    """

    # Extract observed photometry
    cl_colors = [cluster_colors[0]]
    cl_e_colors = [cluster_e_colors[0]]
    if cluster_colors[1] is not None:
        cl_colors.append(cluster_colors[1])
        cl_e_colors.append(cluster_e_colors[1])

    obs_phot = np.array([cluster_mag, *cl_colors])
    obs_e_phot = np.array([cluster_e_mag, *cl_e_colors])

    # Identify stars with invalid photometry or uncertainties
    nan_msk = np.any(~np.isfinite(obs_phot), axis=0)
    nan_msk |= np.any(
        (~np.isfinite(obs_e_phot)) | (obs_e_phot <= 0.0),
        axis=0,
    )

    # Representative uncertainty for each photometric dimension.
    # Use the median to avoid a few poorly measured stars dominating.
    valid_msk = ~nan_msk
    phot_scales = np.nanmedian(
        np.where(valid_msk, obs_e_phot, np.nan),
        axis=1,
    )
    # Protect against pathological cases.
    phot_scales[~np.isfinite(phot_scales) | (phot_scales <= 0.0)] = 1.0
    # Normalize observed photometry
    obs_phot = obs_phot / phot_scales[:, None]

    # Only query stars with valid photometry
    obs_phot_valid = obs_phot[:, valid_msk].T

    # Assign primary and secondary masses for each sampled model
    m12_obs = []
    for isoch in sampled_synthcls:
        # Synthetic photometry
        isoch_phot = isoch[:m_ini_idx] / phot_scales[:, None]

        # Find photometrically closest synthetic stars
        tree = KDTree(isoch_phot.T)
        _, close_stars_idxs = tree.query(obs_phot_valid, k=1)

        # Primary and secondary masses
        masses = isoch[m_ini_idx:, close_stars_idxs]

        # Restore the full observed-star array
        masses_full = np.full(
            (masses.shape[0], obs_phot.shape[1]),
            np.nan,
        )
        masses_full[:, valid_msk] = masses

        m12_obs.append(masses_full)

    # Shape:
    # (N_sampled_models, 2, N_observed_stars)
    m12_obs = np.asarray(m12_obs)

    # Extract primary and secondary masses
    m1_obs = m12_obs[:, 0, :]
    m2_obs = m12_obs[:, 1, :]

    # Primary mass values
    m1_med = np.median(m1_obs, axis=0)
    m1_std = np.std(m1_obs, axis=0)

    # Secondary mass values
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m2_med = np.nanmedian(m2_obs, axis=0)
        m2_std = np.nanstd(m2_obs, axis=0)

    # Secondary mass cannot exceed primary mass
    m2_med = np.minimum(m1_med, m2_med)

    # m2_obs is used to compute the binary probability per observed star
    return m1_med, m1_std, m2_med, m2_std, nan_msk, m2_obs


def get_binar_probabilities(m2_obs: np.ndarray) -> np.ndarray:
    """Binary probability per observed star

    Count how many times the secondary mass of an observed star was assigned a
    value 'not np.nan', i.e.: was identified as a binary system. Dividing this
    value by the number of masses assigned, results in the per observed
    star probability of being a binary system.

    :param m2_obs: Secondary masses assigned to observed stars across sampled models.
    :type m2_obs: np.ndarray

    :return: Array with the binary probability for each observed star
    :rtype: np.ndarray
    """
    binar_prob = np.isfinite(m2_obs).sum(axis=0) / m2_obs.shape[0]

    return binar_prob


def poisson_binomial_percentiles(
    probs: np.ndarray, percentiles: tuple = (50, 16, 84)
) -> np.ndarray:
    """Percentiles of the number of successes in independent trials with different
    probabilities (Poisson-binomial distribution).

    The probability mass function is computed exactly by adding one trial at a
    time, ``O(N**2)`` operations for ``N`` trials. Each percentile is the smallest
    number of successes whose cumulative probability reaches it.

    :param probs: Probability of success of each trial, in ``[0, 1]``
    :type probs: np.ndarray
    :param percentiles: Percentiles to return, in ``[0, 100]``
    :type percentiles: tuple

    :return: Number of successes for each percentile
    :rtype: np.ndarray
    """
    pmf = np.zeros(probs.size + 1)
    pmf[0] = 1.0
    for i, p in enumerate(probs):
        # The right-hand side is evaluated (with the previous values) before the
        # assignment
        pmf[1 : i + 2] = pmf[1 : i + 2] * (1.0 - p) + pmf[: i + 1] * p
        pmf[0] *= 1.0 - p
    cdf = np.cumsum(pmf)
    idx = np.searchsorted(cdf, np.asarray(percentiles) / 100.0)
    # Guard against a total probability slightly below 1 from rounding
    return np.minimum(idx, probs.size).astype(float)


def get_bss_probabilities(
    cluster_mag: np.ndarray,
    cluster_colors: list[np.ndarray | None],
    turn_off_points: dict,
    N_sampled_synthcls: int,
    col_offset: float = -0.05,
    percentile: float = 5,
    dmag_init: float = 0.25,
    dmag_step: float = 0.05,
    min_stars: int = 5,
) -> np.ndarray:
    """Estimate the probability of each observed star being a blue straggler star
    (BSS) based on the turn-off point of the isochrones generated from the sampled
    models.

    For each turn-off point ``(to_col, to_mag)``, the BSS region contains the stars
    bluer than ``to_col + col_offset`` and brighter than ``to_mag``. The color shift
    detaches the region from the isochrone. The color limit is capped at the
    ``percentile`` of the observed colors around ``to_mag``, so that the region never
    extends into the observed main sequence. If the turn-off points were given by the
    user (``turn_off_points["manual"]``), they are used as given: no color shift and
    no cap are applied.

    :param cluster_mag: Observed cluster magnitudes.
    :type cluster_mag: np.ndarray
    :param cluster_colors: Observed cluster colors.
    :type cluster_colors: list[np.ndarray | None]
    :param turn_off_points: Turn-off points of the sampled models, see
        :py:func:`turn_off_points` and :py:func:`turn_off_points_manual`
    :type turn_off_points: dict
    :param N_sampled_synthcls: Number of sampled synthetic clusters.
    :type N_sampled_synthcls: int
    :param col_offset: Offset applied to the turn-off color to detach the BSS region
        from the isochrone
    :type col_offset: float
    :param percentile: Percentile used to estimate the turn-off color from observed
        stars.
    :type percentile: float
    :param dmag_init: Initial half-width of the magnitude interval around the
        detected turn-off magnitude.
    :type dmag_init: float
    :param dmag_step: Increment applied to the magnitude interval half-width.
    :type dmag_step: float
    :param min_stars: Minimum number of stars required in the magnitude interval.
    :type min_stars: int

    :return: Array with the BSS probability for each observed star
    :rtype: np.ndarray
    """
    to_col = np.array(turn_off_points["to_col"])
    to_mag = np.array(turn_off_points["to_mag"])
    manual = turn_off_points.get("manual", False)
    cl_color = cluster_colors[turn_off_points["color_idx"]]

    bss_probs = np.zeros_like(cluster_mag)

    # The magnitude range below can not contain more stars than those with a valid
    # magnitude. Cap the required number to avoid an infinite loop
    min_stars = min(min_stars, int(np.isfinite(cluster_mag).sum()))
    if min_stars == 0:
        # No stars with valid magnitudes, no BSS can be identified
        return bss_probs

    for to_c, to_m in zip(to_col, to_mag):
        if manual:
            col_lim = to_c
        else:
            # Define a magnitude range around the identified turn-off magnitude.
            # Iteratively increase the magnitude range until at least `min_stars` are
            # found within it.
            dmag = dmag_init
            while True:
                mag_range = (cluster_mag > (to_m - dmag)) & (
                    cluster_mag < (to_m + dmag)
                )
                if mag_range.sum() >= min_stars:
                    break
                dmag += dmag_step

            # This approach is more robust to incorrect parameter values assigned to
            # the observed cluster
            to_col_cl = np.nanpercentile(cl_color[mag_range], percentile)
            col_lim = min(to_col_cl, to_c + col_offset)

        # Assign BSS probabilities
        msk = (cl_color < col_lim) & (cluster_mag < to_m)
        bss_probs[msk] += 1

    # Normalize probability
    bss_probs /= N_sampled_synthcls

    # No nan probabilities
    bss_probs[np.isnan(bss_probs)] = 0.0

    return bss_probs


def get_bss_masses(
    cluster_mag: np.ndarray,
    cluster_color: np.ndarray,
    extrinsic_pars: dict,
    ext_law: str,
    ext_coefs: list | np.ndarray,
    rng: np.random.Generator,
    met_age_dict: dict,
    theor_tracks: np.ndarray,
    N_models: int = 100,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Identify the masses of blue straggler stars (BSS) by finding the closest
    point on every young isochrone/model realization.

    :param cluster_mag: Photometric magnitudes of the stars in the cluster.
    :type cluster_mag: np.ndarray
    :param cluster_color: Photometric colors of the stars in the cluster.
    :type cluster_color: np.ndarray
    :param extrinsic_pars: Values and STDDEVs of the extrinsic parameters, and the
        valid range of each parameter (``values, stddevs, limits``)
    :type extrinsic_pars: dict
    :param ext_law: Extinction law.
    :type ext_law: str
    :param ext_coefs: List of extinction coefficients.
    :type ext_coefs: list | np.ndarray
    :param rng: Random number generator.
    :type rng: np.random.Generator
    :param met_age_dict: Dictionary of metallicity and log(age) values
    :type met_age_dict: dict
    :param theor_tracks: Array of isochrones
    :type theor_tracks: np.ndarray
    :param N_models: Number of models to use in the estimation of BSS masses
    :type N_models: int

    :return: Median and STDDEVS of BSS mass estimated for each observed star.
    :rtype: tuple[np.ndarray, np.ndarray]
    """
    from . import synth_cluster_priv as scp
    from . import synth_model as sm

    # Extract young isochrones. Shape of 'isoch_youngs_arrs' is (Na, Nm, Nc, Ni):
    # Na: Number of ages
    # Nm: Number of metallicities
    # Nc: 3 (magnitude, color, initial mass)
    # Ni: Number of interpolated points
    isoch_youngs_arrs = []
    for al, agev in enumerate(met_age_dict["loga"]):
        if agev <= 7:
            mag_col = theor_tracks[:, al, :, :][:, :2, :]
            mini = theor_tracks[:, al, :, :][:, -1, :]
            isoch_youngs_arrs.append(
                np.concatenate((mag_col, mini[:, None, :]), axis=1)
            )

    # Dummy values
    binar_flag = False
    m_ini_idx = 99
    dr_frac = np.array([])
    # Set DR=0 to avoid disrupting the isochrone.
    dr = 0.0

    # Sample extrinsic parameter values
    Na, Nm = [int(_) for _ in np.shape(isoch_youngs_arrs)[:2]]
    N_samples = max(1, int(N_models / (Na * Nm) + 1))
    # Truncated to the valid ranges of the parameters
    sampled_models = {
        k: scp.truncnorm_sample(
            rng,
            f_val,
            extrinsic_pars["stddevs"][k],
            extrinsic_pars["limits"][k][0],
            extrinsic_pars["limits"][k][1],
            N_samples,
        )
        for k, f_val in extrinsic_pars["values"].items()
    }

    n_stars = cluster_mag.size
    valid = np.isfinite(cluster_mag) & np.isfinite(cluster_color)
    star_mag = np.ascontiguousarray(cluster_mag[valid], dtype=float)
    star_color = np.ascontiguousarray(cluster_color[valid], dtype=float)

    # One column per isochrone/model realization, valid stars only.
    masses = []

    for isoch_age in isoch_youngs_arrs:
        for isochrone in isoch_age:
            for rv, av, dm in zip(
                sampled_models["Rv"],
                sampled_models["Av"],
                sampled_models["dm"],
            ):
                # Move theoretical isochrone.
                isoch_moved = sm.move_isochrone(
                    isochrone.copy(),  # copy to avoid overwriting in move_isochrone()
                    binar_flag,
                    m_ini_idx,
                    dm,
                )
                # Only the magnitude and first color are present in the
                # (magnitude, color, initial mass) array: drop the coefficients of
                # a second color, which would otherwise be applied to the mass row
                isoch = sm.extinction(
                    ext_law,
                    ext_coefs[:2],
                    dr_frac,
                    m_ini_idx,
                    binar_flag,
                    av,
                    dr,
                    rv,
                    isoch_moved,
                )

                isoch_mag = np.asarray(isoch[0], dtype=float)
                isoch_color = np.asarray(isoch[1], dtype=float)
                isoch_mass = np.asarray(isoch[2])

                # Closest point on THIS isochrone for each star
                closest_idx = closest_isoch_point(
                    star_mag, star_color, isoch_mag, isoch_color
                )
                masses.append(isoch_mass[closest_idx])

    # (N_valid_stars, N_realizations)
    bss_masses = np.column_stack(masses)

    bss_mass_median = np.full(n_stars, np.nan)
    bss_mass_stddev = np.full(n_stars, np.nan)
    if np.isfinite(bss_masses).all():
        # Fast path (the usual case): no NaN masses. Same results as the
        # nan-aware functions below, ~3x faster
        bss_mass_median[valid] = np.median(bss_masses, axis=1)
        if bss_masses.shape[1] > 1:
            bss_mass_stddev[valid] = np.std(bss_masses, axis=1)
        return bss_mass_median, bss_mass_stddev

    # Extract medians and STDDEVs
    n_valid = np.sum(np.isfinite(bss_masses), axis=1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        med = np.nanmedian(bss_masses, axis=1)
        std = np.nanstd(bss_masses, axis=1)
    # 0 valid masses: median = NaN, stddev = NaN
    # 1 valid mass: median = that mass, stddev = NaN
    # ≥2 valid masses: median and stddev are calculated normally.
    med[n_valid == 0] = np.nan
    std[n_valid < 2] = np.nan
    bss_mass_median[valid] = med
    bss_mass_stddev[valid] = std

    return bss_mass_median, bss_mass_stddev


def closest_isoch_point(
    star_mag: np.ndarray,
    star_color: np.ndarray,
    isoch_mag: np.ndarray,
    isoch_color: np.ndarray,
) -> np.ndarray:
    """Index of the closest isochrone point (Euclidean distance in the
    magnitude-color plane) for each star.

    Equivalent to ``np.argmin(d2, axis=1)`` with ``d2`` the full
    ``(N_stars, N_points)`` matrix of squared distances, without building it. With
    numba installed the isochrone is sorted by magnitude and, for each star, the scan
    moves outwards from the star's magnitude and stops once the magnitude difference
    alone exceeds the best distance found (ties resolved to the lowest index, as
    ``np.argmin``). Otherwise a KD-tree is used.

    :param star_mag: Magnitudes of the stars (all finite).
    :type star_mag: np.ndarray
    :param star_color: Colors of the stars (all finite).
    :type star_color: np.ndarray
    :param isoch_mag: Magnitudes of the isochrone points.
    :type isoch_mag: np.ndarray
    :param isoch_color: Colors of the isochrone points.
    :type isoch_color: np.ndarray

    :return: Index of the closest isochrone point for each star.
    :rtype: np.ndarray
    """
    if _HAS_NUMBA:
        order = np.argsort(isoch_mag, kind="stable")
        return _closest_isoch_point_nb(
            star_mag,
            star_color,
            np.ascontiguousarray(isoch_mag[order]),
            np.ascontiguousarray(isoch_color[order]),
            order.astype(np.int64),
        )

    tree = KDTree(np.column_stack((isoch_mag, isoch_color)))
    _, idx = tree.query(np.column_stack((star_mag, star_color)))
    return np.asarray(idx)


if _HAS_NUMBA:

    @njit(cache=True)
    def _closest_isoch_point_nb(star_mag, star_color, iso_mag_s, iso_col_s, iso_idx_s):
        """Pruned nearest-point search. ``iso_mag_s`` is sorted in increasing
        order and ``iso_idx_s`` holds the original index of each sorted point.
        Squared distances are computed exactly as in the NumPy version, so the
        selected points are identical."""
        n = star_mag.size
        Ni = iso_mag_s.size
        out = np.empty(n, np.int64)
        for i in range(n):
            m = star_mag[i]
            c = star_color[i]
            j = np.searchsorted(iso_mag_s, m)
            best = np.inf
            b = Ni
            # Brighter-to-fainter half: |m - mag| grows with k
            k = j
            while k < Ni:
                dm2 = (m - iso_mag_s[k]) ** 2
                if dm2 > best:
                    break
                d2 = dm2 + (c - iso_col_s[k]) ** 2
                if d2 < best or (d2 == best and iso_idx_s[k] < b):
                    best = d2
                    b = iso_idx_s[k]
                k += 1
            # Other half: |m - mag| grows as k decreases
            k = j - 1
            while k >= 0:
                dm2 = (m - iso_mag_s[k]) ** 2
                if dm2 > best:
                    break
                d2 = dm2 + (c - iso_col_s[k]) ** 2
                if d2 < best or (d2 == best and iso_idx_s[k] < b):
                    best = d2
                    b = iso_idx_s[k]
                k -= 1
            # All-NaN isochrone: fall back to index 0 (as np.argmin would)
            out[i] = b if b < Ni else 0
        return out

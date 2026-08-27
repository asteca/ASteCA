import warnings

import numpy as np
from scipy.signal import savgol_filter
from scipy.spatial import KDTree


def to_point_find(
    x: np.ndarray,
    y: np.ndarray,
    isochs_model: str,
    fallback_min_idx: int,
    window: int = 21,
    polyorder: int = 3,
    deriv: int = 1,
    eps_default: float = 1e-6,
    eps_mist: float = 1e-3,
    confirm_min: int = 3,
    confirm_frac: float = 1 / 3,
) -> tuple[float, float, int]:
    """Estimate turn-off point for an isochrone.

    The function applies a Savitzky-Golay filter to estimate the first derivative of
    the isochrone and identifies the first sustained increasing region.

    :param x: One-dimensional array of isochrone color values.
    :type x: np.ndarray
    :param y: One-dimensional array of isochrone magnitudes values.
    :type y: np.ndarray
    :param isochs_model: Name of the isochrone model.
    :type isochs_model: str
    :param fallback_min_idx: Minimum acceptable index for the detected turn-off.
    :type fallback_min_idx: int
    :param window: Window length for the Savitzky-Golay filter.
        - Smaller window:
            - Less smoothing
            - More sensitive to local features and noise
            - Can detect short-lived changes in the isochrone
            - More likely to produce an early/unstable turn-off
        - Larger window:
            - More smoothing
            - Ignores small-scale structure
            - Requires a broader, more persistent change in color
            - Produces a more robust but potentially later/less precise detection
    :type window: int
    :param polyorder: Polynomial order for the Savitzky-Golay filter.
        - Lower polyorder (e.g. 2):
            - More restrictive local fit
            - Smoother derivative
            - Less sensitive to complicated curvature
        - Higher polyorder(e.g. 4 or 5):
            - Can reproduce more complicated local shapes
            - Potentially more sensitive to real curvature
            - Also potentially more sensitive to noise
    :type polyorder: int
    :param deriv: Derivative order for the Savitzky-Golay filter.
    :type deriv: int
    :param eps_default: Default derivative threshold for an increasing trend.
    :type eps_default: float
    :param eps_mist: Derivative threshold used for the ``mist`` model.
    :type eps_mist: float
    :param confirm_min: Minimum number of consecutive derivative values required.
    :type confirm_min: int
    :param confirm_frac: Fraction of ``window`` used to define consecutive values
        required to confirm the trend.
    :type confirm_frac: float

    :return: Color and magnitude for the turn-off point.
    :rtype: tuple[float, float, int]
    """
    n = len(x)
    eps = eps_default

    # The MIST isochrones require a much larger epsilon value to identify the turn-off
    if isochs_model.lower() == "mist":
        eps = eps_mist

    if n < window:
        to_idx = n - 1
    else:
        if window % 2 == 0:
            window += 1
        dx = savgol_filter(x, window_length=window, polyorder=polyorder, deriv=deriv)
        confirm = max(confirm_min, int(window * confirm_frac))
        to_idx = n - 1
        for i in range(n - confirm):
            if np.all(dx[i : i + confirm] > eps):
                to_idx = i
                break

    # If no sustained positive trend is found, use the leftmost point in the isochrone.
    # This is a fallback mechanism to ensure that a semi-reasonable turn-off point is
    # assigned if the found value is too far down the isochrone
    if to_idx < fallback_min_idx:
        to_idx = int(np.argmin(x))

    # Extract the color and magnitude values at the detected turn-off index
    to_col, to_mag = float(x[to_idx]), float(y[to_idx])

    return to_col, to_mag, to_idx


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


def get_bss_probabilities(
    cluster_mag: np.ndarray,
    cluster_colors: list[np.ndarray | None],
    turn_off_points: dict,
    N_sampled_synthcls: int,
    mag_offset_1: float,
    mag_offset_2: float,
    col_offset_1: float,
    col_offset_2: float,
    percentile: float = 5,
    dmag_init: float = 0.25,
    dmag_step: float = 0.05,
    min_stars: int = 5,
) -> np.ndarray:
    """Estimate the probability of each observed star being a blue straggler star
    (BSS) based on the turn-off point of the isochrones generated from the sampled
    models.

    :param cluster_mag: Observed cluster magnitudes.
    :type cluster_mag: np.ndarray
    :param cluster_colors: Observed cluster colors.
    :type cluster_colors: list[np.ndarray | None]
    :param turn_off_points: Turn-off region limits estimated for sampled models.
    :type turn_off_points: dict
    :param N_sampled_synthcls: Number of sampled synthetic clusters.
    :type N_sampled_synthcls: int
    :param mag_offset_1: Offset applied to the turn-off magnitude to define the
        first BSS region
    :type mag_offset_1: float
    :param mag_offset_2: Offset applied to the turn-off magnitude to define the
        second BSS region
    :type mag_offset_2: float
    :param col_offset_1: Offset applied to the turn-off color to define the
        first BSS region
    :type col_offset_1: float
    :param col_offset_2: Offset applied to the turn-off color to define the
        second BSS region
    :type col_offset_2: float
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
    cl_color = cluster_colors[turn_off_points["color_idx"]]

    # Define magnitude limits around the detected turn-off magnitude
    to_mag_1 = to_mag + mag_offset_1
    to_mag_2 = to_mag + mag_offset_2

    bss_probs = np.zeros_like(cluster_mag)
    for to_c, to_m, to_m1, to_m2 in zip(to_col, to_mag, to_mag_1, to_mag_2):
        # Define a magnitude range around the identified turn-off magnitude.
        # Iteratively increase the magnitude range until at least `min_stars` are found
        # within it.
        dmag = dmag_init
        while True:
            mag_range = (cluster_mag > (to_m - dmag)) & (cluster_mag < (to_m + dmag))
            if mag_range.sum() >= min_stars:
                break
            dmag += dmag_step

        # This approach is more robust to incorrect parameter values assigned to
        # the observed cluster
        to_col_cl = np.nanpercentile(cl_color[mag_range], percentile)
        to_c1 = min(to_col_cl, to_c + col_offset_1)
        to_c2 = min(to_col_cl, to_c + col_offset_2)

        # Assign BSS probabilities
        msk = ((cl_color < to_c1) & (cluster_mag < to_m1)) | (
            (cl_color < to_c2) & (cluster_mag < to_m2)
        )
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
    :param extrinsic_pars: Dictionary containing model parameters.
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
    rand_norm = np.array([])
    rand_unif = np.array([])
    DR_distribution = ""
    # Set DR=0 to avoid disrupting the isochrone.
    dr = 0.0

    # Sample extrinsic parameter values
    Na, Nm = [int(_) for _ in np.shape(isoch_youngs_arrs)[:2]]
    N_samples = max(1, int(N_models / (Na * Nm) + 1))
    sampled_models = {
        k: rng.normal(
            f_val,
            extrinsic_pars["stddevs"][k],
            N_samples,
        )
        for k, f_val in extrinsic_pars["values"].items()
    }

    n_stars = cluster_mag.size

    # One column per isochrone/model realization.
    masses = []

    valid = np.isfinite(cluster_mag) & np.isfinite(cluster_color)

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
                isoch = sm.extinction(
                    ext_law,
                    ext_coefs,
                    rand_norm,
                    rand_unif,
                    DR_distribution,
                    m_ini_idx,
                    binar_flag,
                    av,
                    dr,
                    rv,
                    isoch_moved,
                )

                isoch_mag = np.asarray(isoch[0])
                isoch_color = np.asarray(isoch[1])
                isoch_mass = np.asarray(isoch[2])

                # Distance from every cluster star to every isochrone point.
                d2 = (cluster_mag[valid, None] - isoch_mag[None, :]) ** 2 + (
                    cluster_color[valid, None] - isoch_color[None, :]
                ) ** 2

                # Closest point on THIS isochrone for each star.
                closest_idx = np.argmin(d2, axis=1)

                mass = np.full(n_stars, np.nan)
                mass[valid] = isoch_mass[closest_idx]
                masses.append(mass)

    # (N_stars, N_realizations)
    bss_masses = np.column_stack(masses)

    # Extract medians and STDDEVs
    n_valid = np.sum(np.isfinite(bss_masses), axis=1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        bss_mass_median = np.nanmedian(bss_masses, axis=1)
        bss_mass_stddev = np.nanstd(bss_masses, axis=1)
    # 0 valid masses: median = NaN, stddev = NaN
    # 1 valid mass: median = that mass, stddev = NaN
    # ≥2 valid masses: median and stddev are calculated normally.
    bss_mass_median[n_valid == 0] = np.nan
    bss_mass_stddev[n_valid < 2] = np.nan

    return bss_mass_median, bss_mass_stddev

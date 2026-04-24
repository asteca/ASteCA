import warnings

import numpy as np
from scipy.signal import savgol_filter
from scipy.spatial import KDTree


def to_point_find(
    x: np.ndarray,
    y: np.ndarray,
    isochs_model: str,
    window: int = 21,
    polyorder: int = 3,
    deriv: int = 1,
    eps_default: float = 1e-6,
    eps_mist: float = 1e-3,
    confirm_min: int = 3,
    confirm_frac: float = 1 / 3,
    fallback_min_idx: int = 10,
) -> tuple[float, float]:
    """Estimate turn-off point for an isochrone.

    The function applies a Savitzky-Golay filter to estimate the first derivative of
    the isochrone and identifies the first sustained increasing region.

    :param x: One-dimensional array of isochrone color values.
    :type x: np.ndarray
    :param y: One-dimensional array of isochrone magnitudes corresponding to ``x``.
    :type y: np.ndarray
    :param isochs_model: Name of the isochrone model.
    :type isochs_model: str
    :param window: Window length for the Savitzky-Golay filter.
    :type window: int
    :param polyorder: Polynomial order for the Savitzky-Golay filter.
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
    :param fallback_min_idx: Minimum acceptable index for the detected turn-off.
    :type fallback_min_idx: int

    :return: Color and magnitude for the turn-off point.
    :rtype: tuple[float, float]
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

    # If no sustained positive trend is found, assign to the index that corresponds to
    # the middle point in magnitude. This is a fallback mechanism to ensure that a
    # reasonable turn-off point is assigned when the found value is too close to the
    # bottom of the isochrone
    if to_idx < fallback_min_idx:
        # Assign to the index that corresponds to the middle point in magnitude
        y_mid_point = 0.5 * (max(y) + min(y))
        to_idx = np.argmin(np.abs(y - y_mid_point))

    # Extract the color and magnitude values at the detected turn-off index
    to_col, to_mag = float(x[to_idx]), float(y[to_idx])

    return to_col, to_mag


def get_stellar_masses(
    cluster_mag: np.ndarray,
    cluster_colors: list,
    m_ini_idx: int,
    sampled_synthcls: list,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Assign primary and secondary (synthetic) masses to each observed star,
    for each synthetic model generated.

    :param cluster_mag: Observed magnitude.
    :type cluster_mag: np.ndarray
    :param cluster_colors: Observed colors.
    :type cluster_colors: list
    :param m_ini_idx: Index of the initial mass column
    :type m_ini_idx: int
    :param sampled_synthcls: Sampled synthetic cluster data.
    :type sampled_synthcls: list

    :return: Primary and secondary masses values (median + stddev), binary
        probability per observed star, and a mask for NaN values.
    :rtype: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]
    """

    # Extract observed photometry
    cl_colors = [cluster_colors[0]]
    if cluster_colors[1] is not None:
        cl_colors.append(cluster_colors[1])
    obs_phot = np.array([cluster_mag] + [_ for _ in cl_colors])
    # Replace nans in mag and colors to avoid crashing KDTree()
    nan_msk = np.full(obs_phot.shape[1], False)
    for ophot in obs_phot:
        nan_msk = nan_msk | np.isnan(ophot)
    obs_phot[:, nan_msk] = -10.0
    obs_phot = obs_phot.T

    # Assign primary and secondary (synthetic) masses to each observed star,
    # for each synthetic model generated
    m12_obs = []
    # weights = []
    for isoch in sampled_synthcls:
        # Indexes of the photometrically closest synthetic stars to the observed stars
        tree = KDTree(isoch[:m_ini_idx].T)
        dist, close_stars_idxs = tree.query(obs_phot, k=1)

        # # Store inverse normalized distance used as weights
        # inv_dist = 1 / dist
        # inv_dist_norm = inv_dist / inv_dist.max()
        # weights.append(inv_dist_norm)

        # The secondary mass is stored after the primary mass, hence the ':'
        m12_obs.append(isoch[m_ini_idx:, close_stars_idxs])

    # Extract primary and secondary masses
    m12_obs = np.array(m12_obs)
    m1_obs = m12_obs[:, 0, :]
    m2_obs = m12_obs[:, 1, :]

    # Primary mass values (median + stddev)
    m1_med = np.median(m1_obs, 0)
    # # Weighted average is similar to median
    # m1_med = np.average(m1_obs, 0, weights)
    m1_std = np.std(m1_obs, 0)

    # Secondary mass values (median + stddev)
    # Hide 'All-nan slice' warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m2_med = np.nanmedian(m2_obs, 0)
        m2_std = np.nanstd(m2_obs, 0)
    # m2 can not be larger than m1
    m2_med = np.min([m1_med, m2_med], 0)

    # m2_obs is used to compute the binary probability per observed star
    return m1_med, m1_std, m2_med, m2_std, nan_msk, m2_obs


def get_binar_probabilities(m2_obs: np.ndarray, N_sampled_synthcls: int) -> np.ndarray:
    """Binary probability per observed star

    Count how many times the secondary mass of an observed star was assigned a
    value 'not np.nan', i.e.: was identified as a binary system. Dividing this
    value by the number of synthetic models used, results in the per observed
    star probability of being a binary system.

    :param m2_obs: Secondary masses assigned to observed stars across sampled models.
    :type m2_obs: np.ndarray
    :param N_sampled_synthcls: Number of sampled synthetic clusters.
    :type N_sampled_synthcls: int
    :return: Array with the binary probability for each observed star
    :rtype: np.ndarray
    """
    binar_prob = (~np.isnan(m2_obs)).sum(0) / N_sampled_synthcls

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

    return bss_probs

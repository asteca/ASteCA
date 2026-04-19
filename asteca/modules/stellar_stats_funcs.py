import warnings

import numpy as np
from scipy.signal import savgol_filter
from scipy.spatial import KDTree


def to_point_find(
    x: np.ndarray,
    y: np.ndarray,
    cluster_mag: np.ndarray,
    cluster_color: np.ndarray,
    isochs_model: str,
    window: int = 21,
    polyorder: int = 3,
    deriv: int = 1,
    dmag_init: float = 0.25,
    dmag_step: float = 0.05,
    min_stars: int = 5,
    eps_default: float = 1e-6,
    eps_mist: float = 1e-3,
    confirm_min: int = 3,
    confirm_frac: float = 1 / 3,
    fallback_min_idx: int = 10,
    mag_offset: float = 0.5,
    col_offset_1: float = 0.1,
    col_offset_2: float = 0.05,
    percentile: float = 5,
) -> tuple[float, float, float, float]:
    """
    Identifies the first index where a signal begins a sustained monotonic increase.
    In this case, it identifies the color value where it begins to increase,
    i.e. the turn-off point.

    The function applies a Savitzky-Golay filter to estimate the first derivative
    of the input signal. It then searches for the first sequence of length `confirm`
    where all estimated derivatives are strictly positive.

    Args:
        x (np.ndarray): 1D array of isochrone color values (independent variable).
        y (np.ndarray): 1D array of isochrone magnitudes corresponding to `x`.
        cluster_mag (np.ndarray): Observed cluster magnitudes used to refine the
            turn-off region selection.
        cluster_color (np.ndarray): Observed cluster colors corresponding to
            `cluster_mag`.
        isochs_model (str): Name of the isochrone model. Used to adjust the
            derivative threshold for specific models (e.g. "mist").

        window (int, optional): Window length for the Savitzky-Golay filter.
            Must be a positive odd integer. If even, it is incremented by 1.
        polyorder (int, optional): Polynomial order used by the Savitzky-Golay filter.
            Must be smaller than `window`.
        deriv (int, optional): Order of the derivative computed by the
            Savitzky-Golay filter.
        dmag_init (float, optional): Initial half-width of the magnitude interval
            around the detected turn-off magnitude used to select cluster stars.
        dmag_step (float, optional): Increment applied to the magnitude interval
            half-width until at least `min_stars` are enclosed.
        min_stars (int, optional): Minimum number of observed cluster stars required
            within the magnitude interval around the turn-off point.
        eps_default (float, optional): Default minimum derivative value required to
            consider the signal locally increasing.
        eps_mist (float, optional): Derivative threshold used when the isochrone
            model is "mist".
        confirm_min (int, optional): Minimum number of consecutive derivative values
            required to confirm a sustained positive trend.
        confirm_frac (float, optional): Fraction of `window` used to determine the
            number of consecutive derivative values required to confirm the trend.
        fallback_min_idx (int, optional): Minimum acceptable index for the detected
            turn-off. If the detected index is smaller than this value, a fallback
            based on the magnitude midpoint is used.
        mag_offset (float, optional): Offset applied to the detected magnitude to
            define the upper and lower magnitude limits returned.
        col_offset_1 (float, optional): Offset applied to the detected color to
            define the first returned color limit.
        col_offset_2 (float, optional): Offset applied to the detected color to
            define the second returned color limit.
        percentile (float, optional): Percentile of the observed cluster colors
            within the magnitude interval used to estimate the turn-off color.

    Returns:
        tuple[float, float, float, float]:
            to_col_1 : Lower color bound for the turn-off region.
            to_col_2 : Upper color bound for the turn-off region.
            to_mag_1 : Upper magnitude bound for the turn-off region.
            to_mag_2 : Lower magnitude bound for the turn-off region.
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
    to_col, to_mag = x[to_idx], y[to_idx]

    # Define magnitude limits around the detected turn-off magnitude
    to_mag_1 = to_mag + mag_offset
    to_mag_2 = to_mag - mag_offset

    # Define a magnitude range around the identified turn-off magnitude.
    # Iteratively increase the magnitude range until at least `min_stars` are found
    # within it.
    dmag = dmag_init
    while True:
        mag_range = (cluster_mag > (to_mag - dmag)) & (cluster_mag < (to_mag + dmag))
        if mag_range.sum() >= min_stars:
            break
        dmag += dmag_step

    # This approach is more robust to incorrect parameter values assigned to
    # the observed cluster
    to_col_cl = np.nanpercentile(cluster_color[mag_range], percentile)
    to_col_1 = min(to_col_cl, to_col - col_offset_1)
    to_col_2 = min(to_col_cl, to_col - col_offset_2)

    return to_col_1, to_col_2, to_mag_1, to_mag_2


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


def get_binar_probabilities(m2_obs, N_sampled_synthcls) -> np.ndarray:
    """Binary probability per observed star

    Count how many times the secondary mass of an observed star was assigned a
    value 'not np.nan', i.e.: was identified as a binary system. Dividing this
    value by the number of synthetic models used, results in the per observed
    star probability of being a binary system.

    :return: Array with the binary probability for each observed star
    :rtype: np.ndarray
    """
    binar_prob = (~np.isnan(m2_obs)).sum(0) / N_sampled_synthcls

    return binar_prob


def get_bss_probabilities(
    cluster_mag, cluster_colors, turn_off_points, N_sampled_synthcls
) -> np.ndarray:
    """Estimate the probability of each observed star being a blue straggler star
    (BSS) based on the turn-off point of the isochrones generated from the sampled
    models.

    :return: Array with the BSS probability for each observed star
    :rtype: np.ndarray
    """
    color_idx = turn_off_points["color_idx"]
    to_col_1 = np.array(turn_off_points["to_col_1"])
    to_col_2 = np.array(turn_off_points["to_col_2"])
    to_mag_1 = np.array(turn_off_points["to_mag_1"])
    to_mag_2 = np.array(turn_off_points["to_mag_2"])

    color = cluster_colors[color_idx]
    mag = cluster_mag

    # Assign BSS probabilities
    bss_probs = np.zeros_like(mag)
    for c1, c2, m1, m2 in zip(to_col_1, to_col_2, to_mag_1, to_mag_2):
        msk = ((color < c1) & (mag < m1)) | ((color < c2) & (mag < m2))
        bss_probs[msk] += 1

    # Normalize probability
    bss_probs /= N_sampled_synthcls

    return bss_probs

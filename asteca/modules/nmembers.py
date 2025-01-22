import warnings

import numpy as np
from astropy.stats import RipleysKEstimator

from . import cluster_priv as cp


def density_nmembs(
    x: np.ndarray, y: np.ndarray, center: tuple[float, float], radius: float
) -> int:
    """Estimate the number of cluster members based on a density calculation.

    :param x: Array of x-coordinates.
    :type x: np.ndarray
    :param y: Array of y-coordinates.
    :type y: np.ndarray
    :param center: Center coordinates (x, y).
    :type center: tuple[float, float]
    :param radius: Radius of the cluster.
    :type radius: float

    :return: Estimated number of cluster members.
    :rtype: int
    """
    dist = np.sqrt((x - center[0]) ** 2 + (y - center[1]) ** 2)
    msk_in_rad = dist <= radius
    N_cl_region = msk_in_rad.sum()

    N_total = len(x)
    N_field = N_total - N_cl_region

    dim_x = np.percentile(x, (1, 99))
    dim_y = np.percentile(y, (1, 99))
    A_total = (dim_x[1] - dim_x[0]) * (dim_y[1] - dim_y[0])

    A_cl_region = np.pi * radius**2
    A_field = A_total - A_cl_region
    dens_field = N_field / A_field

    N_field_in_cl_region = int(dens_field * A_cl_region)
    n_memb = int(N_cl_region - N_field_in_cl_region)

    # # n_all = []
    # for rad in np.linspace(0.1*radius, radius, 100):
    #     dist = np.sqrt((x - center[0])**2 + (y - center[1])**2)
    #     msk_in_rad = dist <= rad
    #     N_cl_region = msk_in_rad.sum()

    #     A_cl_region = np.pi*rad**2
    #     A_field = A_total - A_cl_region
    #     dens_field = N_field / A_field

    #     n_field = int(dens_field * A_cl_region)
    #     n_memb = int(N_cl_region - n_field)
    #     if n_memb < 0:
    #         break

    #     # n_all.append([n_memb, n_field])
    #     print(round(rad, 3), N_cl_region, n_field, n_memb)

    return n_memb


def ripley_nmembs(
    x: np.ndarray,
    y: np.ndarray,
    pmRA: np.ndarray,
    pmDE: np.ndarray,
    plx: np.ndarray,
    vpd_c: tuple[float, float],
    plx_c: float,
    N_clust: int = 50,
    N_extra: int = 5,
    N_step: int = 10,
) -> int:
    """Estimate the number of cluster members using Ripley's K-function.

    :param x: Array of x-coordinates.
    :type x: np.ndarray
    :param y: Array of y-coordinates.
    :type y: np.ndarray
    :param pmRA: Array of proper motion in right ascension.
    :type pmRA: np.ndarray
    :param pmDE: Array of proper motion in declination.
    :type pmDE: np.ndarray
    :param plx: Array of parallax values.
    :type plx: np.ndarray
    :param vpd_c: Center coordinates in proper motion (pmRA, pmDE).
    :type vpd_c: tuple[float, float]
    :param plx_c: Center parallax value.
    :type plx_c: float
    :param N_clust: Initial number of stars to consider as a cluster, defaults to 50
    :type N_clust: int
    :param N_extra: Number of extra iterations to perform if the initial
        clustering fails, defaults to 5
    :type N_extra: int
    :param N_step: Step size for increasing the number of cluster stars, defaults to 10
    :type N_step: int

    :return: Estimated number of cluster members.
    :rtype: int
    """
    rads, Kest, C_thresh_N = init_ripley(x, y)

    # Obtain the ordered indexes of the distances to the (pmra, pmde, plx) center
    cents_3d = np.array([list(vpd_c) + [plx_c]])
    data_3d = np.array([pmRA, pmDE, plx]).T
    d_pm_plx_idxs = cp.get_Nd_dists(cents_3d, data_3d)

    # Select those clusters where the stars are different enough from a
    # random distribution
    xy = np.array([x, y]).T
    idx_survived = ripley_core(rads, Kest, C_thresh_N, d_pm_plx_idxs, xy, N_clust)

    # If the default clustering number did not work, try a few
    # more values with an increasing number of cluster stars
    if not idx_survived:
        for _ in range(N_extra):
            N_clust_surv = int(N_clust + (_ + 1) * N_step)
            idx_survived = ripley_core(
                rads, Kest, C_thresh_N, d_pm_plx_idxs, xy, N_clust_surv
            )
            # Break out when (if) any value selected stars
            if len(idx_survived) > 0:
                break

    N_survived = len(idx_survived)

    return N_survived


def init_ripley(
    lon: np.ndarray, lat: np.ndarray
) -> tuple[np.ndarray, RipleysKEstimator, float]:
    """Initialize Ripley's K-function estimator.

    https://rdrr.io/cran/spatstat/man/Kest.html
    "For a rectangular window it is prudent to restrict the r values to a
    maximum of 1/4 of the smaller side length of the rectangle
    (Ripley, 1977, 1988; Diggle, 1983)"

    :param lon: Array of longitude values.
    :type lon: np.ndarray
    :param lat: Array of latitude values.
    :type lat: np.ndarray

    :return: Radii, Ripley's K-estimator, and threshold value.
    :rtype: tuple[np.ndarray, RipleysKEstimator, float]
    """
    xmin, xmax = lon.min(), lon.max()
    ymin, ymax = lat.min(), lat.max()
    area = (xmax - xmin) * (ymax - ymin)
    Kest = RipleysKEstimator(area=area, x_max=xmax, y_max=ymax, x_min=xmin, y_min=ymin)

    # Ripley's rule of thumb
    thumb = 0.25 * min((xmax - xmin), (ymax - ymin))
    # Large sample rule
    rho = len(lon) / area
    large = np.sqrt(1000 / (np.pi * rho))

    lmin = min(thumb, large)
    rads = np.linspace(lmin * 0.01, lmin, 100)  # HARDCODED

    C_thresh_N = 1.68 * np.sqrt(area)  # HARDCODED

    return rads, Kest, C_thresh_N


def ripley_core(
    rads: np.ndarray,
    Kest: RipleysKEstimator,
    C_thresh_N: float,
    d_pm_plx_idxs: np.ndarray,
    xy: np.ndarray,
    N_clust: int,
    N_break: int = 5,
) -> list[int]:
    """Core function that estimates the number of members based on Ripley's K-function.

    :param rads: Array of radii.
    :type rads: np.ndarray
    :param Kest: Ripley's K-estimator.
    :type Kest: RipleysKEstimator
    :param C_thresh_N: Threshold value.
    :type C_thresh_N: float
    :param d_pm_plx_idxs: Ordered indexes of the distances to the (pmra, pmde, plx)
     center.
    :type d_pm_plx_idxs: np.ndarray
    :param xy: Array of (x, y) coordinates.
    :type xy: np.ndarray
    :param N_clust: Number of stars to consider as a cluster.
    :type N_clust: int
    :param N_break: Number of breaks before stopping the loop, defaults to 5
    :type N_break: int

    :return: List of indexes of the survived stars.
    :rtype: list[int]
    """
    N_total_stars = xy.shape[0]

    # Define K-function threshold
    C_thresh = C_thresh_N / N_clust

    idx_survived = []
    N_break_count, step_old = 0, 0
    for step in np.arange(N_clust, N_total_stars, N_clust):
        # Ring of stars around the VPD+Plx centers
        msk_ring = d_pm_plx_idxs[step_old:step]
        # Obtain their Ripely K estimator
        C_s = rkfunc(xy[msk_ring], rads, Kest)

        if not np.isnan(C_s):
            # This group of stars survived
            if C_s >= C_thresh:
                idx_survived += list(msk_ring)
            else:
                # Increase break condition
                N_break_count += 1
        if N_break_count > N_break:
            break
        step_old = step
    return idx_survived


def rkfunc(xy: np.ndarray, rads: np.ndarray, Kest: RipleysKEstimator) -> float:
    """Test how similar this cluster's (x, y) distribution is compared
    to a uniform random distribution using Ripley's K.

    https://stats.stackexchange.com/a/122816/10416

    :param xy: Array of (x, y) coordinates.
    :type xy: np.ndarray
    :param rads: Array of radii.
    :type rads: np.ndarray
    :param Kest: Ripley's K-estimator.
    :type Kest: RipleysKEstimator
    :return: Ripley's K-function value.
    :rtype: float
    """
    # Avoid large memory consumption if the data array is too big
    # if xy.shape[0] > 5000:
    #     mode = "none"
    # else:
    #     mode = 'translation'

    # Hide RunTimeWarning
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        L_t = Kest.Lfunction(xy, rads, mode="translation")

    # Catch all-nans. Avoid 'RuntimeWarning: All-NaN slice encountered'
    if np.isnan(L_t).all():
        C_s = np.nan
    else:
        C_s = np.nanmax(abs(L_t - rads))

    return C_s

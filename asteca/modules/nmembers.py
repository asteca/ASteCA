import warnings
from typing import Literal

import numpy as np
from astropy.stats import RipleysKEstimator
from scipy.spatial import KDTree

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
    xy_center: tuple[float, float],
    vpd_c: tuple[float, float],
    plx_c: float,
    N_clust_min: int,
    N_clust_max: int,
    N_clust: int = 50,
    N_extra: int = 5,
    N_step: int = 10,
    N_close: int = 10,
    X_stars_max: float = 5,
    N_stars_max: int = 10_000,
) -> np.ndarray:
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
    :param xy_center: Center coordinates in (x,y)
    :type xy_center: tuple[float, float]
    :param vpd_c: Center coordinates in proper motion (pmRA, pmDE).
    :type vpd_c: tuple[float, float]
    :param plx_c: Center parallax value.
    :type plx_c: float
    :param N_clust_min: Minimum number of stars to consider as a cluster
    :type N_clust_min: int
    :param N_clust_max: Maximum number of stars to consider as a cluster
    :type N_clust_max: int
    :param N_clust: Initial number of stars to consider as a cluster, defaults to 50
    :type N_clust: int
    :param N_extra: Number of extra iterations to perform if the initial
        clustering fails, defaults to 5
    :type N_extra: int
    :param N_step: Step size for increasing the number of cluster stars, defaults to 10
    :type N_step: int
    :param N_close: Number of stars to select if no star survives the Ripley cleaning
    :type N_close: int
    :param X_stars_max: Multiplicative term to determine the maximum number of stars
        that should pass the first filter
    :type X_stars_max: float
    :param N_stars_max: Absolute maximum number of stars that should pass the first filter
    :type N_stars_max: int

    :return: Indexes of the estimated cluster members.
    :rtype: np.ndarray
    """
    # Generate input data array
    X = np.array(
        [
            x,
            y,
            pmRA,
            pmDE,
            plx,
        ]
    )
    idx_clean, X_no_nan = cp.reject_nans(X)
    # Unpack input data with no 'nans'
    lon, lat, pmRA, pmDE, plx = X_no_nan
    # These arrays are not used, fill with dummy values
    e_pmRA, e_pmDE, e_plx = [np.empty(len(lon))] * 3

    # Remove the most obvious field stars
    N_filter_max = int(min(N_clust_max * X_stars_max, N_stars_max))
    idx_clean, x, y, pmRA, pmDE, plx, _, _, _ = cp.first_filter(
        idx_clean,
        vpd_c,
        plx_c,
        lon,
        lat,
        pmRA,
        pmDE,
        plx,
        e_pmRA,
        e_pmDE,
        e_plx,
        N_filter_max,
    )

    # Initialize Ripley's K-function estimator
    rads, Kest, C_thresh_N = init_ripley(x, y)

    # Obtain the ordered indexes of the distances to the (pmra, pmde, plx) center
    cents_3d = np.array([list(vpd_c) + [plx_c]])
    d_pm_plx_idxs = cp.get_Nd_dists(cents_3d, np.array([pmRA, pmDE, plx]).T)

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

    # Check number of surviving stars
    if len(idx_survived) == 0:
        # If no stars survived, select the N_close closest to the xy center
        dist_to_center = np.linalg.norm(xy - xy_center, axis=1)
        idx_survived = np.argsort(dist_to_center)[:N_close]
        return idx_clean[idx_survived]

    elif len(idx_survived) < N_clust_min:
        # No need for further cleaning
        return idx_clean[idx_survived]

    # Apply local density based cleaning
    dist_to_cent, local_density, fdens = estimate_field_density(
        xy_center, x, y, idx_survived
    )
    idx_survived = apply_density_clean(idx_survived, dist_to_cent, local_density, fdens)

    return idx_clean[idx_survived]


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

    # See 'Tests of randomness for spatial point patterns', Ripley (1979)
    C_thresh_N = 1.68 * np.sqrt(area)

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


def rkfunc(
    xy: np.ndarray,
    rads: np.ndarray,
    Kest: RipleysKEstimator,
    mode: Literal["none", "translation", "ohser", "ripley"] = "translation",
) -> float | np.floating:
    """Test how similar this cluster's (x, y) distribution is compared
    to a uniform random distribution using Ripley's K.

    To avoid large memory consumption if the data array is too big, use:
    mode = "none"

    https://stats.stackexchange.com/a/122816/10416

    :param xy: Array of (x, y) coordinates.
    :type xy: np.ndarray
    :param rads: Array of radii.
    :type rads: np.ndarray
    :param Kest: Ripley's K-estimator.
    :type Kest: RipleysKEstimator
    :param mode: Mode for Ripley's K-estimator.
    :type mode: str

    :return: Ripley's K-function value.
    :rtype: float | np.floating
    """

    # Hide RunTimeWarning
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        L_t = Kest.Lfunction(xy, rads, mode=mode)

    # Catch all-nans. Avoid 'RuntimeWarning: All-NaN slice encountered'
    if np.isnan(L_t).all():
        C_s = np.nan
    else:
        C_s = np.nanmax(abs(L_t - rads))

    return C_s


def estimate_field_density(
    xy_center: tuple[float, float],
    x: np.ndarray,
    y: np.ndarray,
    idx: list[int],
    N_k: int = 5,
    p_bounds: tuple[float, float, float] = (75, 85, 95),
    large_cluster_N: int = 250,
    fdens_ratio_bounds: tuple[float, float] = (0.2, 0.7),
) -> tuple[np.ndarray, np.ndarray, float]:
    """
    Estimate local stellar density and field density level.

    Local stellar density is estimated from the distance to the ``N_k``-th nearest
    neighbor. A field density level is inferred from outer regions of the cluster,
    and stars with lower densities or large distances from the center are rejected.

    The 'large_cluster_N' and 'fdens_ratio_bounds' values are selected so that for
    clusters with at least 250 identified members, if the 95th percentile field density
    is between 20% and 70% of the 75th field density, then the 95th percentile is used.
    This helps to select more members for large clusters, while also avoiding
    over-cleaning for the remaining clusters. The 20% lower bound is there to catch
    clusters with a very low number of stars at the edges.

    :param xy_center: Cartesian coordinates of the cluster center ``(x_c, y_c)``.
    :type xy_center: tuple(float, float)
    :param x: X coordinates of all stars.
    :type x: ndarray
    :param y: Y coordinates of all stars.
    :type y: ndarray
    :param idx: Indices of stars used in the density estimation.
    :type idx: list(int)
    :param N_k: Neighbor order used to estimate local density.
    :type N_k: int
    :param p_bounds: Percentiles of the distance distribution used to estimate
        the field density region ``(p_low, p_mid, p_high)``.
    :type p_bounds: tuple(float, float, float)
    :param large_cluster_N: Minimum number of members required to apply the
        alternative field density estimate.
    :type large_cluster_N: int
    :param fdens_ratio_bounds: Accepted ratio range between outer and inner
        field densities that triggers the alternative estimate.
    :type fdens_ratio_bounds: tuple(float, float)

    :return:
        dist_to_cent : distances to cluster center
        local_density : normalized local densities
        fdens : estimated field density level
    :rtype: tuple(ndarray, ndarray, float)
    """

    data = np.array([x[idx], y[idx]]).T

    dist_to_cent = np.linalg.norm(data - xy_center, axis=1)

    # Local density from N_k nearest neighbors
    tree = KDTree(data)
    dists, _ = tree.query(data, k=N_k + 1)
    kth_nn_dist = np.asarray(dists)[:, N_k]
    local_density = 1 / kth_nn_dist
    local_density /= local_density.max()

    # Field density estimate from outer region
    p_low, p_mid, p_high = np.percentile(dist_to_cent, p_bounds)
    msk_mid = (dist_to_cent > p_low) & (dist_to_cent < p_mid)
    fdens = np.median(local_density[msk_mid])

    # Alternative estimate for large clusters
    fdens_outer = np.median(local_density[dist_to_cent > p_high])
    if fdens > 0:
        ratio = fdens_outer / fdens
    else:
        ratio = np.inf
    if len(idx) > large_cluster_N:
        if fdens_ratio_bounds[0] < ratio < fdens_ratio_bounds[1]:
            fdens = fdens_outer

    return dist_to_cent, local_density, fdens


def apply_density_clean(
    idx: list[int],
    dist_to_cent: np.ndarray,
    local_density: np.ndarray,
    fdens: float,
    band_init: tuple[float, float] = (0.95, 1.05),
    band_step: float = 0.05,
    band_min_count: int = 5,
    rad_percentile: float = 5,
) -> list[int]:
    """
    Reject stars with low density or large distance from the center.

    The method attempts to identify the radius at which the cluster density
    becomes indistinguishable from the field density. Stars beyond this radius
    or with densities below the estimated field level are rejected.

    :param idx: Indices of stars to clean.
    :type idx: list(int)
    :param dist_to_cent: Distance of each star to the cluster center.
    :type dist_to_cent: ndarray
    :param local_density: Normalized local stellar density.
    :type local_density: ndarray
    :param fdens: Estimated field density level.
    :type fdens: float
    :param band_init: Initial multiplicative bounds around the field density.
    :type band_init: tuple(float, float)
    :param band_step: Increment applied to widen the density band if needed.
    :type band_step: float
    :param band_min_count: Minimum number of stars required in the density band.
    :type band_min_count: int
    :param rad_percentile: Percentile used to estimate limiting radius.
    :type rad_percentile: float

    :return: Updated indices of surviving stars.
    :rtype: list(int)
    """
    # Determine radius where density ~ field density
    f_low, f_high = band_init
    while True:
        msk_fdens = (local_density > fdens * f_low) & (local_density < fdens * f_high)

        if msk_fdens.sum() >= band_min_count:
            break
        # If no stars where found, make the band larger
        f_low -= band_step
        f_high += band_step

    # The radius is the rad_percentile of the stars in the band around fdens
    rad = np.percentile(dist_to_cent[msk_fdens], rad_percentile)

    # The final mask rejects stars beyond the radius and with a local density
    # below the field density
    msk = (local_density > fdens) & (dist_to_cent < rad)
    if msk.sum() == 0:
        return idx

    return list(np.array(idx)[msk])

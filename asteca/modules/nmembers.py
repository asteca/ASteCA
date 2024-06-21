import warnings
import numpy as np
from astropy.stats import RipleysKEstimator
from . import cluster_priv as cp


def density_nmembs(x, y, center, radius):
    """ """
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

    n_field = int(dens_field * A_cl_region)
    n_memb = int(N_cl_region - n_field)

    # n_all = []
    # ra, dec = frame_arr[0], frame_arr[1]
    # for radius in np.linspace(0.001, 1, 100):
    #     dist = np.sqrt((ra - center[0])**2 + (dec - center[1])**2)
    #     msk = dist <= radius
    #     cl_region = frame_arr[2:, msk]
    #     N_cl_region = cl_region.shape[1]

    #     A_cl_region = np.pi*radius**2
    #     A_field = A_total - A_cl_region
    #     dens_field = N_field / A_field

    #     n_field = int(dens_field * A_cl_region)
    #     n_memb = int(N_cl_region - n_field)
    #     if n_memb < 0:
    #         break

    #     n_all.append([n_memb, n_field])
    #     print(round(radius, 3), N_cl_region, n_field, n_memb)
    # breakpoint()

    return n_memb


def ripley_nmembs(
    x,
    y,
    pmRA,
    pmDE,
    plx,
    xy_c,
    vpd_c,
    plx_c,
    N_clust=50,
    N_extra=5,
    N_step=10,
    N_break=5,
    # dims_norm,
    # centers_ex={}
):
    """
    Estimate the number of cluster members
    """
    # extra_cls_dict, dims_msk = prep_extra_cl_dict(centers_ex)

    rads, Kest, C_thresh_N = init_ripley(x, y)

    # idx_survived = ripley_survive(
    #     lon,
    #     lat,
    #     pmRA,
    #     pmDE,
    #     plx,
    #     xy_c,
    #     vpd_c,
    #     plx_c,
    #     N_clust,
    #     N_extra,
    #     N_step,
    #     N_break,
    #     C_thresh_N,
    #     rads,
    #     Kest,
    #     # extra_cls_dict,
    #     # dims_msk,
    #     # dims_norm,
    # )
    # N_survived = len(idx_survived)

    cents_3d = np.array([list(vpd_c) + [plx_c]])
    data_3d = np.array([pmRA, pmDE, plx]).T
    # Ordered indexes according to smallest distances to 'cents_3d'
    d_pm_plx_idxs = cp.get_Nd_dists(cents_3d, data_3d)
    xy = np.array([x, y]).T
    N_stars = xy.shape[0]

    def core(N_clust_surv):
        """
        This is the core function that estimates the number of members
        based on Ripley's K-function.
        """
        # Define K-function threshold
        C_thresh = C_thresh_N / N_clust_surv
        idx_survived = []
        N_break_count, step_old = 0, 0
        for step in np.arange(N_clust_surv, N_stars, N_clust_surv):
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

    # Select those clusters where the stars are different enough from a
    # random distribution
    idx_survived = core(N_clust)
    if not idx_survived:
        # If the default clustering number did not work, try a few
        # more values with an increasing number of cluster stars
        for _ in range(N_extra):
            N_clust_surv = int(N_clust + (_ + 1) * N_step)
            idx_survived = core(N_clust_surv)
            # Break out when (if) any value selected stars
            if len(idx_survived) > 0:
                break

    # if N_survived < N_clust_min:
    #     warnings.warn("The estimated number of cluster members is " + f"<{N_clust_min}")
    #     return N_clust_min

    # Filter extra clusters in frame (if any)
    # msk = np.array(idx_survived)
    # idx_survived = filter_cls_in_frame(
    #     lon[msk],
    #     lat[msk],
    #     pmRA[msk],
    #     pmDE[msk],
    #     plx[msk],
    #     xy_c,
    #     vpd_c,
    #     plx_c,
    #     msk,
    #     extra_cls_dict,
    #     dims_msk,
    #     N_clust_min,
    #     dims_norm,
    # )

    # DONT REMEMBER WHY I REMOVED THIS BLOCK, POOR PERFORMANCE MOST LIKELY. NEED TO TEST
    # # Filter by (lon, lat) KDE
    # kde_probs = self.kde_probs(lon, lat, idx_survived, msk)
    # if kde_probs is not None:
    #     kde_prob_cut = np.percentile(kde_probs, 25)
    #     msk = kde_probs > kde_prob_cut
    #     idx_survived = idx_survived[msk]

    N_survived = len(idx_survived)

    # if N_survived > N_clust_max:
    #     warnings.warn("The estimated number of cluster members is " + f">{N_clust_max}")
    #     # # Select the maximum number of stars from those closest to the
    #     # # center
    #     # data_norm, cents_norm = mp.get_dims_norm(
    #     #     N_clust_min,
    #     #     dims_norm,
    #     #     lon,
    #     #     lat,
    #     #     pmRA,
    #     #     pmDE,
    #     #     plx,
    #     #     xy_c,
    #     #     vpd_c,
    #     #     plx_c,
    #     #     idx_survived,
    #     # )
    #     # d_idxs = cp.get_Nd_dists(cents_norm, data_norm[idx_survived])
    #     # idx_survived = idx_survived[d_idxs][:N_clust_max]
    #     return N_clust_max

    return N_survived


def init_ripley(lon, lat):
    """
    https://rdrr.io/cran/spatstat/man/Kest.html
    "For a rectangular window it is prudent to restrict the r values to a
    maximum of 1/4 of the smaller side length of the rectangle
    (Ripley, 1977, 1988; Diggle, 1983)"
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


def rkfunc(xy, rads, Kest):
    """
    Test how similar this cluster's (x, y) distribution is compared
    to a uniform random distribution using Ripley's K.
    https://stats.stackexchange.com/a/122816/10416

    Parameters
    ----------
    xy : TYPE
        Description
    rads : TYPE
        Description
    Kest : TYPE
        Description

    Returns
    -------
    TYPE
        Description
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

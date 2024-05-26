import warnings
import numpy as np
from astropy.stats import RipleysKEstimator
from scipy import spatial

# from scipy.stats import gaussian_kde # NEEDS TEST, 05/24
from . import cluster_priv as cp


def fastMP(
    X,
    xy_c,
    vpd_c,
    plx_c,
    fixed_centers,
    N_cluster,
    N_clust_min,
    N_clust_max,
    centers_ex,
    N_resample,
):
    """ """
    # HARDCODED
    N_break = max(50, int(N_resample * 0.05))
    # HARDCODED

    # Prepare dictionary of parameters for extra clusters in frame (if any)
    extra_cls_dict, dims_msk = prep_extra_cl_dict(centers_ex)

    # Remove 'nan' values
    N_all = X.shape[1]
    idx_clean, X_no_nan = cp.reject_nans(X)
    # Unpack input data with no 'nans'
    lon, lat, pmRA, pmDE, plx, e_pmRA, e_pmDE, e_plx = X_no_nan

    # Estimate initial center
    xy_c, vpd_c, plx_c = get_center(
        xy_c,
        vpd_c,
        plx_c,
        fixed_centers,
        N_cluster,
        N_clust_min,
        lon,
        lat,
        pmRA,
        pmDE,
        plx,
    )
    cents_init = [xy_c, vpd_c, plx_c]

    # Remove the most obvious field stars to speed up the process
    idx_clean, lon, lat, pmRA, pmDE, plx, e_pmRA, e_pmDE, e_plx = first_filter(
        N_clust_max,
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
    )

    # Prepare Ripley's K data
    rads, Kest, C_thresh_N = init_ripley(lon, lat)

    # Initiate here as None, will be estimated after the members estimation
    # using a selection of probable members
    dims_norm = None

    # Estimate the number of members
    if N_cluster is None:
        N_survived, st_idx = estimate_nmembs(
            N_clust_min,
            N_clust_max,
            extra_cls_dict,
            dims_msk,
            lon,
            lat,
            pmRA,
            pmDE,
            plx,
            xy_c,
            vpd_c,
            plx_c,
            C_thresh_N,
            rads,
            Kest,
            dims_norm,
        )
    else:
        N_survived = int(N_cluster)
        st_idx = None

        cents_3d = np.array([list(vpd_c) + [plx_c]])
        data_3d = np.array([pmRA, pmDE, plx]).T
        # Ordered indexes according to smallest distances to 'cents_3d'
        d_pm_plx_idxs = cp.get_Nd_dists(cents_3d, data_3d)
        st_idx = d_pm_plx_idxs[:N_survived]

    # Define here 'dims_norm' value used for data normalization
    data_5d = np.array([lon, lat, pmRA, pmDE, plx]).T
    cents_5d = np.array([xy_c + vpd_c + [plx_c]])
    data_mvd = data_5d - cents_5d
    if st_idx is not None:
        dims_norm = 2 * np.nanmedian(abs(data_mvd[st_idx]), 0)
    else:
        dims_norm = 2 * np.nanmedian(abs(data_mvd), 0)

    idx_selected = []
    N_runs, N_05_old, prob_old, break_check = 0, 0, 1, 0
    for _ in range(N_resample + 1):
        # Sample data
        s_pmRA, s_pmDE, s_plx = data_sample(pmRA, pmDE, plx, e_pmRA, e_pmDE, e_plx)

        # Data normalization
        data_5d, cents_5d = get_dims_norm(
            N_clust_min,
            dims_norm,
            lon,
            lat,
            s_pmRA,
            s_pmDE,
            s_plx,
            xy_c,
            vpd_c,
            plx_c,
            st_idx,
        )

        # Indexes of the sorted 5D distances to the estimated center
        d_idxs = cp.get_Nd_dists(cents_5d, data_5d)

        # Star selection
        st_idx = d_idxs[:N_survived]

        # Filter extra clusters in frame (if any)
        st_idx = filter_cls_in_frame(
            lon[st_idx],
            lat[st_idx],
            pmRA[st_idx],
            pmDE[st_idx],
            plx[st_idx],
            xy_c,
            vpd_c,
            plx_c,
            st_idx,
            extra_cls_dict,
            dims_msk,
            N_clust_min,
            dims_norm,
        )

        # Re-estimate centers using the selected stars
        if len(st_idx) > N_clust_min:
            xy_c, vpd_c, plx_c = get_center(
                xy_c,
                vpd_c,
                plx_c,
                fixed_centers,
                N_cluster,
                N_clust_min,
                lon[st_idx],
                lat[st_idx],
                pmRA[st_idx],
                pmDE[st_idx],
                plx[st_idx],
            )

            idx_selected += list(st_idx)
            N_runs += 1

        # Convergence check
        if idx_selected:
            break_check, prob_old, N_05_old = get_break_check(
                break_check, N_runs, idx_selected, prob_old, N_05_old
            )
            if break_check > N_break:
                break

    # Assign final probabilities
    probs_final = assign_probs(N_all, idx_clean, idx_selected, N_runs)
    # Change '0' probabilities using linear relation
    probs_final = probs_0(N_clust_min, dims_norm, X, cents_init, probs_final)

    if break_check > N_break:
        print(f"Convergence reached at {N_runs} runs.")
    else:
        print(f"Maximum number of runs reached: {N_resample}.")
    print(f"Estimated number of members: {N_survived}")

    return probs_final


def get_break_check(break_check, N_runs, idx_selected, prob_old, N_05_old):
    """ """
    counts = np.unique(idx_selected, return_counts=True)[1]
    probs = counts / N_runs
    msk = probs > 0.5  # HARDCODED
    N_05 = msk.sum()
    if N_05 > 2:
        prob_mean = np.mean(probs[msk])
        delta_probs = abs(prob_mean - prob_old)
        if N_05 == N_05_old or delta_probs < 0.001:  # HARDCODED
            break_check += 1
        else:
            # Reset
            break_check = 0
        prob_old, N_05_old = prob_mean, N_05

    return break_check, prob_old, N_05_old


def prep_extra_cl_dict(centers_ex):
    """
    The parameter 'centers_ex' must be a list of dictionaries, one
    dictionary for each extra cluster in the frame. The dictionaries
    must have at most three keys, 'xy', 'pms', 'plx', each with a list
    containing the center value(s) in those dimensions.

    Examples:

    centers_ex = [{'xy': [105.39, 0.9], 'plx': [1.3]}]

    centers_ex = [
        {'xy': [105.39, 0.9], 'pms': [3.5, -0.7]},
        {'xy': [0.82, -4.5], 'pms': [3.5, -0.7], 'plx': [3.5]}
    ]
    """
    extra_cls_dict = {"run_flag": False}

    if centers_ex is None:
        return extra_cls_dict, None

    extra_cls_dict["run_flag"] = True

    # Set the distribution of dimensions
    dims_msk = {
        "xy": np.array([0, 1]),
        "pms": np.array([2, 3]),
        "plx": np.array([4]),
        "xy_pms": np.array([0, 1, 2, 3]),
        "xy_plx": np.array([0, 1, 4]),
        "pms_plx": np.array([2, 3, 4]),
        "xy_pms_plx": np.array([0, 1, 2, 3, 4]),
    }

    # Read centers of the extra clusters in frame
    dims = ["xy", "pms", "plx"]
    dim_keys, cents = [], []
    for extra_cl in centers_ex:
        dims_ex = extra_cl.keys()
        dims_id, centers = "", []
        for d in dims:
            if d in dims_ex:
                center = extra_cl[d]
                if center is not None:
                    dims_id += "_" + d
                    centers += center
        # Store dimensions ids and centers for each extra cluster
        dim_keys.append(dims_id[1:])
        cents.append(centers)

    extra_cls_dict["dim_keys"] = dim_keys
    extra_cls_dict["centers"] = cents

    return extra_cls_dict, dims_msk


def get_center(
    xy_c,
    vpd_c,
    plx_c,
    fixed_centers,
    N_cluster,
    N_clust_min,
    lon,
    lat,
    pmRA,
    pmDE,
    plx,
):
    """Return 5-dimensional center values, within the given constrains"""

    # Skip process if all centers are fixed
    if (
        fixed_centers is True
        and xy_c is not None
        and vpd_c is not None
        and plx_c is not None
    ):
        x_c, y_c = xy_c
        pmra_c, pmde_c = vpd_c
        plx_c = plx_c
        return [x_c, y_c], [pmra_c, pmde_c], plx_c

    # Estimate initial center
    x_c, y_c, pmra_c, pmde_c, plx_c = cp.get_5D_center(
        lon, lat, pmRA, pmDE, plx, xy_c, vpd_c, plx_c, N_cluster, N_clust_min
    )
    # Re-write values if necessary
    if fixed_centers is True:
        if xy_c is not None:
            x_c, y_c = xy_c
        if vpd_c is not None:
            pmra_c, pmde_c = vpd_c
        if plx_c is not None:
            plx_c = plx_c

    return [x_c, y_c], [pmra_c, pmde_c], plx_c


def first_filter(
    N_clust_max,
    idx_all,
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
    plx_cut=0.5,
    v_kms_max=5,
    pm_max=3,
    N_times=5,
):
    """
    plx_cut: Parallax value that determines the cut in the filter rule
    between using 'v_kms_max' or 'pm_max'
    v_kms_max:
    pm_max:
    N_times: number of times used to multiply 'N_clust_max' to determine
    how many stars to keep
    """
    # Remove obvious field stars
    pms = np.array([pmRA, pmDE]).T
    pm_rad = spatial.distance.cdist(pms, np.array([vpd_c])).T[0]
    msk1 = (plx > plx_cut) & (pm_rad / (abs(plx) + 0.0001) < v_kms_max)
    msk2 = (plx <= plx_cut) & (pm_rad < pm_max)
    # Stars to keep
    msk = msk1 | msk2

    # Update arrays
    lon, lat, pmRA, pmDE, plx = lon[msk], lat[msk], pmRA[msk], pmDE[msk], plx[msk]
    e_pmRA, e_pmDE, e_plx = e_pmRA[msk], e_pmDE[msk], e_plx[msk]
    # Update indexes of surviving elements
    idx_all = idx_all[msk]

    if msk.sum() < N_clust_max * N_times:
        return idx_all, lon, lat, pmRA, pmDE, plx, e_pmRA, e_pmDE, e_plx

    # Sorted indexes of distances to pms+plx center
    cents_3d = np.array([list(vpd_c) + [plx_c]])
    data_3d = np.array([pmRA, pmDE, plx]).T
    d_pm_plx_idxs = cp.get_Nd_dists(cents_3d, data_3d)
    # Indexes of stars to keep and reject based on their distance
    idx_acpt = d_pm_plx_idxs[: int(N_clust_max * N_times)]

    # Update arrays
    lon, lat, pmRA, pmDE, plx = (
        lon[idx_acpt],
        lat[idx_acpt],
        pmRA[idx_acpt],
        pmDE[idx_acpt],
        plx[idx_acpt],
    )
    e_pmRA, e_pmDE, e_plx = e_pmRA[idx_acpt], e_pmDE[idx_acpt], e_plx[idx_acpt]
    # Update indexes of surviving elements
    idx_all = idx_all[idx_acpt]

    return idx_all, lon, lat, pmRA, pmDE, plx, e_pmRA, e_pmDE, e_plx


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


def estimate_nmembs(
    N_clust_min,
    N_clust_max,
    extra_cls_dict,
    dims_msk,
    lon,
    lat,
    pmRA,
    pmDE,
    plx,
    xy_c,
    vpd_c,
    plx_c,
    C_thresh_N,
    rads,
    Kest,
    dims_norm,
    kde_prob_cut=0.25,
    N_clust=50,
    N_extra=5,
    N_step=10,
    N_break=5,
):
    """
    Estimate the number of cluster members
    """
    idx_survived = ripley_survive(
        lon,
        lat,
        pmRA,
        pmDE,
        plx,
        xy_c,
        vpd_c,
        plx_c,
        N_clust,
        N_extra,
        N_step,
        N_break,
        C_thresh_N,
        rads,
        Kest,
        extra_cls_dict,
        dims_msk,
        N_clust_min,
        dims_norm,
    )
    N_survived = len(idx_survived)

    if N_survived < N_clust_min:
        warnings.warn("The estimated number of cluster members is " + f"<{N_clust_min}")
        return N_clust_min, None

    # Filter extra clusters in frame (if any)
    msk = np.array(idx_survived)
    idx_survived = filter_cls_in_frame(
        lon[msk],
        lat[msk],
        pmRA[msk],
        pmDE[msk],
        plx[msk],
        xy_c,
        vpd_c,
        plx_c,
        msk,
        extra_cls_dict,
        dims_msk,
        N_clust_min,
        dims_norm,
    )

    # NOT SURE WHY I REMOVED THIS BLOCK, POOR PERFORMANCE MOST LIKELY. NEED TO TEST
    # # Filter by (lon, lat) KDE
    # kde_probs = self.kde_probs(lon, lat, idx_survived, msk)
    # if kde_probs is not None:
    #     kde_prob_cut = np.percentile(kde_probs, 25)
    #     msk = kde_probs > kde_prob_cut
    #     idx_survived = idx_survived[msk]

    N_survived = len(idx_survived)

    if N_survived < N_clust_min:
        warnings.warn("The estimated number of cluster members is " + f"<{N_clust_min}")
        return N_clust_min, None

    if N_survived > N_clust_max:
        warnings.warn("The estimated number of cluster members is " + f">{N_clust_max}")
        # Select the maximum number of stars from those closest to the
        # center
        data_norm, cents_norm = get_dims_norm(
            N_clust_min,
            dims_norm,
            lon,
            lat,
            pmRA,
            pmDE,
            plx,
            xy_c,
            vpd_c,
            plx_c,
            idx_survived,
        )
        d_idxs = cp.get_Nd_dists(cents_norm, data_norm[idx_survived])
        idx_survived = idx_survived[d_idxs][:N_clust_max]
        return N_clust_max, idx_survived

    return N_survived, idx_survived


def ripley_survive(
    lon,
    lat,
    pmRA,
    pmDE,
    plx,
    xy_c,
    vpd_c,
    plx_c,
    N_clust,
    N_extra,
    N_step,
    N_break,
    C_thresh_N,
    rads,
    Kest,
    extra_cls_dict,
    dims_msk,
    N_clust_min,
    dims_norm,
):
    """
    Process data to identify the indexes of stars that survive the
    Ripley's K filter
    """
    cents_3d = np.array([list(vpd_c) + [plx_c]])
    data_3d = np.array([pmRA, pmDE, plx]).T
    # Ordered indexes according to smallest distances to 'cents_3d'
    d_pm_plx_idxs = cp.get_Nd_dists(cents_3d, data_3d)
    xy = np.array([lon, lat]).T
    N_stars = xy.shape[0]

    # Identify stars associated to extra clusters in frame (if any)
    msk = np.arange(0, len(lon))
    idx_survived_init = filter_cls_in_frame(
        lon,
        lat,
        pmRA,
        pmDE,
        plx,
        xy_c,
        vpd_c,
        plx_c,
        msk,
        extra_cls_dict,
        dims_msk,
        N_clust_min,
        dims_norm,
    )

    def core(N_clust_surv, idx_survived_init):
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

            # Remove stars associated to other clusters
            msk_extra_cls = np.isin(msk_ring, idx_survived_init)
            msk = msk_ring[msk_extra_cls]
            # If only a small percentage survived, discard
            if len(msk) < int(N_clust_surv * 0.25):  # HARDCODED
                C_s = 0

            if not np.isnan(C_s):
                # This group of stars survived
                if C_s >= C_thresh:
                    idx_survived += list(msk)
                else:
                    # Increase break condition
                    N_break_count += 1
            if N_break_count > N_break:
                break
            step_old = step
        return idx_survived

    # Select those clusters where the stars are different enough from a
    # random distribution
    idx_survived = core(N_clust, idx_survived_init)
    if not idx_survived:
        # If the default clustering number did not work, try a few
        # more values with an increasing number of cluster stars
        for _ in range(N_extra):
            N_clust_surv = int(N_clust + (_ + 1) * N_step)
            idx_survived = core(N_clust_surv, idx_survived_init)
            # Break out when (if) any value selected stars
            if len(idx_survived) > 0:
                break

    return idx_survived


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


def data_sample(pmRA, pmDE, plx, e_pmRA, e_pmDE, e_plx):
    """Gaussian random sample"""
    data_3 = np.array([pmRA, pmDE, plx])
    grs = np.random.normal(0.0, 1.0, data_3.shape[1])
    data_err = np.array([e_pmRA, e_pmDE, e_plx])
    return data_3 + grs * data_err


def get_dims_norm(
    N_clust_min, dims_norm, lon, lat, pmRA, pmDE, plx, xy_c, vpd_c, plx_c, msk
):
    """
    Normalize dimensions using twice the median of the selected probable
    members.
    """
    data_5d = np.array([lon, lat, pmRA, pmDE, plx]).T
    cents_5d = np.array([xy_c + vpd_c + [plx_c]])

    if msk is None or len(msk) < N_clust_min:
        return data_5d, cents_5d

    data_mvd = data_5d - cents_5d

    if dims_norm is None:
        dims_norm = 2 * np.nanmedian(abs(data_mvd[msk]), 0)
        data_norm = data_mvd / dims_norm
        # from sklearn.preprocessing import RobustScaler
        # data_norm = RobustScaler().fit(data_5d).transform(data_5d)
    else:
        data_norm = data_mvd / dims_norm

    cents_norm = np.array([[0.0, 0.0, 0.0, 0.0, 0.0]])

    return data_norm, cents_norm


def filter_cls_in_frame(
    lon,
    lat,
    pmRA,
    pmDE,
    plx,
    xy_c,
    vpd_c,
    plx_c,
    idx_survived,
    extra_cls_dict,
    dims_msk,
    N_clust_min,
    dims_norm,
):
    """
    Filter extra clusters in frame (if any)
    """
    # If there are no extra clusters to remove, skip
    if extra_cls_dict["run_flag"] is False:
        return idx_survived

    # Add centers to the end of these lists to normalize their values too
    dims_plus_cents = [list(_) for _ in (lon, lat, pmRA, pmDE, plx)]
    for cent in extra_cls_dict["centers"]:
        # For each dimension, some won't be present
        for i in range(5):
            try:
                dims_plus_cents[i].append(cent[i])
            except IndexError:
                dims_plus_cents[i].append(np.nan)
    lon2, lat2, pmRA2, pmDE2, plx2 = dims_plus_cents

    # Data normalization
    msk = np.full(len(lon2), True)
    data, cents = get_dims_norm(
        N_clust_min, dims_norm, lon2, lat2, pmRA2, pmDE2, plx2, xy_c, vpd_c, plx_c, msk
    )
    data = data.T

    # Extract normalized centers for extra clusters
    new_cents = []
    for i in range(len(extra_cls_dict["centers"]), 0, -1):
        new_cents.append([_ for _ in data[:, -i] if not np.isnan(_)])
    # Remove center values from the 'data' array
    data = data[:, : -len(new_cents)]

    # Find the distances to the center for all the combinations of data
    # dimensions in the extra clusters in frame
    dists = {}
    for dims in list(set(extra_cls_dict["dim_keys"])):
        msk = dims_msk[dims]
        dists[dims] = cp.get_Nd_dists(cents[:, msk], data[msk, :].T, True)

    # Identify stars closer to the cluster's center than the centers of
    # extra clusters
    msk_d = np.full(len(idx_survived), True)
    for i, cent_ex in enumerate(new_cents):
        dims_ex = extra_cls_dict["dim_keys"][i]
        msk = dims_msk[dims_ex]
        # Distance to this extra cluster's center
        dists_ex = cp.get_Nd_dists(np.array([cent_ex]), data[msk, :].T, True)

        # If the distance to the selected cluster center is smaller than
        # the distance to this extra cluster's center, keep the star
        msk_d = msk_d & (dists[dims_ex] <= dists_ex)

    # Never return less than N_clust_min stars
    if msk_d.sum() < N_clust_min:
        return idx_survived

    return idx_survived[msk_d]


def assign_probs(N_all, idx_clean, idx_selected, N_runs):
    """Assign final probabilities for all stars

    N_all: Total number of stars in input frame
    idx_clean: indexes of stars that survived the removal of 'nans' and
    of stars that are clear field stars
    """
    # Initial -1 probabilities for *all* stars
    probs_final = np.zeros(N_all) - 1

    # Number of processed stars (ie: not rejected as nans)
    N_stars = len(idx_clean)
    # Initial zero probabilities for the processed stars
    probs_all = np.zeros(N_stars)

    if idx_selected:
        # Estimate probabilities as averages of counts
        values, counts = np.unique(idx_selected, return_counts=True)
        probs = counts / N_runs
        # Store probabilities for processed stars
        probs_all[values] = probs
    else:
        warnings.warn("No stars were identified as possible members")

    # Assign the estimated probabilities to the processed stars
    probs_final[idx_clean] = probs_all

    return probs_final


def probs_0(N_clust_min, dims_norm, X, cents_init, probs_final, p_min=0.1):
    """
    To all stars with prob=0 assign a probability from 0 to p_min
    that follows a linear relation associated to their 5D distance to the
    initial defined center
    """
    # Stars with '0' probabilities
    msk_0 = probs_final == 0.0
    # If no stars with prob=0, nothing to do
    if msk_0.sum() == 0:
        return probs_final

    # Original full data
    lon, lat, pmRA, pmDE, plx = X[:5]
    xy_c, vpd_c, plx_c = cents_init

    # Data normalization for all the stars
    msk = np.full((len(lon)), True)
    data_5d, cents_5d = get_dims_norm(
        N_clust_min, dims_norm, lon, lat, pmRA, pmDE, plx, xy_c, vpd_c, plx_c, msk
    )
    # 5D distances to the estimated center
    dists = cp.get_Nd_dists(cents_5d, data_5d, True)

    # Select 'p_min' as the minimum probability between (0., 0.5)
    msk_0_5 = (probs_final > 0.0) & (probs_final < 0.5)
    if msk_0_5.sum() > 1:
        p_min = probs_final[msk_0_5].min()

    # Linear relation for: (0, d_max), (p_min, d_min)
    d_min, d_max = dists[msk_0].min(), dists[msk_0].max()
    m, h = (d_min - d_max) / p_min, d_max
    probs_0_v = (dists[msk_0] - h) / m

    # Assign new probabilities to 'msk_0' stars
    probs_final[msk_0] = probs_0_v

    return probs_final

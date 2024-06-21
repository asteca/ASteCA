import warnings
import numpy as np
from scipy import spatial

# from scipy.stats import gaussian_kde # NEEDS TEST, 05/24
from . import cluster_priv as cp
from . import membership_priv as mp


def fastMP(
    X,
    xy_c,
    vpd_c,
    plx_c,
    fixed_centers,
    N_cluster,
    N_clust_min,
    N_clust_max,
    # centers_ex,
    N_resample,
):
    """ """
    # HARDCODED
    N_break = max(50, int(N_resample * 0.05))
    # HARDCODED

    # Remove 'nan' values
    N_all = X.shape[1]
    idx_clean, X_no_nan = cp.reject_nans(X)
    # Unpack input data with no 'nans'
    lon, lat, pmRA, pmDE, plx, e_pmRA, e_pmDE, e_plx = X_no_nan

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

    # Initiate here as None, will be estimated after the members estimation
    # using a selection of probable members
    # dims_norm = None

    # # Estimate the number of members
    # if N_cluster is None:
    #     N_survived, st_idx = estimate_nmembs(
    # else:
    # N_survived = int(N_cluster)
    # st_idx = None

    # Select initial set of stars to estimate the 'dims_norm' parameter
    cents_3d = np.array([list(vpd_c) + [plx_c]])
    data_3d = np.array([pmRA, pmDE, plx]).T
    # Ordered indexes according to smallest distances to 'cents_3d'
    d_pm_plx_idxs = cp.get_Nd_dists(cents_3d, data_3d)
    st_idx = d_pm_plx_idxs[:N_cluster]
    # Define here 'dims_norm' value used for data normalization
    data_5d = np.array([lon, lat, pmRA, pmDE, plx]).T
    cents_5d = np.array([xy_c + vpd_c + [plx_c]])
    data_mvd = data_5d - cents_5d
    # if st_idx is not None:
    dims_norm = 2 * np.nanmedian(abs(data_mvd[st_idx]), 0)
    # else:
    #     dims_norm = 2 * np.nanmedian(abs(data_mvd), 0)

    idx_selected = []
    N_runs, N_05_old, prob_old, break_check = 0, 0, 1, 0
    for _ in range(N_resample + 1):
        # Sample data
        s_pmRA, s_pmDE, s_plx = data_sample(pmRA, pmDE, plx, e_pmRA, e_pmDE, e_plx)

        # Data normalization
        data_5d, cents_5d = mp.get_dims_norm(
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
        st_idx = d_idxs[:N_cluster]

        # # Filter extra clusters in frame (if any)
        # st_idx = filter_cls_in_frame(
        #     lon[st_idx],
        #     lat[st_idx],
        #     pmRA[st_idx],
        #     pmDE[st_idx],
        #     plx[st_idx],
        #     xy_c,
        #     vpd_c,
        #     plx_c,
        #     st_idx,
        #     extra_cls_dict,
        #     dims_msk,
        #     N_clust_min,
        #     dims_norm,
        # )

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
        print(f"Convergence reached at {N_runs} runs")
    else:
        print(f"Maximum number of runs reached: {N_resample}")

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


def data_sample(pmRA, pmDE, plx, e_pmRA, e_pmDE, e_plx):
    """Gaussian random sample"""
    data_3 = np.array([pmRA, pmDE, plx])
    grs = np.random.normal(0.0, 1.0, data_3.shape[1])
    data_err = np.array([e_pmRA, e_pmDE, e_plx])
    return data_3 + grs * data_err


def assign_probs(N_all, idx_clean, idx_selected, N_runs):
    """Assign final probabilities for all stars

    N_all: Total number of stars in input frame
    idx_clean: indexes of stars that survived the removal of 'nans' and
    of stars that are clear field stars
    """
    # Initial 0 probabilities for *all* stars
    probs_final = np.zeros(N_all)

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
    data_5d, cents_5d = mp.get_dims_norm(
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

import numpy as np
from scipy import spatial
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
    rng,
    N_resample,
):
    """ """

    # Remove 'nan' values
    N_all = X.shape[1]
    idx_clean, X_no_nan = cp.reject_nans(X)
    # Unpack input data with no 'nans'
    lon, lat, pmRA, pmDE, plx, e_pmRA, e_pmDE, e_plx = X_no_nan

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

    st_idx = None
    cents_5d = np.array([[0.0, 0.0, 0.0, 0.0, 0.0]])
    N_stars = len(idx_clean)
    probs_all = np.zeros(N_stars)
    prob_old_arr = np.zeros(N_stars)
    N_break = 50
    for r in range(N_resample):
        # Sample data
        s_pmRA, s_pmDE, s_plx = data_sample(rng, pmRA, pmDE, plx, e_pmRA, e_pmDE, e_plx)

        # Data normalization
        data_5d = get_dims_norm(
            N_cluster,
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

        # Re-estimate centers using the selected stars
        xy_c, vpd_c, plx_c = get_center(
            xy_c,
            vpd_c,
            plx_c,
            fixed_centers,
            N_clust_min,
            N_clust_max,
            lon[st_idx],
            lat[st_idx],
            pmRA[st_idx],
            pmDE[st_idx],
            plx[st_idx],
        )

        probs_all[st_idx] += 1
        probs = probs_all / (r + 1)
        msk = probs > 0.5
        # Check that all P>0.5 probabilities converged to 1%
        if (abs(prob_old_arr[msk] - probs[msk]) < 0.01).all() and r > N_break:
            break
        else:
            prob_old_arr = np.array(probs)

    if r < N_resample:
        print(f"Convergence reached at {r+1} runs")
    else:
        print(f"Maximum number of runs reached: {N_resample}")

    probs_final = np.zeros(N_all)
    probs_final[idx_clean] = probs

    return probs_final


def get_dims_norm(N_cluster, lon, lat, pmRA, pmDE, plx, xy_c, vpd_c, plx_c, st_idx):
    """
    Normalize dimensions using twice the median of the selected probable
    members.
    """
    data_5d = np.array([lon, lat, pmRA, pmDE, plx]).T
    cents_5d = np.array([xy_c + vpd_c + [plx_c]])
    data_mvd = data_5d - cents_5d

    if st_idx is None:
        # Initial 'dims_norm' estimation
        cents_3d = np.array([list(vpd_c) + [plx_c]])
        data_3d = np.array([pmRA, pmDE, plx]).T
        # Ordered indexes according to smallest distances to 'cents_3d'
        d_pm_plx_idxs = cp.get_Nd_dists(cents_3d, data_3d)
        st_idx = d_pm_plx_idxs[:N_cluster]

    # This is the old way of normalizing the dimensions. Does not work well when the
    # frame is not square (e.g.: when the declination is very large and transforming
    # the frame to galactic coordinates visibly rotates it)
    # dims_norm = 2 * np.nanmedian(abs(data_mvd[st_idx]), 0)

    # Use the IQR
    dims_norm = np.ptp(np.nanpercentile(data_mvd[st_idx, :], [25, 75], axis=0), axis=0)

    data_norm = data_mvd / dims_norm
    return data_norm


def get_center(
    xy_c,
    vpd_c,
    plx_c,
    fixed_centers,
    N_clust_min,
    N_clust_max,
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
        lon, lat, pmRA, pmDE, plx, xy_c, vpd_c, plx_c, N_clust_min, N_clust_max
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


def data_sample(rng, pmRA, pmDE, plx, e_pmRA, e_pmDE, e_plx):
    """Gaussian random sample"""
    data_3 = np.array([pmRA, pmDE, plx])
    grs = rng.normal(0.0, 1.0, data_3.shape[1])
    data_err = np.array([e_pmRA, e_pmDE, e_plx])
    return data_3 + grs * data_err

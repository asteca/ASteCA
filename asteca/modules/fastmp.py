import numpy as np
from scipy.stats import iqr

from . import cluster_priv as cp


def fastMP(
    rng: np.random.Generator,
    xy_c: tuple[float, float],
    vpd_c: tuple[float, float],
    plx_c: float,
    N_cluster: int,
    fixed_centers: bool,
    X: np.ndarray,
    N_resample: int,
) -> tuple[str, np.ndarray]:
    """Perform a fast iterative Monte Carlo process to identify cluster members.

    :param rng: Random number generator for resampling.
    :type rng: np.random.Generator
    :param xy_c: Initial center coordinates (longitude, latitude).
    :type xy_c: tuple[float, float]
    :param vpd_c: Initial center proper motion values (pmRA, pmDE).
    :type vpd_c: tuple[float, float]
    :param plx_c: Initial center parallax value.
    :type plx_c: float
    :param N_cluster: Number of stars to select for cluster identification.
    :type N_cluster: int
    :param fixed_centers: If True, keep the centers fixed during the iterative process.
    :type fixed_centers: bool
    :param X: Input data array containing coordinates, proper motions, and parallax
        values.
    :type X: np.ndarray
    :param N_resample: Number of resampling iterations.
    :type N_resample: int

    :returns: A tuple containing:
        - The output message (str)
        - Array of final membership probabilities for each star (np.ndarray)
    :rtype: tuple[str, np.ndarray]
    """

    # Remove 'nan' values
    N_all = X.shape[1]
    idx_clean, X_no_nan = cp.reject_nans(X)
    # Unpack input data with no 'nans'
    lon, lat, pmRA, pmDE, plx, e_pmRA, e_pmDE, e_plx = X_no_nan

    # Remove the most obvious field stars. This has no effect on the results but
    # speeds up the process enormously
    N_filter_max = min(N_cluster * 10, 10_000)  # FIXED VALUE
    idx_clean, lon, lat, pmRA, pmDE, plx, e_pmRA, e_pmDE, e_plx = cp.first_filter(
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

    # Use the IQR to estimate the normalization constants
    data_iqr = iqr([lon, lat, pmRA, pmDE, plx], 1)

    st_idx = None
    N_stars = len(idx_clean)
    probs_all = np.zeros(N_stars)
    prob_old_arr = np.zeros(N_stars)
    data_err = np.array([e_pmRA, e_pmDE, e_plx])
    N_break = 50
    r, probs = 0, []
    for r in range(N_resample):
        # Sample data in parallax and proper motions
        grs = rng.normal(0.0, 1.0, N_stars)
        s_pmRA, s_pmDE, s_plx = np.array([pmRA, pmDE, plx]) + grs * data_err

        # Normalize sampled data. Move it to 0 center first otherwise the
        # 'get_Nd_dists' call will fail
        data_5d = np.array([lon, lat, s_pmRA, s_pmDE, s_plx]).T - np.array(
            [list(xy_c) + list(vpd_c) + [plx_c]]
        )
        data_5d /= data_iqr

        # Indexes of the sorted 5D distances to the (moved) centers
        cents_5d = np.array([[0.0, 0.0, 0.0, 0.0, 0.0]])
        d_idxs = cp.get_Nd_dists(cents_5d, data_5d)

        # Star selection
        st_idx = d_idxs[:N_cluster]

        # Only re-estimate the IQR for the first set of selected members
        if r == 0:
            data_iqr = iqr(
                [lon[st_idx], lat[st_idx], pmRA[st_idx], pmDE[st_idx], plx[st_idx]], 1
            )

        if fixed_centers is False:
            # Re-estimate centers using the selected stars
            x_c, y_c, pmra_c, pmde_c, plx_c = cp.get_knn_5D_center(
                lon[st_idx],
                lat[st_idx],
                pmRA[st_idx],
                pmDE[st_idx],
                plx[st_idx],
                xy_c,
                vpd_c,
                plx_c,
            )
            xy_c, vpd_c = (x_c, y_c), (pmra_c, pmde_c)

        probs_all[st_idx] += 1
        probs = probs_all / (r + 1)
        msk = probs > 0.5
        # Check that all P>0.5 probabilities converged to 1%
        if (abs(prob_old_arr[msk] - probs[msk]) < 0.01).all() and r > N_break:
            break
        else:
            prob_old_arr = np.array(probs)

    if r < N_resample:
        out_mssg = f"Convergence reached at {r + 1} runs"
    else:
        out_mssg = f"Maximum number of runs reached: {N_resample}"

    probs_final = np.zeros(N_all)
    probs_final[idx_clean] = probs

    return out_mssg, probs_final

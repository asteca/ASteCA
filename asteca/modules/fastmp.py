import numpy as np

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

    st_idx = None
    cents_5d = np.array([[0.0, 0.0, 0.0, 0.0, 0.0]])
    N_stars = len(idx_clean)
    probs_all = np.zeros(N_stars)
    prob_old_arr = np.zeros(N_stars)
    N_break = 50
    r, probs = 0, []
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

        if fixed_centers is False:
            # Re-estimate centers using the selected stars
            x_c, y_c, pmra_c, pmde_c, plx_c = cp.get_knn_5D_center(
                lon, lat, pmRA, pmDE, plx, xy_c, vpd_c, plx_c
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


def get_dims_norm(
    N_cluster: int,
    lon: np.ndarray,
    lat: np.ndarray,
    pmRA: np.ndarray,
    pmDE: np.ndarray,
    plx: np.ndarray,
    xy_c: tuple[float, float],
    vpd_c: tuple[float, float],
    plx_c: float,
    st_idx: np.ndarray | None,
) -> np.ndarray:
    """Normalize spatial, proper motion, and parallax data using the interquartile
    range of the selected probable members.

    :param N_cluster: Number of stars to select for cluster identification. Only used
     when st_idx is None the first time the for block is run
    :type N_cluster: int
    :param lon: Longitude values.
    :type lon: np.ndarray
    :param lat: Latitude values.
    :type lat: np.ndarray
    :param pmRA: Proper motions in right ascension.
    :type pmRA: np.ndarray
    :param pmDE: Proper motions in declination.
    :type pmDE: np.ndarray
    :param plx: Parallax values.
    :type plx: np.ndarray
    :param xy_c: Center coordinates (longitude, latitude).
    :type xy_c: tuple[float, float]
    :param vpd_c: Center proper motion values (pmRA, pmDE).
    :type vpd_c: tuple[float, float]
    :param plx_c: Center parallax value.
    :type plx_c: float
    :param st_idx: Indices of the selected stars.
    :type st_idx: np.ndarray | None

    :returns: Normalized 5D data.
    :rtype: np.ndarray
    """
    data_5d = np.array([lon, lat, pmRA, pmDE, plx]).T
    cents_5d = np.array([list(xy_c) + list(vpd_c) + [plx_c]])
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

    # Use the IQR to estimate the normalization constant
    dims_norm = np.ptp(np.nanpercentile(data_mvd[st_idx, :], [25, 75], axis=0), axis=0)

    data_norm = data_mvd / dims_norm
    return data_norm


def data_sample(
    rng: np.random.Generator,
    pmRA: np.ndarray,
    pmDE: np.ndarray,
    plx: np.ndarray,
    e_pmRA: np.ndarray,
    e_pmDE: np.ndarray,
    e_plx: np.ndarray,
) -> np.ndarray:
    """Generate a Gaussian random sample of proper motions and parallax.

    :param rng: Random number generator for sampling.
    :type rng: np.random.Generator
    :param pmRA: Proper motions in right ascension.
    :type pmRA: np.ndarray
    :param pmDE: Proper motions in declination.
    :type pmDE: np.ndarray
    :param plx: Parallax values of the stars.
    :type plx: np.ndarray
    :param e_pmRA: Errors in proper motions (pmRA).
    :type e_pmRA: np.ndarray
    :param e_pmDE: Errors in proper motions (pmDE).
    :type e_pmDE: np.ndarray
    :param e_plx: Errors in parallax.
    :type e_plx: np.ndarray

    :returns: Gaussian random sampled data based on proper motions and parallax.
    :rtype: np.ndarray
    """
    data_3 = np.array([pmRA, pmDE, plx])
    grs = rng.normal(0.0, 1.0, data_3.shape[1])
    data_err = np.array([e_pmRA, e_pmDE, e_plx])
    data_err = data_3 + grs * data_err
    return data_err

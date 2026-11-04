import numpy as np
from scipy.stats import iqr

from . import cluster_priv as cp

# Maximum number of stars to consider for the first filter. This is a safeguard
# against very large datasets that would make the process too slow.
# The first filter is applied to remove the most obvious field stars, which
# has no effect on the results but speeds up the process enormously.
N_STARS_MAX = 10_000

# The 5D centers are set to zero because the data is normalized by subtracting
# the centers before computing the distances. This ensures that the distances are
# calculated relative to the center of the cluster, which is now at the origin
# (0, 0, 0, 0, 0) in the normalized space. The normalization allows for a
# consistent comparison of distances across different dimensions and scales.
CENTS_5D = np.array([[0.0, 0.0, 0.0, 0.0, 0.0]])


def fastMP(
    rng: np.random.Generator,
    xy_c: tuple[float, float],
    vpd_c: tuple[float, float],
    plx_c: float,
    N_cluster: int,
    fixed_centers: bool,
    X: np.ndarray,
    N_resample: int,
    mag: np.ndarray | None = None,
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
    :param mag: Photometric magnitudes (e.g.: Gaia's G). If given, a reference
        selection without error resampling is obtained first, and the final
        probabilities are adjusted so that a ``P>=0.5`` cut follows the magnitude
        distribution of this reference selection, see :py:func:`mag_dist_correct`.
        This removes the bias towards stars with small uncertainties (ie: brighter
        stars). If ``None``, no adjustment is applied.
    :type mag: np.ndarray | None

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
    N_filter_max = min(N_cluster * 10, N_STARS_MAX)
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

    data_err = np.array([e_pmRA, e_pmDE, e_plx])
    data = (lon, lat, pmRA, pmDE, plx, data_err)

    # Use the IQR to estimate the normalization constants
    data_iqr = iqr([lon, lat, pmRA, pmDE, plx], 1)

    out_mssg, probs = resample_loop(
        rng, xy_c, vpd_c, plx_c, N_cluster, fixed_centers, data, data_iqr, N_resample
    )

    if mag is not None:
        # Reference selection with no error resampling. This gives the magnitude
        # distribution that the members should follow
        ref_msk, ref_dist = extract_reference_stars(
            xy_c, vpd_c, plx_c, N_cluster, data, data_iqr
        )
        # Extract the magnitudes of the clean stars
        mag_clean = mag[idx_clean]
        # Adjust the probabilities to follow the magnitude distribution of the
        # reference selection
        probs = mag_dist_correct(probs, ref_msk, ref_dist, mag_clean)

    probs_final = np.zeros(N_all)
    probs_final[idx_clean] = probs

    return out_mssg, probs_final


def extract_reference_stars(
    xy_c: tuple[float, float],
    vpd_c: tuple[float, float],
    plx_c: float,
    N_cluster: int,
    data: tuple,
    data_iqr: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Extract a reference selection of ``N_cluster`` stars without error resampling.

    The ``N_cluster`` stars closest to the (fixed) centers are selected using the
    IQR of all the stars as normalization constants. The IQR is then re-estimated
    with these stars and the selection is repeated with the new IQR, as done in the
    resampling loop.

    :param xy_c: Center coordinates (longitude, latitude).
    :type xy_c: tuple[float, float]
    :param vpd_c: Center proper motion values (pmRA, pmDE).
    :type vpd_c: tuple[float, float]
    :param plx_c: Center parallax value.
    :type plx_c: float
    :param N_cluster: Number of stars to select.
    :type N_cluster: int
    :param data: Tuple with (lon, lat, pmRA, pmDE, plx, data_err); ``data_err`` is
        not used.
    :type data: tuple
    :param data_iqr: IQR of the 5 dimensions for all the stars, used as initial
        normalization constants.
    :type data_iqr: np.ndarray

    :returns: A tuple containing:
        - Mask of the stars in the reference selection (np.ndarray)
        - Normalized 5D distances of all the stars to the centers (np.ndarray)
    :rtype: tuple[np.ndarray, np.ndarray]
    """
    lon, lat, pmRA, pmDE, plx, _ = data
    N_stars = len(lon)

    # Select the N_cluster stars closest to the centers using the IQR of all the stars
    data_5d = norm_data(lon, lat, pmRA, pmDE, plx, xy_c, vpd_c, plx_c, data_iqr)
    st_idx = cp.get_Nd_dists(CENTS_5D, data_5d)[:N_cluster]

    # Re-estimate the IQR with the selected stars and select again
    data_iqr = iqr(
        [lon[st_idx], lat[st_idx], pmRA[st_idx], pmDE[st_idx], plx[st_idx]], 1
    )
    data_5d = norm_data(lon, lat, pmRA, pmDE, plx, xy_c, vpd_c, plx_c, data_iqr)
    st_idx = cp.get_Nd_dists(CENTS_5D, data_5d)[:N_cluster]

    # Compute the normalized 2D distances of all stars to the centers
    ref_dist = np.linalg.norm(data_5d, axis=1)

    # Mask of the stars in the reference selection
    ref_msk = np.zeros(N_stars, bool)
    ref_msk[st_idx] = True

    return ref_msk, ref_dist


def resample_loop(
    rng: np.random.Generator,
    xy_c: tuple[float, float],
    vpd_c: tuple[float, float],
    plx_c: float,
    N_cluster: int,
    fixed_centers: bool,
    data: tuple,
    data_iqr: np.ndarray,
    N_resample: int,
) -> tuple[str, np.ndarray]:
    """Iteratively resample the data and select the stars closest to the centers.

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
    :param data: Tuple with (lon, lat, pmRA, pmDE, plx, data_err), where
        ``data_err`` is the array of (e_pmRA, e_pmDE, e_plx) uncertainties.
    :type data: tuple
    :param data_iqr: IQR of the 5 dimensions for all the stars, used as
        normalization constants in the first iteration. After it, they are
        re-estimated with the selected stars.
    :type data_iqr: np.ndarray
    :param N_resample: Maximum number of resampling iterations.
    :type N_resample: int

    :returns: A tuple containing:
        - The output message (str)
        - Array of membership probabilities for each star (np.ndarray)
    :rtype: tuple[str, np.ndarray]
    """
    lon, lat, pmRA, pmDE, plx, data_err = data
    e_pmRA, e_pmDE, e_plx = data_err
    N_stars = len(lon)

    st_idx = None
    probs_all = np.zeros(N_stars)
    prob_old_arr = np.zeros(N_stars)
    N_break = 50
    r, probs = 0, np.zeros(N_stars)
    converged = False
    # Squared normalized (lon, lat) distances to the center. These are not
    # resampled, so they are only re-computed when the IQR or centers change
    d2_xy = None
    for r in range(N_resample):
        # Sample data in parallax and proper motions, with independent draws for
        # each dimension
        grs_pmRA, grs_pmDE, grs_plx = rng.normal(0.0, 1.0, (3, N_stars))

        if d2_xy is None:
            d2_xy = ((lon - xy_c[0]) / data_iqr[0]) ** 2 + (
                (lat - xy_c[1]) / data_iqr[1]
            ) ** 2
        # Squared normalized 5D distances to the (moved) centers. Squared
        # distances select the same stars, and avoid the square root
        d2 = (
            d2_xy
            + ((pmRA + grs_pmRA * e_pmRA - vpd_c[0]) / data_iqr[2]) ** 2
            + ((pmDE + grs_pmDE * e_pmDE - vpd_c[1]) / data_iqr[3]) ** 2
            + ((plx + grs_plx * e_plx - plx_c) / data_iqr[4]) ** 2
        )

        # Star selection: the 'N_cluster' closest stars. Their order is not
        # needed so a partition is used instead of a full sort
        if N_cluster >= N_stars:
            st_idx = np.arange(N_stars)
        else:
            st_idx = np.argpartition(d2, N_cluster - 1)[:N_cluster]

        # Only re-estimate the IQR for the first set of selected members
        if r == 0:
            data_iqr = iqr(
                [lon[st_idx], lat[st_idx], pmRA[st_idx], pmDE[st_idx], plx[st_idx]], 1
            )
            d2_xy = None

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
            d2_xy = None

        # Update the probabilities of the selected stars
        probs_all[st_idx] += 1
        # Compute the probabilities for all stars
        probs = probs_all / (r + 1)

        # Check for convergence
        msk = probs > 0.5
        # Check that all P>0.5 probabilities converged to 1%
        if (abs(prob_old_arr[msk] - probs[msk]) < 0.01).all() and r > N_break:
            converged = True
            break
        else:
            # 'probs' is a new array in every iteration, no copy required
            prob_old_arr = probs

    if converged:
        out_mssg = f"Convergence reached at {r + 1} runs"
    else:
        out_mssg = f"Maximum number of runs reached: {N_resample}"

    return out_mssg, np.array(probs)


def norm_data(
    lon: np.ndarray,
    lat: np.ndarray,
    pmRA: np.ndarray,
    pmDE: np.ndarray,
    plx: np.ndarray,
    xy_c: tuple[float, float],
    vpd_c: tuple[float, float],
    plx_c: float,
    data_iqr: np.ndarray,
) -> np.ndarray:
    """Move the 5D data to a 0 center (otherwise the 'get_Nd_dists' call will fail)
    and normalize it by the IQR.

    :param lon: Longitude of the stars.
    :type lon: np.ndarray
    :param lat: Latitude of the stars.
    :type lat: np.ndarray
    :param pmRA: Proper motions in right ascension.
    :type pmRA: np.ndarray
    :param pmDE: Proper motions in declination.
    :type pmDE: np.ndarray
    :param plx: Parallaxes.
    :type plx: np.ndarray
    :param xy_c: Center coordinates (longitude, latitude).
    :type xy_c: tuple[float, float]
    :param vpd_c: Center proper motion values (pmRA, pmDE).
    :type vpd_c: tuple[float, float]
    :param plx_c: Center parallax value.
    :type plx_c: float
    :param data_iqr: IQR of the 5 dimensions, used as normalization constants.
    :type data_iqr: np.ndarray

    :returns: Array of shape (N_stars, 5) with the normalized data.
    :rtype: np.ndarray
    """
    data_5d = np.array([lon, lat, pmRA, pmDE, plx]).T - np.array(
        [list(xy_c) + list(vpd_c) + [plx_c]]
    )
    return data_5d / data_iqr


def mag_dist_correct(
    probs: np.ndarray,
    ref_msk: np.ndarray,
    ref_dist: np.ndarray,
    mag_clean: np.ndarray,
    bin_width: float = 0.5,
    p_promote: float = 0.51,
    p_demote: float = 0.49,
) -> np.ndarray:
    """Adjust the probabilities to follow the magnitude distribution of a reference
    selection.

    Faint members have large uncertainties, so the resampling scatters them away
    from the centers and they are selected in a smaller fraction of resamples than
    bright members. As a result, a cut on the probabilities (e.g., P>=0.5)
    under-represents the faint end. The reference selection (obtained with no error
    resampling) does not suffer from this bias.

    For each magnitude bin, the number of reference stars ``n_ref`` is compared with
    the number of stars with ``P>=0.5`` ``n_sel``:

    - ``n_ref == n_sel``: no changes
    - ``n_ref > n_sel``: the ``n_ref - n_sel`` stars with ``P<0.5`` and the smallest
      reference distances are promoted to ``P=p_promote``
    - ``n_ref < n_sel``: the ``n_sel - n_ref`` stars with ``P>=0.5`` and the largest
      reference distances are demoted to ``P=p_demote``

    After this, a ``P>=0.5`` cut selects ``n_ref`` stars in every bin. The reference
    distance is the normalized 5D distance of each star to the centers, obtained
    with no error resampling. Stars with ``nan`` magnitudes are not binned (they do
    not count towards ``n_ref`` or ``n_sel``) and their probabilities are not
    modified. The ``probs`` array is modified in place.

    :param probs: Membership probabilities of the stars.
    :type probs: np.ndarray
    :param ref_msk: Mask of the stars in the reference selection.
    :type ref_msk: np.ndarray
    :param ref_dist: Reference 5D distances of the stars to the centers.
    :type ref_dist: np.ndarray
    :param mag_clean: Magnitudes of the stars.
    :type mag_clean: np.ndarray
    :param bin_width: Width of the magnitude bins.
    :type bin_width: float
    :param p_promote: Probability assigned to the promoted stars (must be ``>=0.5``).
    :type p_promote: float
    :param p_demote: Probability assigned to the demoted stars (must be ``<0.5``).
    :type p_demote: float

    :returns: Array of adjusted probabilities.
    :rtype: np.ndarray
    """
    # Only stars with valid magnitudes are binned (and possibly modified)
    idx_valid = np.flatnonzero(np.isfinite(mag_clean))
    if len(idx_valid) == 0:
        return probs
    mag_valid = mag_clean[idx_valid]

    # Bins covering all the stars with valid magnitudes
    edges = np.arange(mag_valid.min(), mag_valid.max() + bin_width, bin_width)
    if len(edges) < 2:
        edges = np.array([mag_valid.min(), mag_valid.min() + bin_width])
    N_bins = len(edges) - 1

    # Index of the bin for each star. Use np.clip to avoid out-of-bounds
    bin_idx = np.clip(np.digitize(mag_valid, edges) - 1, 0, N_bins - 1)

    # Number of reference stars (with valid magnitudes) per bin
    n_ref = np.bincount(bin_idx[ref_msk[idx_valid]], minlength=N_bins)

    for b in np.unique(bin_idx):
        i_b = idx_valid[bin_idx == b]
        sel = probs[i_b] >= 0.5
        dN = n_ref[b] - sel.sum()

        if dN > 0:
            # Promote the non-selected stars closest to the centers
            cand = i_b[~sel]
            probs[cand[np.argsort(ref_dist[cand])[:dN]]] = p_promote
        elif dN < 0:
            # Demote the selected stars farthest from the centers
            cand = i_b[sel]
            probs[cand[np.argsort(-ref_dist[cand])[:-dN]]] = p_demote

    return probs

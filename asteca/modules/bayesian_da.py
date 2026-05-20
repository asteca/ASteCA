import warnings

import numpy as np


def bayesian_mp(
    N_cluster: int,
    frame_arr: np.ndarray,
    e_frame_arr: np.ndarray,
    center: tuple[float, float],
    radius: float,
    bayesda_runs: int,
    rng: np.random.Generator,
) -> tuple[str, np.ndarray]:
    """Bayesian field decontamination algorithm.

    See Perren et al. (2015) for details.

    A larger radius appears to translate to more estimated members using the
    "density" method (which) is expected, but a smaller number of stars with
    P>0.5 (and a cleaner final sequence)

    :param N_cluster: Number of cluster members.
    :type N_cluster: int
    :param frame_arr: Array with the data.
    :type frame_arr: np.ndarray
    :param e_frame_arr: Array with the errors.
    :type e_frame_arr: np.ndarray
    :param center: Center of the cluster.
    :type center: tuple[float, float]
    :param radius: Radius of the cluster.
    :type radius: float
    :param bayesda_runs: Number of iterations.
    :type bayesda_runs: int
    :param rng: Random number generator.
    :type rng: np.random.Generator

    :raises ValueError: If no field region stars are found.

    :return: Message and array with the membership probabilities.
    :rtype: tuple[str, np.ndarray]
    """
    cl_region, e_cl_region, fl_region_all, e_fl_region_all, cl_reg_idxs = get_regions(
        frame_arr, e_frame_arr, center, radius
    )
    if fl_region_all.size == 0:
        raise ValueError(
            "No field region stars found. Consider reducing the radius or increasing"
            " the frame size."
        )

    N_cl_region = cl_region.shape[1]
    n_field = max(5, N_cl_region - N_cluster)

    arr_shape = (N_cl_region, N_cluster, cl_region.shape[0])
    size_Mb = np.prod(arr_shape) * np.dtype(np.float64).itemsize / (1024 * 1024)
    if size_Mb > 1000:
        warnings.warn(
            f"\nThe array generated will be larger than {size_Mb:.0f} Mb."
            + " Consider reducing the size\nof the frame or the radius."
        )

    # Normalize data
    cl_region, e_cl_region2 = dataNorm(cl_region, e_cl_region)
    fl_region_all, e_fl_region2_all = dataNorm(fl_region_all, e_fl_region_all)

    # Prepare cluster-region data for broadcasting: (N_cl_region, 1, n_dims)
    cl_region_T = cl_region.T[:, None]
    e_cl_region2_T = e_cl_region2.T[:, None]

    # Precompute cluster-region self-likelihood matrix (N_cl_region x N_cl_region)
    L_cl = _likelihood_matrix(cl_region_T, e_cl_region2_T, cl_region, e_cl_region2)

    prob_old_arr = np.zeros(N_cl_region)
    # Probabilities for all stars in the cluster region.
    sum_cl_probs = np.zeros(N_cl_region)

    N_break = 50
    r, probs = 0, []
    for r in range(bayesda_runs):
        # Select cluster members according to current probability estimates
        if N_cl_region > N_cluster:
            if r == 0:
                # Initial run
                p = rng.choice(N_cl_region, N_cluster, replace=False)
            else:
                p = rng.choice(
                    N_cl_region,
                    N_cluster,
                    replace=False,
                    p=sum_cl_probs / sum_cl_probs.sum(),
                )
        else:
            p = np.ones(N_cl_region, dtype=int)

        # Generate a random field region and compute its likelihood
        fl_region, e_fl_region2 = generate_field_region(
            rng, fl_region_all, e_fl_region2_all, n_field
        )
        # Compare cluster region with this field region
        fl_lkl = _likelihood_sum_chunked(
            cl_region_T, e_cl_region2_T, fl_region, e_fl_region2
        )
        # Cluster likelihood from precomputed matrix (fast sum over selected columns)
        cl_lkl = L_cl[:, p].sum(axis=1)

        # Bayesian probability for each star
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            # Bayesian probability for each star within the cluster region.
            bayes_prob = 1.0 / (1.0 + (fl_lkl / cl_lkl))

        # Replace possible nan values with 0.
        bayes_prob[np.isnan(bayes_prob)] = 0.0
        sum_cl_probs += bayes_prob

        probs = sum_cl_probs / (r + 1)
        msk = probs > 0.5
        # Check that all P>0.5 probabilities converged to 1%
        if (abs(prob_old_arr[msk] - probs[msk]) < 0.01).all() and r > N_break:
            break
        prob_old_arr = probs.copy()

    if r < bayesda_runs:
        out_mssg = f"Convergence reached at {r + 1} runs"
    else:
        out_mssg = f"Maximum number of runs reached: {bayesda_runs}"

    probs_final = np.zeros(frame_arr.shape[1])
    probs_final[cl_reg_idxs] = probs

    return out_mssg, probs_final


def get_regions(
    frame_arr: np.ndarray,
    e_frame_arr2: np.ndarray,
    center: tuple[float, float],
    radius: float,
    N_max: int = 2000,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Identify stars inside/outside the cluster region.

    :param frame_arr: Array with the data.
    :type frame_arr: np.ndarray
    :param e_frame_arr2: Array with the errors.
    :type e_frame_arr2: np.ndarray
    :param center: Center of the cluster.
    :type center: tuple[float, float]
    :param radius: Radius of the cluster.
    :type radius: float
    :param N_max: Maximum number of stars in the cluster region.
    :type N_max: int

    :return: Arrays with the cluster region, errors, field region, errors and indexes.
    :rtype: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]
    """
    ra, dec = frame_arr[0], frame_arr[1]
    dist_sq = (ra - center[0]) ** 2 + (dec - center[1]) ** 2
    radius_sq = radius**2

    msk_in_rad = dist_sq <= radius_sq
    if msk_in_rad.sum() > N_max:
        warnings.warn(
            f"\nNumber of stars within the radius ({msk_in_rad.sum()}) is greater than"
            f" the accepted maximum ({N_max}). Selecting {N_max} stars closest to center."
            " Consider reducing the radius or employing another method to improve the performance."
        )
        # Select the N_cluster stars closest to the center
        closest_idxs = np.argsort(dist_sq)[:N_max]
        msk_in_rad = np.zeros_like(msk_in_rad, dtype=bool)
        msk_in_rad[closest_idxs] = True

    cl_reg_idxs = np.where(msk_in_rad)[0]

    # Cluster region inside radius
    cl_region = frame_arr[2:, msk_in_rad]
    e_cl_region2 = e_frame_arr2[:, msk_in_rad]

    # Field region outside radius
    fl_region_all = frame_arr[2:, ~msk_in_rad]
    e_fl_region2_all = e_frame_arr2[:, ~msk_in_rad]

    return cl_region, e_cl_region2, fl_region_all, e_fl_region2_all, cl_reg_idxs


def generate_field_region(
    rng: np.random.Generator,
    fl_region_all: np.ndarray,
    e_fl_region2_all: np.ndarray,
    n_field: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Generate random field regions.

    :param rng: Random number generator.
    :type rng: np.random.Generator
    :param fl_region_all: Array with the field region data.
    :type fl_region_all: np.ndarray
    :param e_fl_region2_all: Array with the field region errors.
    :type e_fl_region2_all: np.ndarray
    :param n_field: Number of field stars.
    :type n_field: int

    :return: Arrays with the field region and errors.
    :rtype: tuple[np.ndarray, np.ndarray]
    """
    idxs = rng.choice(fl_region_all.shape[1], n_field, replace=False)
    return fl_region_all[:, idxs], e_fl_region2_all[:, idxs]


def dataNorm(
    arr: np.ndarray, e_arr: np.ndarray, sigma_max: float = 4.0
) -> tuple[np.ndarray, np.ndarray]:
    """Mask 'sigma_max' sigma outliers and normalize arrays.

    This function masks outliers based on a given sigma threshold,
    and normalizes the input arrays.

    :param arr: Array with the data.
    :type arr: np.ndarray
    :param e_arr: Array with the errors.
    :type e_arr: np.ndarray
    :param sigma_max: Sigma threshold for outlier masking.
    :type sigma_max: float

    :return: Normalized data and scaled errors.
    :rtype: tuple[np.ndarray, np.ndarray]
    """
    for tarr in (arr, e_arr):
        med = np.nanmedian(tarr, axis=1, keepdims=True)
        std = np.nanstd(tarr, axis=1, keepdims=True)
        msk = (tarr < med - sigma_max * std) | (tarr > med + sigma_max * std)
        tarr[msk] = np.nan

    dmin = np.nanmin(arr, axis=1)  # shape (d,)
    data_norm_T = arr.T - dmin  # (N, d)
    dmax = np.nanmax(data_norm_T, axis=0)  # (d,)
    data_norm_T /= dmax

    e_scaled2_T = (e_arr.T / dmax) ** 2
    return data_norm_T.T, e_scaled2_T.T


def _likelihood_matrix(
    cl_region_T: np.ndarray,
    e_cl_region2_T: np.ndarray,
    region: np.ndarray,
    e_region2: np.ndarray,
) -> np.ndarray:
    r"""Obtain the likelihood, for each star in the cluster region ('cl_reg_prep'),
    of being a member of the region passed ('region').

    This is basically the core of the 'tolstoy' likelihood with some added
    weights.

    .. math::
        L_i = \sum_{j=1}^{N_r}
                 \frac{1}{\sqrt{\prod_{k=1}^d \sigma_{ijk}^2}}\;\;
                    exp \left[-\frac{1}{2} \sum_{k=1}^d
                       \frac{(q_{ik}-q_{jk})^2}{\sigma_{ijk}^2} \right ]

    where:
        - i: cluster region star
        - j: field region star
        - k: data dimension
        - :math:`L_i`: likelihood for star i in the cluster region
        - :math:`N_r`: number of stars in field region
        - d: number of data dimensions
        - :math:`\sigma_{ijk}^2`: sum of squared uncertainties for stars i,j in dimension k
        - :math:`q_{ik}`: data for star i in dimension k
        - :math:`q_{jk}`: data for star j in dimension k

    :param cl_region_T: Array with the cluster region data.
    :type cl_region_T: np.ndarray
    :param e_cl_region2_T: Array with the cluster region errors.
    :type e_cl_region2_T: np.ndarray
    :param region: Array with the field region data.
    :type region: np.ndarray
    :param e_region2: Array with the field region errors.
    :type e_region2: np.ndarray

    :return: Array with the likelihoods.
    :rtype: np.ndarray
    """
    diff = cl_region_T - region.T[None, :]  # (N_cl, N_r, d)
    e_sum2 = e_cl_region2_T + e_region2.T[None, :]

    # Handle 'nan' values. THIS IS IMPORTANT
    diff[np.isnan(diff)] = 0.0
    e_sum2[np.isnan(e_sum2)] = 1.0

    # In-place square avoids a full-size temporary
    np.square(diff, out=diff)
    diff /= e_sum2
    Dsum = diff.sum(axis=-1)

    # Clip max values in place. THIS IS IMPORTANT
    np.clip(Dsum, None, 50.0, out=Dsum)
    np.log(e_sum2, out=e_sum2)
    Dsum += e_sum2.sum(axis=-1)

    return np.exp(-0.5 * Dsum)


def _likelihood_sum_chunked(
    cl_region_T: np.ndarray,
    e_cl_region2_T: np.ndarray,
    region: np.ndarray,
    e_region2: np.ndarray,
    chunk_size: int = 500,
) -> np.ndarray:
    """Summed likelihood per cluster star, chunked over field stars."""
    n_f = region.shape[1]
    if n_f <= chunk_size:
        L = _likelihood_matrix(cl_region_T, e_cl_region2_T, region, e_region2)
        return np.nansum(L, axis=-1)

    out = np.zeros(cl_region_T.shape[0])
    for start in range(0, n_f, chunk_size):
        end = min(start + chunk_size, n_f)
        chunk = region[:, start:end]
        e_chunk = e_region2[:, start:end]
        L_chunk = _likelihood_matrix(cl_region_T, e_cl_region2_T, chunk, e_chunk)
        out += np.nansum(L_chunk, axis=-1)
    return out

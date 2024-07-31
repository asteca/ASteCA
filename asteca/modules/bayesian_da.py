import numpy as np
import warnings


def bayesian_mp(frame_arr, e_frame_arr, center, radius, N_cluster, bayesda_runs):
    """
    Bayesian field decontamination algorithm. See Perren et al. (2015)

    A larger radius appears to translate to more estimated members using the
    "density" method (which) is expected, but a smaller number of stars with
    P>0.5 (and a cleaner final sequence)
    """
    cl_region, e_cl_region, fl_region_all, e_fl_region_all, cl_reg_idxs = get_regions(
        frame_arr, e_frame_arr, center, radius
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

    # Proper shape for the likelihood function
    cl_region_T, e_cl_region2_T = cl_region.T[:, None], e_cl_region2.T[:, None]

    # Initial null probabilities for all stars in the cluster region.
    prob_old_arr = np.zeros(N_cl_region)
    # Probabilities for all stars in the cluster region.
    sum_cl_probs = np.zeros(N_cl_region)

    N_break = 50
    for r in range(bayesda_runs):
        # Select stars from the cluster region according to their
        # associated probabilities so far.
        if N_cl_region > N_cluster:
            if r == 0:
                # Initial run
                p = np.random.choice(N_cl_region, N_cluster, replace=False)
            else:
                p = np.random.choice(
                    N_cl_region,
                    N_cluster,
                    replace=False,
                    p=sum_cl_probs / sum_cl_probs.sum(),
                )
        else:
            p = np.ones(N_cl_region, dtype=int)

        # Generate a random field region
        fl_region, e_fl_region2 = generate_field_region(
            fl_region_all, e_fl_region2_all, n_field
        )
        # Compare cluster region with this field region
        fl_lkl = likelihood(cl_region_T, e_cl_region2_T, fl_region, e_fl_region2)
        # Compare cluster region with the most probable cluster members (so far)
        cl_lkl = likelihood(
            cl_region_T, e_cl_region2_T, cl_region[:, p], e_cl_region2[:, p]
        )

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
        else:
            prob_old_arr = np.array(probs)

    if r < bayesda_runs:
        print(f"Convergence reached at {r+1} runs")
    else:
        print(f"Maximum number of runs reached: {bayesda_runs}")

    probs_final = np.zeros(frame_arr.shape[1])
    probs_final[cl_reg_idxs] = probs

    return probs_final


def get_regions(frame_arr, e_frame_arr2, center, radius):
    """Identify stars inside/outside the cluster region"""
    ra, dec = frame_arr[0], frame_arr[1]
    dist = np.sqrt((ra - center[0]) ** 2 + (dec - center[1]) ** 2)
    msk_in_rad = dist <= radius
    cl_reg_idxs = np.arange(len(ra))[msk_in_rad]

    # Cluster region inside radius
    cl_region = frame_arr[2:, msk_in_rad]
    e_cl_region2 = e_frame_arr2[:, msk_in_rad]

    # Field region outside radius
    fl_region_all = frame_arr[2:, ~msk_in_rad]
    e_fl_region2_all = e_frame_arr2[:, ~msk_in_rad]

    return cl_region, e_cl_region2, fl_region_all, e_fl_region2_all, cl_reg_idxs


def generate_field_region(fl_region_all, e_fl_region2_all, n_field):
    """Generate random field regions."""
    idxs = np.random.choice(fl_region_all.shape[1], n_field, replace=False)

    fl_region = fl_region_all[:, idxs]
    e_fl_region2 = e_fl_region2_all[:, idxs]

    return fl_region, e_fl_region2


def dataNorm(arr, e_arr, sigma_max=4.0):
    """
    Mask 'sigma_max' sigma outliers (particularly important when PMs are used),
    and normalize arrays.
    """
    for tarr in (arr, e_arr):
        for dim in tarr:
            med, std = np.nanmedian(dim), np.nanstd(dim)
            msk = np.logical_or(
                dim < med - sigma_max * std, dim > med + sigma_max * std
            )
            dim[msk] = np.nan

    # Minimum values for all arrays
    dmin = np.nanmin(arr, 1)
    data_norm = arr.T - dmin
    dmax = np.nanmax(data_norm, 0)
    data_norm /= dmax

    # Scale errors
    e_scaled = e_arr.T / dmax
    # Square all uncertainties here
    e_scaled2 = e_scaled**2

    return data_norm.T, e_scaled2.T


def likelihood(cl_region_T, e_cl_region2_T, region, e_region2):
    """
    Obtain the likelihood, for each star in the cluster region ('cl_reg_prep'),
    of being a member of the region passed ('region').

    This is basically the core of the 'tolstoy' likelihood with some added
    weights.

    L_i = \sum_{j=1}^{N_r}
             \frac{1}{\sqrt{\prod_{k=1}^d \sigma_{ijk}^2}}\;\;
                exp \left[-\frac{1}{2} \sum_{k=1}^d
                   \frac{(q_{ik}-q_{jk})^2}{\sigma_{ijk}^2} \right ]

    where
    i: cluster region star
    j: field region star
    k: data dimension
    L_i: likelihood for star i in the cluster region
    N_r: number of stars in field region
    d: number of data dimensions
    \sigma_{ijk}^2: sum of squared uncertainties for stars i,j in dimension k
    q_{ik}: data for star i in dimension k
    q_{jk}: data for star j in dimension k

    """
    arr_diff = cl_region_T - region.T[None, :]
    e_sum2 = e_cl_region2_T + e_region2.T[None, :]

    # Handle 'nan' values. THIS IS IMPORTANT
    arr_diff[np.isnan(arr_diff)] = 0.0
    e_sum2[np.isnan(e_sum2)] = 1.0

    Dsum = ((arr_diff) ** 2 / e_sum2).sum(axis=-1)
    # Clip max values in place. THIS IS IMPORTANT
    np.clip(Dsum, a_min=None, a_max=50.0, out=Dsum)

    sum_M_j = np.exp(-0.5 * Dsum) / np.sqrt(np.prod(e_sum2, axis=-1))
    Lkl = np.nansum(sum_M_j, axis=-1)

    return Lkl

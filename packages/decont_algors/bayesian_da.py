
import numpy as np
import warnings
from .. import update_progress


def main(
    colors, plx_col, pmx_col, pmy_col, rv_col, bayesda_runs, bayesda_dflag,
        cl_region, field_regions):
    '''
    Bayesian field decontamination algorithm.
    '''
    print('Applying Bayesian DA ({} runs)'.format(bayesda_runs))

    # cl_region = [[id, x, y, mags, e_mags, cols, e_cols, kine, ek], [], ...]
    # len(cl_region) = number of stars inside the cluster's radius.
    # len(cl_region[_][3]) = number of magnitudes defined.
    # len(field_regions) = number of field regions.
    # len(field_regions[i]) = number of stars inside field region 'i'.

    # Remove data dimensions.
    mags, cols, kinem, e_mags, e_cols, e_kinem = rmDimensions(
        cl_region, colors, plx_col, pmx_col, pmy_col, rv_col, bayesda_dflag)
    # Magnitudes and colors (and their errors) for all stars in the cluster
    # region, stored with the appropriate format.
    cl_reg_prep, w_cl = reg_data(
        len(cl_region), mags, cols, kinem, e_mags, e_cols, e_kinem)
    # Normalize data.
    cl_reg_prep, N_msk_cl = dataNorm(cl_reg_prep)

    # Likelihoods between all field regions and the cluster region.
    fl_likelihoods, N_msk_fr = [], 0
    for fl_region in field_regions:
        n_fl = len(fl_region)

        # Obtain likelihood, for each star in the cluster region, of
        # being a field star.
        mags, cols, kinem, e_mags, e_cols, e_kinem = rmDimensions(
            fl_region, colors, plx_col, pmx_col, pmy_col, rv_col,
            bayesda_dflag)
        fl_reg_prep, w_fl = reg_data(
            n_fl, mags, cols, kinem, e_mags, e_cols, e_kinem)
        fl_reg_prep, N_msk = dataNorm(fl_reg_prep)
        N_msk_fr += N_msk

        fl_likelihoods.append([n_fl, likelihood(
            fl_reg_prep, w_fl, cl_reg_prep, w_cl)])

    if N_msk_cl != 0 or N_msk_fr != 0:
        print("Masked data (outliers): N_cl={}, N_frs={}".format(
            N_msk_cl, N_msk_fr))

    # Create copy of the cluster region to be shuffled below.
    clust_reg_shuffle, w_cl_shuffle = cl_reg_prep[:], w_cl[:]

    # Initial null probabilities for all stars in the cluster region.
    prob_avrg_old = np.zeros(len(cl_region))
    # Probabilities for all stars in the cluster region.
    runs_fields_probs = np.zeros(len(cl_region))

    # Run 'bayesda_runs*fl_likelihoods' times.
    N_total = 0
    for run_num in range(bayesda_runs):
        # Iterate through all the 'field stars' regions that were populated.
        for n_fl, fl_lkl in fl_likelihoods:

            if n_fl < len(cl_region):
                # TODO DEPRECATED June 2019
                # # Randomly shuffle the stars within the cluster region.
                # p = np.random.permutation(len(clust_reg_shuffle))
                # clust_reg_shuffle, w_cl_shuffle = clust_reg_shuffle[p],\
                #     w_cl_shuffle[p]
                # # Remove n_fl random stars from the cluster region and
                # # obtain the likelihoods for each star in this "cleaned"
                # # cluster region.
                # cl_lkl = likelihood(
                #     bayesda_weights, clust_reg_shuffle[n_fl:],
                #     w_cl_shuffle[n_fl:], cl_reg_prep, w_cl)

                # Select stars from the cluster region according to their
                # associated probabilities.
                n_memb = len(clust_reg_shuffle) - n_fl
                if n_memb > 0:
                    # Identify first run.
                    if N_total > 0:
                        # Select stars according to their probabilities so far.
                        p = np.random.choice(
                            len(clust_reg_shuffle), n_memb, replace=False,
                            p=runs_fields_probs / runs_fields_probs.sum())
                    else:
                        p = np.random.choice(
                            len(clust_reg_shuffle), n_memb, replace=False)
                else:
                    p = np.arange(len(clust_reg_shuffle))
                clust_reg_shuffle_nmemb, w_cl_shuffle_nmemb =\
                    clust_reg_shuffle[p], w_cl_shuffle[p]
                # cluster region.
                cl_lkl = likelihood(
                    clust_reg_shuffle_nmemb, w_cl_shuffle_nmemb, cl_reg_prep,
                    w_cl)
            else:
                # If there are *more* field region stars than the total of
                # stars within the cluster region (highly contaminated
                # cluster), assign zero likelihood of being a true member to
                # all stars within the cluster region.
                cl_lkl = np.ones(len(cl_region)) * 1e-7

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                # Bayesian probability for each star within the cluster region.
                bayes_prob = 1. / (1. + (fl_lkl / cl_lkl))
            # Replace possible nan values with 0.
            bayes_prob[np.isnan(bayes_prob)] = 0.
            N_total += 1
            runs_fields_probs += bayes_prob

        # Check if probabilities converged. If so, break out.
        prob_avrg_old, break_flag = break_check(
            prob_avrg_old, runs_fields_probs, bayesda_runs, run_num, N_total)
        if break_flag:
            print("| MPs converged (run {})".format(run_num))
            break
        update_progress.updt(bayesda_runs, run_num + 1)

    # Average all Bayesian membership probabilities into a single value for
    # each star inside 'cl_region'.
    memb_probs_cl_region = runs_fields_probs / N_total

    return memb_probs_cl_region


def rmDimensions(
    region, colors, plx_col, pmx_col, pmy_col, rv_col,
        bayesda_dflag):
    """
    Remove data dimensions turned off by the user.
    """
    region_z = list(zip(*region))

    # Put each magnitude, color, and kinematic parameter into a separate list.
    # Remove dimensions according to input flags.
    mags, cols, kinem = [], [], []
    if bayesda_dflag[0] == 'y':
        mags = list(zip(*region_z[3]))
    for i, flag in enumerate(bayesda_dflag[1:len(colors) + 1]):
        if flag == 'y':
            cols.append(list(zip(*region_z[5]))[i])
    for i, k_d in enumerate((plx_col, pmx_col, pmy_col, rv_col)):
        if k_d is not False and bayesda_dflag[1 + len(colors) + i] == 'y':
            kinem.append(list(zip(*region_z[7]))[i])

    # Uncertainties.
    e_mags, e_cols, e_kinem = [], [], []
    if bayesda_dflag[0] == 'y':
        e_mags = list(zip(*region_z[4]))
    for i, flag in enumerate(bayesda_dflag[1:len(colors) + 1]):
        if flag == 'y':
            e_cols.append(list(zip(*region_z[6]))[i])
    for i, k_d in enumerate((plx_col, pmx_col, pmy_col, rv_col)):
        if k_d is not False and bayesda_dflag[1 + len(colors) + i] == 'y':
            e_kinem.append(list(zip(*region_z[8]))[i])

    # Remove kinematic dimensions where *all* the elements are 'nan'.
    e_kinem = [
        e_kinem[i] for i, _ in enumerate(kinem) if not np.isnan(_).all()]
    kinem = [_ for _ in kinem if not np.isnan(_).all()]

    return mags, cols, kinem, e_mags, e_cols, e_kinem


def reg_data(N_reg, mags, cols, kinem, e_mags, e_cols, e_kinem):
    """
    Generate list with photometric data in the correct format, and obtain the
    "dimensional" weights used by the likelihood.
    """

    # Combine photometry and uncertainties.
    data = np.array(mags + cols + kinem)
    # Uncertainties are squared in dataNorm()
    e_data = np.array(e_mags + e_cols + e_kinem)
    # Generate array with the appropriate format.
    data_err = np.stack((data, e_data)).T

    # Total number of information dimensions.
    d_T = len(mags) + len(cols) + len(kinem)
    d_info = np.zeros(N_reg)
    for m in mags:
        d_info += ~np.isnan(m)
    for c in cols:
        d_info += ~np.isnan(c)
    for k in kinem:
        d_info += ~np.isnan(k)
    # Final "dimensional information" weight. Equals '1.' if the star
    # contains valid data in all the defined information dimensions. Otherwise
    # it is a smaller float, down to zero when the star has no valid data.
    wi = d_info / d_T
    # wi = np.ones(len(region))

    return data_err, wi


def dataNorm(data_arr):
    """
    Mask 4 sigma outliers (particularly important when PMs are used), and
    normalize arrays.
    """

    # Mask outliers (np.nan).
    N_msk = 0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # Number of dimensions
        _, d, _ = data_arr.shape
        for i in range(d):
            # For both data and uncertainties
            for j in range(2):
                med, std = np.nanmedian(data_arr[:, i, j]),\
                    np.nanstd(data_arr[:, i, j])
                msk = np.logical_or(
                    data_arr[:, i, j] < med - 4. * std,
                    data_arr[:, i, j] > med + 4. * std)
                data_arr[:, i, j][msk] = np.nan
                N_msk += sum(msk)

    # Minimum values for all arrays
    dmin = np.nanmin(data_arr[:, :, 0], 0)
    data_norm = data_arr[:, :, 0] - dmin
    dmax = np.nanmax(data_norm, 0)
    data_norm /= dmax

    # Scale errors
    e_scaled = data_arr[:, :, 1] / dmax
    # Square errors
    e_scaled = np.square(e_scaled)
    # Avoid divide by zero in likelihood function.
    e_scaled[e_scaled == 0.] = np.nan

    # Combine into proper shape
    data_norm = np.array([data_norm.T, e_scaled.T]).T

    return data_norm, N_msk


def likelihood(region, w_j, cl_reg_prep, w_i):
    """
    Obtain the likelihood, for each star in the cluster region ('cl_reg_prep'),
    of being a member of the region passed ('region').

    This is basically the core of the 'tolstoy' likelihood with some added
    weights.

    L_i = w_i \sum_{j=1}^{N_r}
             \frac{w_j}{\sqrt{\prod_{k=1}^d \sigma_{ijk}^2}}\;\;
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
    w_i: data dimensions weight for star i
    w_j: data dimensions weight for star j

    """
    # Data difference (cluster_region - region), for all dimensions.
    data_dif = cl_reg_prep[:, None, :, 0] - region[None, :, :, 0]
    # Sum of squared errors, for all dimensions.
    sigma_sum = cl_reg_prep[:, None, :, 1] + region[None, :, :, 1]

    # Handle 'nan' values.
    data_dif[np.isnan(data_dif)] = 0.
    sigma_sum[np.isnan(sigma_sum)] = 1.

    # Sum for all dimensions.
    Dsum = (np.square(data_dif) / sigma_sum).sum(axis=-1)
    # This makes the code substantially faster.
    np.clip(Dsum, a_min=None, a_max=50., out=Dsum)

    # Product of summed squared sigmas.
    sigma_prod = np.prod(sigma_sum, axis=-1)

    # All elements inside summatory.
    sum_M_j = w_j * np.exp(-0.5 * Dsum) / np.sqrt(sigma_prod)

    # Sum for all stars in this 'region'.
    sum_M = w_i * np.sum(sum_M_j, axis=-1)
    # np.clip(sum_M, a_min=1e-7, a_max=None, out=sum_M)

    return sum_M


def break_check(prob_avrg_old, runs_fields_probs, runs, run_num, N_total):
    """
    Check if DA converged to MPs within a 0.1% tolerance, for all stars inside
    the cluster region.
    """
    # Average all probabilities.
    prob_avrg = runs_fields_probs / N_total

    # Set flag.
    break_flag = False

    # Check if probabilities changed less than 0.1% with respect to the
    # previous iteration.
    if np.allclose(prob_avrg_old, prob_avrg, 0.001):
        # Check that at least 10% of iterations have passed.
        if run_num >= max(1, int(0.1 * runs)):
            # Arrays are equal within tolerance and enough iterations have
            # passed. Break out.
            break_flag = True

    if break_flag is False:
        # Store new array in old one and proceed to new iteration.
        prob_avrg_old = prob_avrg

    return prob_avrg_old, break_flag

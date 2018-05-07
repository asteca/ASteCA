
import numpy as np
from .. import update_progress


def break_check(prob_avrg_old, runs_fields_probs, runs, run_num):
    """
    Check if DA converged to MPs within a 0.1% tolerance, for all stars inside
    the cluster region.
    """
    # Average all probabilities.
    prob_avrg = np.asarray(runs_fields_probs).mean(1).mean(0)

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


def likelihood(region, w_r, cl_reg_prep, w_c):
    """
    Obtain the likelihood, for each star in the cluster region, of being a
    member of the region passed.

    This is basically the core of the 'tolstoy' likelihood.

    See https://github.com/asteca/ASteCA/issues/352 for an explanation of the
    changes made since the old version.
    """
    # Photometric difference (cluster_region - region), for all dimensions.
    phot_dif = np.array(cl_reg_prep)[:, None, :, 0] -\
        np.array(region)[None, :, :, 0]
    # Sum of squared photometric errors, for all dimensions.
    sigma_sum = np.array(cl_reg_prep)[:, None, :, 1] +\
        np.array(region)[None, :, :, 1]

    phot_dif[np.isnan(phot_dif)] = 0.
    sigma_sum[np.isnan(sigma_sum)] = 1.

    # Sum for all photometric dimensions.
    Dsum = (np.square(phot_dif) / sigma_sum).sum(axis=-1)
    # Product of summed squared sigmas.
    sigma_prod = np.prod(sigma_sum, axis=-1)
    # All elements inside summatory.
    sum_M_j = w_r * np.exp(-0.5 * Dsum) / np.sqrt(sigma_prod)
    # Sum for all stars in this 'region'.
    sum_M = np.sum(sum_M_j, axis=-1)
    sum_M[sum_M <= 1e-7] = 1e-7
    sum_M = sum_M * w_c

    return sum_M


def reg_data(region):
    """
    Generate list with photometric data in the correct format, and obtain the
    "dimensional" weights used by the likelihood.
    """
    region_z = list(zip(*region))
    # Square errors.
    sq_em, sq_ec = np.square(region_z[4]), np.square(region_z[6])
    # Put each magnitude and color into a separate list.
    mags, cols = list(zip(*region_z[3])), list(zip(*region_z[5]))
    e_mags, e_cols = list(zip(*sq_em)), list(zip(*sq_ec))

    # Combine photometry and errors.
    phot = np.array(mags + cols)
    ephot = np.array(e_mags + e_cols)
    # Generate array with the appropriate format.
    phot_err = np.stack((phot, ephot)).T

    # Total number of information dimensions.
    d_T = len(mags) + len(cols)
    d_info = np.zeros(len(region))
    for m in mags:
        d_info += ~np.isnan(m)
    for c in cols:
        d_info += ~np.isnan(c)
    # Final "dimensional information" weight. Equals 1. if the star
    # contains valid data in ell the defined information dimensions. Otherwise
    # it is a smaller float, up to zero when the star has no valid data.
    wi = d_info / d_T

    return phot_err, wi


def main(cl_region, field_regions, bayesda_runs):
    '''
    Bayesian field decontamination algorithm.
    '''
    print('Applying Bayesian DA ({} runs).'.format(bayesda_runs))

    # cl_region = [[id, x, y, mags, e_mags, cols, e_cols], [], [], ...]
    # len(cl_region) = number of stars inside the cluster's radius.
    # len(cl_region[_][3]) = number of magnitudes defined.
    # len(field_regions) = number of field regions.
    # len(field_regions[i]) = number of stars inside field region 'i'.

    # This list holds one sub-list per run. Each of those sub-lists holds N
    # sub-sub-lists (where N = len(field_region)) with the membership
    # probabilities assigned to each star in the 'cluster_region'.
    #
    # runs_fields_probs = [A_1, A_2, ..., A_runs]
    # A_i = [B_1, B_2, ..., B_N] ; N = len(field_regions)
    # B_j = [p_1, p_2, ..., p_n] ; n = len(cl_region)
    # A,B --> arrays ; p --> floats (Bayesian probabilities)
    runs_fields_probs = []

    # Initial null probabilities for all stars in the cluster region.
    prob_avrg_old = np.zeros(len(cl_region))

    # Magnitudes and colors (and their errors) for all stars in the cluster
    # region, stored with the appropriate format.
    cl_reg_prep, w_cl = reg_data(cl_region)
    # Likelihoods between all field regions and the cluster region.
    fl_likelihoods = []
    for fl_region in field_regions:
        # Obtain likelihood, for each star in the cluster region, of
        # being a field star.
        fl_reg_prep, w_fl = reg_data(fl_region)
        # Number of stars in field region.
        n_fl = len(fl_region)
        fl_likelihoods.append([n_fl, likelihood(
            fl_reg_prep, w_fl, cl_reg_prep, w_cl)])

    # Create copy of the cluster region to be shuffled below.
    clust_reg_shuffle, w_cl_shuffle = cl_reg_prep[:], w_cl[:]

    # Run 'bayesda_runs' times.
    for run_num in range(bayesda_runs):
        field_reg_probs = []
        # Iterate through all the 'field stars' regions that were populated.
        for n_fl, fl_lkl in fl_likelihoods:

            if n_fl < len(cl_region):
                # Randomly shuffle the stars within the cluster region.
                p = np.random.permutation(len(clust_reg_shuffle))
                clust_reg_shuffle, w_cl_shuffle = clust_reg_shuffle[p],\
                    w_cl_shuffle[p]
                # Remove n_fl random stars from the cluster region and
                # obtain the likelihoods for each star in this "cleaned"
                # cluster region.
                cl_lkl = likelihood(
                    clust_reg_shuffle[n_fl:], w_cl_shuffle[n_fl:],
                    cl_reg_prep, w_cl)
            else:
                # If there are *more* field region stars than the total of
                # stars within the cluster region (highly contaminated
                # cluster), assign zero likelihood of being a true member to
                # all stars within the cluster region.
                cl_lkl = np.array([1e-7] * len(cl_region))

            # Bayesian probability for each star within the cluster region.
            bayes_prob = 1. / (1. + (fl_lkl / cl_lkl))
            # Now we have 'bayesda_runs' probabilities for each star in
            # 'cl_region' of being an true cluster member (MPs) stored in
            # 'field_reg_probs', for this field region.
            field_reg_probs.append(bayes_prob)

        # Append this list to the list that holds all the runs.
        runs_fields_probs.append(field_reg_probs)

        # Check if probabilities converged. If so, break out.
        prob_avrg_old, break_flag = break_check(
            prob_avrg_old, runs_fields_probs, bayesda_runs, run_num)
        if break_flag:
            print('| MPs converged (run {}).'.format(run_num))
            break
        update_progress.updt(bayesda_runs, run_num + 1)

    # Average all Bayesian membership probabilities into a single value for
    # each star inside 'cl_region'.
    memb_probs_cl_region = np.asarray(runs_fields_probs).mean(1).mean(0)

    return memb_probs_cl_region

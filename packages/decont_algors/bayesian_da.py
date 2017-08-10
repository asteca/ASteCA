
import numpy as np
from .. import update_progress


def break_check(prob_avrg_old, runs_fields_probs, runs, run_num):
    '''
    Check if DA converged to MPs within a 0.1% tolerance, for all stars inside
    the cluster region.
    '''

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


def likelihood(region, cl_phot):
    '''
    Obtain the likelihood, for each star in the cluster region, of being a
    member of the region passed.
    '''
    # Number of stars in this region.
    N = len(region)
    # Store cleaned cluster/field region and full cluster_region as arrays
    # skipping IDs otherwise the whole array is converted to strings.
    r_phot = [zip(*zip(*region)[1:][2]), zip(*zip(*region)[1:][3]),
              zip(*zip(*region)[1:][4]), zip(*zip(*region)[1:][5])]

    # TODO The second [0] index means we are using the first magnitude and
    # color defined. This should be generalized to N magnitudes and M colors
    # eventually.
    # Full cluster region. Square errors in color and magnitude.
    P = [cl_phot[0][0], np.square(cl_phot[1][0]), cl_phot[2][0],
         np.square(cl_phot[3][0])]
    # Cleaned cluster/field region.
    Q = [np.asarray(r_phot[0][0]), np.square(r_phot[1][0]),
         np.asarray(r_phot[2][0]), np.square(r_phot[3][0])]

    # Small value used to replace zeros.
    epsilon = 1e-10
    # For every star in the full cluster region.
    clust_stars_probs = []
    for star in zip(*P):
        # Squares sum of errors.
        e_col_2 = np.maximum(star[3] + Q[3], epsilon)
        e_mag_2 = np.maximum(star[1] + Q[1], epsilon)
        # star[2] & Q[2] = colors
        # star[0] & Q[0] = magnitude
        B = np.square(star[2] - Q[2]) / e_col_2
        C = np.square(star[0] - Q[0]) / e_mag_2
        synth_stars = np.exp(-0.5 * (B + C)) / np.sqrt(e_col_2 * e_mag_2)

        # The likelihood for this cluster star is the sum over all region
        # stars.
        likelihood = synth_stars.sum() / N

        # Use 1e-10 to avoid nan and inf values.
        clust_stars_probs.append(max(likelihood, epsilon))

    return clust_stars_probs


def main(cl_region, field_regions, bayesda_runs):
    '''
    Bayesian field decontamination algorithm.
    '''
    print('Applying decontamination algorithm.')

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
    prob_avrg_old = np.array([0.] * len(cl_region))

    # Extract magnitudes and colors (and their errors) for all stars
    # in the cluster region. The [1:] is meant to leave out the IDs when
    # zipping, otherwise all floats are converted to strings.
    # This is done here to save time.
    cl_phot = [zip(*zip(*cl_region)[1:][2]), zip(*zip(*cl_region)[1:][3]),
               zip(*zip(*cl_region)[1:][4]), zip(*zip(*cl_region)[1:][5])]

    # Run 'runs' times.
    for run_num in range(bayesda_runs):

        # This list will hold the probabilities for each field region.
        field_reg_probs = [[] for _ in field_regions]
        # Iterate through all the 'field stars' regions that were
        # populated.
        for indx, fl_region in enumerate(field_regions):

            # Obtain likelihood, for each star in the cluster region, of
            # being a field star.
            fl_lkl = likelihood(fl_region, cl_phot)
            # Store number of stars in field region.
            n_fl = len(fl_region)

            # Randomly shuffle the stars in the cluster region.
            clust_reg_shuffle = cl_region[:]
            np.random.shuffle(clust_reg_shuffle)
            # Remove n_fl random stars from the cluster region and pass
            # it to the function that obtains the likelihoods for each
            # star in the "cleaned" cluster region, ie: P(B)
            if n_fl < len(cl_region):
                clust_reg_clean = clust_reg_shuffle[n_fl:]
            else:
                # If field region has more stars than the cluster region,
                # don't remove any star. This should not happen though.
                clust_reg_clean = clust_reg_shuffle
            # Obtain likelihood, for each star in the cluster region, of
            # being a true cluster member.
            cl_lkl = likelihood(clust_reg_clean, cl_phot)

            # Obtain Bayesian probability for each star in cl_region.
            n_cl = len(clust_reg_clean)
            p_a, p_b = np.array(cl_lkl), np.array(fl_lkl)
            bayes_prob = 1. / (1. + (n_fl * p_b) / (n_cl * p_a))
            # Store probabilities obtained with this field region.
            field_reg_probs[indx] = bayes_prob

        # Now we have the probabilities of each star in 'cl_region' of
        # being an actual cluster member (membership) in 'field_reg_probs',
        # one sub-list per field region.

        # Append this list to the list that holds all the runs.
        runs_fields_probs.append(field_reg_probs)

        # Check if probabilities converged. If so, break out.
        prob_avrg_old, break_flag = break_check(
            prob_avrg_old, runs_fields_probs, bayesda_runs, run_num)
        if break_flag:
            print('| MPs converged in iteration {}.'.format(run_num))
            break
        update_progress.updt(bayesda_runs, run_num + 1)

    # Average all Bayesian membership probabilities into a single value for
    # each star inside 'cl_region'.
    memb_probs_cl_region = np.asarray(runs_fields_probs).mean(1).mean(0)

    return memb_probs_cl_region

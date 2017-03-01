
import numpy as np


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


def sort_members(memb_lst):
    '''
    Sort this list first by the membership probability from max
    value (1) to min (0) and then by its error values and magnitude value,
    (in that order) from min to max value.
    item[7] is the star's memb_prob and item[3] its magnitude.
    '''
    membership_prob_avrg_sort = sorted(memb_lst, key=lambda item: (-item[7],
                                       item[4], item[6], item[3]))

    return membership_prob_avrg_sort


def mpas(cl_reg_rad, runs_fields_probs):
    """
    Append averaged probabilities to each star inside the cluster radius and
    sort by their values.
    """

    # cl_reg_rad = [star1, star2, star3, ...]
    # starX = [id, x, y, mags, em, cols, ec]
    # mags = [mag1, mag2, mag3, ...]

    # Average all Bayesian membership probabilities into a single value for
    # each star inside 'cl_reg_rad'.
    clust_reg_prob_avrg = np.asarray(runs_fields_probs).mean(1).mean(0)

    # Create new list appending the membership probability to each star inside
    # the cluster radius.
    temp_prob_members = []
    for st_indx, star in enumerate(cl_reg_rad):
        temp_prob_members.append(star + [round(clust_reg_prob_avrg[st_indx],
                                               3)])

    # Stars inside the cluster's radius are now saved in the list
    # 'temp_prob_members' where each item contains the data for each
    # star:
    # starX = [id, x, y, mags, em, cols, ec, memb_prob].

    # Sort members list.
    membership_prob_avrg_sort = sort_members(temp_prob_members)

    return membership_prob_avrg_sort


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


def main(clp, npd, da_params, **kwargs):
    '''
    Bayesian field decontamination algorithm.
    '''

    flag_no_fl_regs, cl_region, field_regions = [
        clp[_] for _ in ['flag_no_fl_regs', 'cl_region', 'field_regions']]
    mode_da, run_n = da_params

    # Check if at least one equal-sized field region was obtained.
    if mode_da in ('auto', 'manual') and flag_no_fl_regs:
        print("  WARNING: no field regions found. Will not\n"
              "  apply decontamination algorithm.")
        mode_da = 'skip'

    flag_decont_skip = False
    # Run algorithm for any of these selections.
    if mode_da in ('auto', 'manual'):

        print('Applying decontamination algorithm.')

        # Set total number of runs.
        runs = 1000 if mode_da == 'auto' else run_n

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
        milestones = [25, 50, 75, 100]
        for run_num in range(runs):

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
            prob_avrg_old, break_flag = break_check(prob_avrg_old,
                                                    runs_fields_probs,
                                                    runs, run_num)
            if break_flag:
                print('  MPs converged in iteration {}.'.format(run_num))
                break

            percentage_complete = (100.0 * (run_num + 1) / runs)
            while len(milestones) > 0 and percentage_complete >= milestones[0]:
                print "  {}%".format(milestones[0])
                # Remove that milestone from the list.
                milestones = milestones[1:]

    elif mode_da == 'read':
        print('Reading membership probabilities from file.')
        # Read IDs from file.
        memb_file = npd['memb_file']
        data = np.genfromtxt(memb_file, dtype=str, unpack=True)
        id_list = data[0].tolist()
        # Read probabilities from file.
        data = np.genfromtxt(memb_file, dtype=float, unpack=True)
        memb_probs = data[7].tolist()

        probs = []
        # Assign probabilities read from file according to the star's IDs.
        # Those stars not present in the list are assigned a very low value.
        for indx, star in enumerate(cl_region):
            if star[0] in id_list:
                # Index of star in file.
                i = id_list.index(star[0])
                # Assign the probability stored in file for this star.
                probs.append(memb_probs[i])
            else:
                probs.append(0.01)

        # Store probabilities in list.
        runs_fields_probs = [[probs]]

    elif mode_da == 'skip':
        print('Assign equal probabilities to all stars in cluster region.')
        # Assign equal probabilities to all stars.
        runs_fields_probs = [[[1.] * len(cl_region)]]
        flag_decont_skip = True

    # Call function to average all probabilities.
    memb_prob_avrg_sort = mpas(cl_region, runs_fields_probs)

    clp['memb_prob_avrg_sort'], clp['flag_decont_skip'] =\
        memb_prob_avrg_sort, flag_decont_skip
    return clp

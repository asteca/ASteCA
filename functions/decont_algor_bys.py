"""
@author: gabriel
"""

import numpy as np
import get_in_params as g


def mpas(cl_reg_rad, runs_fields_probs):
    """
    Append averaged probabilities to each star inside the cluster radius and
    sort by their values.
    """

    # cl_reg_rad = [[id,x,y,mag,e_mag,color1,e_col1], [], [], ...]

    # Average all Bayesian membership probabilities into a single value for
    # each star inside 'cl_reg_rad'.
    clust_reg_prob_avrg = np.asarray(runs_fields_probs).mean(1).mean(0)

    # Create new list appending the membership probability to each star inside
    # the cluster radius.
    temp_prob_members = []
    for st_indx, star in enumerate(cl_reg_rad):
            temp_prob_members.append(star +
            [round(clust_reg_prob_avrg[st_indx], 3)])

    # Stars inside the cluster's radius are now saved in the list
    # 'temp_prob_members' where each item contains the data for each
    # star: [id,x,y,mag,e_mag,color1,e_col1,memb_prob].

    # Sort this list first by the membership probability from max
    # value (1) to min (0) and then by its error values and magnitude value,
    # (in that order) from min to max value.
    # item[7] is the star's memb_prob and item[3] its magnitude.
    membership_prob_avrg_sort = sorted(temp_prob_members,
                                       key=lambda item: (-item[7], item[4],
                                                         item[6], item[3]))

    return membership_prob_avrg_sort


def likelihood(cl_fl_arr, cl_full_arr):
    '''
    Obtain the likelihood/membership probability for each star in the
    region defined as that inside the cluster's radius.
    '''

    # Store cluster/field region and full cluster_region as arrays
    # skipping IDs otherwise the whole array is converted to strings.
    # Cluster/field region.
    #cl_fl_arr = np.array(zip(*zip(*region)[1:]), dtype=float)
    N = len(cl_fl_arr)  # Number of stars in this region.

    # Small value used to replace zeros.
    epsilon = 1e-10

    ## Split full cluster region array.
    #P = np.split(np.array(cl_full_arr), 7, axis=1)
    ## Square errors in color and magnitude.
    #P[4][0] = np.square(P[4][0])  # magnitude error
    #P[6][0] = np.square(P[6][0])  # color1 error
    #P[6][1] = np.square(P[6][1])  # color2 error
    #P = np.hstack(P)

    ## Split array.
    #Q = np.split(np.array(cl_fl_arr), 7, axis=1)
    ## Square photometric errors.
    #Q[4][0] = np.square(Q[4][0])  # magnitude error
    #Q[6][0] = np.square(Q[6][0])  # color1 error
    #Q[6][1] = np.square(Q[6][1])  # color2 error
    ##Q[5] = np.square(Q[5])
    ##Q[3] = np.square(Q[3])

    # For every star in the full cluster region.
    clust_stars_probs = []
    for star in cl_full_arr:
        # Squares sum of errors.
        e_mag_2 = np.square(star[4][0]) + np.square(cl_fl_arr[id_s][4][0])
        e_col_12 = np.square(star[6][0]) + np.square(cl_fl_arr[id_s][6][0])
        e_col_22 = np.square(star[6][1]) + np.square(cl_fl_arr[id_s][6][1])
        # magnitude
        B = np.square(star[3][0] - cl_fl_arr[id_s][3][0]) / e_mag_2
        # color1
        C = np.square(star[5][0] - cl_fl_arr[id_s][5][0]) / e_col_12
        # color2
        D = np.square(star[5][1] - cl_fl_arr[id_s][5][1]) / e_col_22
        synth_stars = np.exp(-0.5 * (B + C + D)) / \
        np.sqrt(e_mag_2 * e_col_12 * e_col_22)

        # The likelihood for this cluster star is the sum over all region
        # stars.
        likelihood = synth_stars.sum() / N

        # Use 1e-10 to avoid nan and inf values.
        clust_stars_probs.append(max(likelihood, epsilon))

    return clust_stars_probs


def bys_da(flag_area_stronger, cl_region, field_regions, memb_file):
    '''
    Bayesian field decontamination algorithm.
    '''

    mode_da, run_n = g.da_params

    # Check if at least one field region was obtained.
    if mode_da in {'auto', 'manual'} and flag_area_stronger:
        print ("  WARNING: no field regions found. Skipping\n"
        "  decontamination algorithm.")
        mode_da = 'skip'

    flag_decont_skip = False
    # Run algorithm for any of these selections.
    if mode_da in {'auto', 'manual'}:

        print 'Applying decontamination algorithm.'

        # Set total number of runs.
        runs = 1000 if mode_da == 'auto' else run_n
        print len(cl_region)

        # Full cluster region without IDs.
        #cl_full_arr = np.array(zip(*zip(*cl_reg_rad)[1:]), dtype=float)

        # cl_region = [[id, x, y, [mags], [e_mags], [cols], [e_cols]], [], ...]
        # len(cl_region) = number of stars inside the cluster's radius.
        # len(field_region[i]) = number of stars inside field region 'i'

        # This list holds one sub-list per run. Each of those sub-lists holds N
        # sub-sub-lists (where N = len(field_region)) with the membership
        # probabilities assigned to each star in the 'cluster_region'.
        #
        # runs_fields_probs = [A_1, A_2, ..., A_runs]
        # A_i = [B_1, B_2, ..., B_N] ; N = len(field_regions)
        # B_j = [p_1, p_2, ..., p_n] ; n = len(cl_reg_rad)
        # A,B --> arrays ; p --> floats (Bayesian probabilities)
        runs_fields_probs = []

        # Run 'runs' times.
        milestones = [25, 50, 75, 100]
        for run_num in range(runs):

            # This list will hold the probabilities for each field region.
            field_reg_probs = [[] for _ in field_regions]
            # Iterate through all the 'field stars' regions that were populated.
            for indx, fl_region in enumerate(field_regions):

                # Obtain likelihoods for each star in the clean cluster region
                # using this field region, ie: P(A)
                reg_decont_fl = likelihood(fl_region, cl_region)
                # Store number of stars in field region.
                n_fl = len(fl_region)

                # Randomly shuffle the stars in the cluster region.
                clust_reg_shuffle = np.random.permutation(cl_region)
                # Remove n_fl random stars from the cluster region and pass
                # it to the function that obtains the likelihoods for each
                # star in the "cleaned" cluster region, ie: P(B)
                if n_fl < len(cl_region):
                    clust_reg_clean = clust_reg_shuffle[n_fl:]
                else:
                    # If field region has more stars than the cluster region,
                    # don't remove any star. This should not happen though.
                    clust_reg_clean = clust_reg_shuffle
                n_cl = len(clust_reg_clean)
                reg_decont_cl = likelihood(clust_reg_clean, cl_region)

                # Obtain Bayesian probability for each star in cl_reg_rad.
                p_a, p_b = np.array(reg_decont_fl), np.array(reg_decont_cl)
                bayes_prob = 1. / (1. + (n_fl * p_a) / (n_cl * p_b))
                # Store probabilities obtained with this field region.
                field_reg_probs[indx] = bayes_prob

            # Now we have the probabilities of each star in 'cl_reg_rad' of
            # being an actual cluster member (membership) in 'field_reg_probs',
            # one sub-list per field region.

            # Append this list to the list that holds all the runs.
            runs_fields_probs.append(field_reg_probs)

            percentage_complete = (100.0 * (run_num + 1) / runs)
            while len(milestones) > 0 and percentage_complete >= milestones[0]:
                print "  {}% done".format(milestones[0])
                # Remove that milestone from the list.
                milestones = milestones[1:]

    elif mode_da == 'read':
        print 'Reading membership probabilities from file.'
        # Read IDs from file.
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
        print 'Assign equal probabilities to all stars in cluster region.'
        # Assign equal probabilities to all stars.
        runs_fields_probs = [[[1.] * len(cl_region)]]
        flag_decont_skip = True

    # Call function to average all probabilities.
    memb_prob_avrg_sort = mpas(cl_region, runs_fields_probs)

    return memb_prob_avrg_sort, flag_decont_skip
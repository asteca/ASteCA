
import bayesian_da
import fixed_da
import read_da


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


def mpas(cl_region, memb_probs_cl_region):
    """
    Append probabilities to each star inside the cluster radius and
    sort by their values.
    """

    # cl_region = [star1, star2, star3, ...]
    # starX = [id, x, y, mags, em, cols, ec]
    # mags = [mag1, mag2, mag3, ...]

    # Create new list appending the membership probability to each star inside
    # the cluster radius.
    temp_prob_members = []
    for st_indx, star in enumerate(cl_region):
        temp_prob_members.append(
            star + [round(memb_probs_cl_region[st_indx], 3)])

    # Stars inside the cluster's radius are now saved in the list
    # 'temp_prob_members' where each item contains the data for each star:
    # starX = [id, x, y, mags, em, cols, ec, memb_prob].

    # Sort members list.
    membership_prob_sort = sort_members(temp_prob_members)

    return membership_prob_sort


def main(clp, npd, da_algor, bayesda_runs, fixedda_port, readda_idcol,
         readda_mpcol, **kwargs):
    """
    Apply selected decontamination algorithm.
    """

    # Check if at least one equal-sized field region was obtained.
    if da_algor == 'bayes' and clp['flag_no_fl_regs']:
        print("  WARNING: no field regions found. Can not apply Bayesian DA.")
        da_algor = 'skip'

    flag_decont_skip = False
    if da_algor == 'skip':
        print('Assign equal probabilities to all stars in cluster region.')
        memb_probs_cl_region = [1.] * len(clp['cl_region'])
        flag_decont_skip = True

    elif da_algor == 'bayes':
        memb_probs_cl_region = bayesian_da.main(
            clp['cl_region'], clp['field_regions'], bayesda_runs)

    elif da_algor == 'fixed':
        memb_probs_cl_region = fixed_da.main(clp['cl_region'], fixedda_port)

    elif da_algor == 'read':
        memb_probs_cl_region = read_da.main(
            clp['cl_region'], npd['memb_file'], readda_idcol, readda_mpcol)

    # Average and sort probabilities.
    memb_prob_avrg_sort = mpas(clp['cl_region'], memb_probs_cl_region)

    clp['memb_prob_avrg_sort'], clp['flag_decont_skip'] =\
        memb_prob_avrg_sort, flag_decont_skip
    return clp

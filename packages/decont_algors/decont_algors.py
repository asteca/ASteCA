
# import numpy as np
from . import bayesian_da
from . import read_da


def main(
    clp, npd, colors, id_ids, plx_col, pmx_col, pmy_col, da_algor,
        bayesda_runs, bayesda_dflag, **kwargs):
    """
    Apply selected decontamination algorithm.
    """

    # Check if at least one equal-sized field region was obtained
    # (used by the Bayesian DA).
    if da_algor == 'y' and clp['flag_no_fl_regs']:
        print("  WARNING: no field regions found. Can not apply Bayesian DA")
        da_algor = 'n'

    flag_decont_skip = False
    if da_algor == 'n':
        print("Assign equal probabilities to all stars in cluster region")
        memb_probs_cl_region = [1.] * len(clp['cl_region'])
        flag_decont_skip = True

    elif da_algor == 'y':
        memb_probs_cl_region = bayesian_da.main(
            colors, plx_col, pmx_col, pmy_col, bayesda_runs,
            bayesda_dflag, clp['cl_region'], clp['field_regions'])

    else:
        memb_probs_cl_region = read_da.main(
            clp['cl_region'], npd['data_file'], id_ids, da_algor)

    # Store for plotting in 'C1' block and adding to the members output file.
    clp['memb_probs_cl_region'] = memb_probs_cl_region

    # DEPRECATED 26/03/22
    # memb_probs_cl_region = filterMPs(memb_probs_cl_region, clp['cl_region'])

    # Append MPs and sort by probabilities.
    memb_prob_avrg_sort = mpas(clp['cl_region'], memb_probs_cl_region)

    clp['memb_prob_avrg_sort'], clp['flag_decont_skip'] =\
        memb_prob_avrg_sort, flag_decont_skip
    return clp


# DEPRECATED 26/03/22
# def filterMPs(memb_probs_cl_region, cl_region):
#     """
#     Pass along MPs only for stars in the *complete* dataset.
#     """

#     ids_c, ids_i = list(zip(*cl_region_c))[0], list(zip(*cl_region_i))[0]
#     memb_probs_cl_region_c = np.zeros(len(ids_c))
#     dict_ids_c = {_: i for i, _ in enumerate(ids_c)}
#     for i, id_i in enumerate(ids_i):
#         try:
#             memb_probs_cl_region_c[dict_ids_c[id_i]] = memb_probs_cl_region[i]
#         except KeyError:
#             pass

#     return memb_probs_cl_region_c


def mpas(cl_region, memb_probs_cl_region):
    """
    Append probabilities to each star inside the cluster radius and
    sort by their values.
    """

    # cl_region = [star1, star2, star3, ...]
    # starX = [id, x, y, mags, em, cols, ec, kine, ek]
    # mags = [mag1, mag2, mag3, ...]

    # Create new list appending the membership probability to each star inside
    # the cluster radius.
    temp_prob_members = []
    for st_indx, star in enumerate(cl_region):
        temp_prob_members.append(
            star + [round(memb_probs_cl_region[st_indx], 3)])

    # Stars inside the cluster's radius are now saved in the list
    # 'temp_prob_members' where each item contains the data for each star:
    # starX = [id, x, y, mags, em, cols, ec, kine, ek, memb_prob].

    # Sort members list.
    membership_prob_sort = sort_members(temp_prob_members)

    return membership_prob_sort


def sort_members(memb_lst):
    '''
    Sort this list first by the membership probability from max
    value (1) to min (0) and then by its main magnitude.

    shape(memb_lst) = (N_stars, 10)
    '''
    membership_prob_avrg_sort = sorted(
        memb_lst, key=lambda item: (-item[9], item[3][0]))

    return membership_prob_avrg_sort

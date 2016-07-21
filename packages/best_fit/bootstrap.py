
import numpy as np
import random
from ..inp import input_params as g
import genetic_algorithm
import obs_clust_prepare


def resample_replacement(obs_clust):
    '''
    Resamples the observed cluster with replacement. Used by the bootstrap
    process.
    '''
    obs_cl = np.array([random.choice(obs_clust) for _ in obs_clust],
                      dtype=float)

    return obs_cl


def main(err_lst, memb_prob_avrg_sort, completeness, ip_list,
         st_dist_mass):
    '''
    Bootstrap process, runs the selected algorithm a number of times each
    time generating a new observed cluster representation through resampling
    with replacement.
    '''

    best_fit_algor, N_b = g.bf_params[1], g.bf_params[-1]

    print 'Begin bootstrap process (%d).' % N_b

    # List that holds the parameters values obtained by the bootstrap
    # process.
    params_boot = []

    milestones = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    # Begin bootstrap block (run a minimum of two times).
    for i in range(N_b):

        # Resample cluster with replacement.
        obs_cl_r = resample_replacement(memb_prob_avrg_sort)
        # Obtain prepared observed cluster according to the likelihood method
        # selected.
        obs_cl = obs_clust_prepare.main(obs_cl_r)

        # Algorithm selected.
        if best_fit_algor == 'genet':
            # Let the GA algor know this call comes from the bootstrap
            # process so it will not print percentages to screen.
            flag_print_perc = False
            params_boot.append(genetic_algorithm.main(
                flag_print_perc, err_lst, obs_cl, completeness, ip_list,
                st_dist_mass)[0])

        percentage_complete = (100.0 * (i + 1) / max(N_b, 2))
        while len(milestones) > 0 and percentage_complete >= milestones[0]:
            print "  {}% done".format(milestones[0])
            # Remove that milestone from the list.
            milestones = milestones[1:]

    # Calculate errors for each parameter.
    isoch_fit_errors = np.std(params_boot, 0)
    # Errors can not be smaller than the largest step in each parameter.
    par_vals = ip_list[1]
    for i, p_er in enumerate(isoch_fit_errors):
        # If any parameter has a single valued range, assign an error of -1.
        if len(par_vals[i]) > 1:
            # Find largest delta in this parameter used values.
            largest_delta = np.diff(par_vals[i]).max()
            # Store the maximum value.
            isoch_fit_errors[i] = max(largest_delta, p_er)
        else:
            isoch_fit_errors[i] = -1.

    return isoch_fit_errors

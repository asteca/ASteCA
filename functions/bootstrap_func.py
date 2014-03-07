# -*- coding: utf-8 -*-
"""
Created on Fri Mar 07 10:50:05 2014

@author: gabriel
"""

import numpy as np
import random
from genetic_algorithm import gen_algor as g_a


def resample_replacement(stars_list):
    '''
    Resamples the observed cluster with replacement. Used by the bootstrap
    process.
    '''
    obs_clust = [random.choice(stars_list) for _ in stars_list]

    return obs_clust


def bootstrap(err_lst, memb_prob_avrg_sort, completeness, ip_list, bf_params,
             sc_params, ga_params, ps_params):
    '''
    Bootstrap process, runs the selected algorithm a number of times each
    time generating a new observed cluster representation through resampling
    with replacement.
    '''

    bf_flag, best_fit_algor, N_b = bf_params
    sys_sel = ps_params[1]

    print 'Start bootstrap process (%d).' % N_b

    # List that holds the parameters values obtained by the bootstrap
    # process.
    params_boot = []

    milestones = [15, 30, 45, 60, 75, 90, 100]
    # Begin bootstrap block (run a minimum of two times).
    for i in range(max(N_b, 2)):

        # Resample cluster with replacement.
        obs_clust = resample_replacement(memb_prob_avrg_sort)
        # Let the GA algor know this call comes from the bootstrap process
        # so it will not print percentages to screen.
        flag_print_perc = False

        if best_fit_algor == 'genet':
            # Genetic algorithm.
            params_boot.append(g_a(flag_print_perc, err_lst, obs_clust,
            completeness, ip_list, sc_params, ga_params, sys_sel)[0])

        percentage_complete = (100.0 * (i + 1) / N_b)
        while len(milestones) > 0 and percentage_complete >= milestones[0]:
            print "  {}% done".format(milestones[0])
            # Remove that milestone from the list.
            milestones = milestones[1:]

    # Calculate errors for each parameter.
    isoch_fit_errors = np.std(params_boot, 0)
    # Errors can not be smaller than the steps in each parameter.
    for i in range(4):
        isoch_fit_errors[i] = max(ps_params[i + 3][2],
        isoch_fit_errors[i])

    return params_boot, isoch_fit_errors
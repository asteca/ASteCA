# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 16:10:05 2014

@author: gabriel
"""

from .._in import get_in_params as g
from obs_clust_prepare import prepare as prep
from genetic_algorithm import gen_algor as g_a
from brute_force_algorithm import brute_force as b_f
from bootstrap_func import bootstrap
from synth_cluster import synth_clust as s_c
from ..errors.error_round import round_sig_fig as rsf
from get_N_IMF import N_IMF as N_imf
from move_isochrone import move_isoch


def synth_cl_plot(ip_list, isoch_fit_params, err_lst, completeness,
    st_dist_mass):
    '''
    For plotting purposes.
    '''
    # Get list of stored isochrones and their parameters.
    isoch_list, param_values = ip_list[0], ip_list[1]
    # Read best fit values for all parameters.
    m, a, e, d, mass, binar_f = isoch_fit_params[0]
    # Find indexes for metallicity and age. If indexes are not found due
    # to some difference in the significant figures, use the indices
    # [0, 0] to prevent the code from halting.
    try:
        m_i, a_i = param_values[0].index(m), param_values[1].index(a)
    except:
        m_i, a_i = [0, 0]
    # Generate shifted best fit isochrone.
    shift_isoch = move_isoch(isoch_list[m_i][a_i][:2], e, d)
    # Generate best fit synthetic cluster.
    synth_clst = s_c(err_lst, completeness, st_dist_mass, isoch_list[m_i][a_i],
        [-1., -1., e, d, mass, binar_f])

    return shift_isoch, synth_clst


def params_errors(ip_list, err_lst, memb_prob_avrg_sort, completeness,
    st_dist_mass, isoch_fit_params):
    '''
    Obtain errors for the fitted parameters.
    '''

    best_fit_algor, N_b = g.bf_params[1], g.bf_params[-1]

    if  best_fit_algor == 'brute':
        # Assign errors as the steps in each parameter.
        isoch_fit_errors = [p_rs[2] for p_rs in ip_list[2]]

    elif best_fit_algor == 'genet':
        if N_b >= 2:
            # Call bootstrap function with resampling to get the uncertainty
            # in each parameter.
            isoch_fit_errors = bootstrap(err_lst, memb_prob_avrg_sort,
                completeness, ip_list, st_dist_mass)
        else:
            print 'Skipping bootstrap process.'
            # No error assignment.
            isoch_fit_errors = [-1.] * len(isoch_fit_params[0])

    # If any parameter has a single valued range, assign an error of -1.
    for i, par_vals in enumerate(ip_list[1]):
        if min(par_vals) == max(par_vals):
            isoch_fit_errors[i] = -1.

    return isoch_fit_errors


def best_fit(err_lst, memb_prob_avrg_sort, completeness, ip_list):
    '''
    Perform a best fitting process to find the cluster's parameters:
    E(B-V), distance (distance modulus), age and metallicity.
    '''

    bf_flag, best_fit_algor, lkl_method, bin_method = g.bf_params[:-1]

    # Check if algorithm should run.
    if bf_flag:

        print 'Searching for optimal parameters.'

        obs_clust = prep(memb_prob_avrg_sort)

        # Obtain mass distribution using the selected IMF. We run it once
        # because the array only depends on the IMF selected.
        st_dist_mass = N_imf()

        # Call algorithm to calculate the likelihoods for the set of
        # isochrones and return the best fitting parameters.
        if best_fit_algor == 'brute':

            print 'Using Brute Force algorithm ({}).'.format(lkl_method + '/' +
            bin_method if lkl_method == 'dolphin' else lkl_method)
            # Brute force algorithm.
            isoch_fit_params = b_f(err_lst, obs_clust, completeness, ip_list,
                st_dist_mass)

        elif best_fit_algor == 'genet':

            print 'Using Genetic Algorithm ({}).'.format(lkl_method + '/' +
            bin_method if lkl_method == 'dolphin' else lkl_method)
            # Genetic algorithm.
            # Let the GA algor know this call comes from the main function
            # so it will print percentages to screen.
            flag_print_perc = True
            isoch_fit_params = g_a(flag_print_perc, err_lst, obs_clust,
                completeness, ip_list, st_dist_mass)

        print ("Best fit params obtained.")

        # Assign errors for each parameter.
        isoch_fit_errors = params_errors(ip_list, err_lst, memb_prob_avrg_sort,
            completeness, st_dist_mass, isoch_fit_params)

        # Generate shifted isochrone and synthetic cluster for plotting.
        # Do this BEFORE rounding the parameter values.
        shift_isoch, synth_clst = synth_cl_plot(ip_list, isoch_fit_params,
            err_lst, completeness, st_dist_mass)

        # Round errors to 1 significant digit and round params values
        # to the corresponding number of significant digits given by
        # the errors.
        isoch_fit_params[0], isoch_fit_errors = rsf(isoch_fit_params[0],
            isoch_fit_errors)

    else:
        # Pass empty lists to make_plots.
        print 'Skipping parameters fitting process.'
        isoch_fit_params, isoch_fit_errors, shift_isoch, synth_clst = \
        [[-1., -1., -1., -1., -1., -1.]], [-1., -1., -1., -1., -1., -1.], [], []

    bf_return = [isoch_fit_params, isoch_fit_errors, shift_isoch, synth_clst]

    return bf_return
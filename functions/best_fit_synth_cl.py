# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 16:10:05 2014

@author: gabriel
"""

from genetic_algorithm import gen_algor as g_a
from brute_force_algorithm import brute_force as b_f
from bootstrap_func import bootstrap
from synth_cluster import synth_clust as s_c
from get_IMF_PDF import IMF_PDF as i_p
from move_isochrone import move_isoch


def best_fit(err_lst, memb_prob_avrg_sort, completeness, ip_list, bf_params,
             sc_params, ga_params, ps_params):
    '''
    Perform a best fitting process to find the cluster's parameters:
    E(B-V), distance (distance modulus), age and metallicity.
    '''

    bf_flag, best_fit_algor, N_b = bf_params
    sys_sel = ps_params[1]

    # Check if algorithm should run.
    if bf_flag:

        print 'Searching for optimal parameters.'

        # Obtain the selected IMF's PDF. We run it once because the array only
        # depends on the IMF selected.
        imf_pdf = i_p(sc_params[0])
        # Replace the name of the IMF by its PDF.
        sc_params[0] = imf_pdf

        # Call algorithm to calculate the likelihoods for the set of
        # isochrones and return the best fitting parameters.
        if best_fit_algor == 'brute':

            print 'Using Brute Force algorithm.'
            # Brute force algorithm.
            isoch_fit_params = b_f(err_lst, memb_prob_avrg_sort, completeness,
                                   ip_list, sc_params, ga_params, sys_sel)
            # Assign errors as the steps in each parameter.
            isoch_fit_errors = [ps_params[i + 3][2] for i in range(4)]
            print 'Best fit params obtained (%0.4f, %0.2f, %0.2f, %0.2f).' % \
            (isoch_fit_params[0][0], isoch_fit_params[0][1],
            isoch_fit_params[0][2], isoch_fit_params[0][3])

        elif best_fit_algor == 'genet':

            print 'Using Genetic Algorithm.'

            # The first pass is done with no resampling to calculate the final
            # values. After that we resample to get the uncertainty in each
            # parameter.

            # Genetic algorithm.
            # Let the GA algor know this call comes from the main function
            # so it will print percentages to screen.
            flag_print_perc = True
            isoch_fit_params = g_a(flag_print_perc, err_lst,
                memb_prob_avrg_sort, completeness, ip_list, sc_params,
                ga_params, sys_sel)
            print 'Best fit params obtained (%0.4f, %0.2f, %0.2f, %0.2f).' % \
            (isoch_fit_params[0][0], isoch_fit_params[0][1],
            isoch_fit_params[0][2], isoch_fit_params[0][3])

            # Call bootstrap function.
            params_boot, isoch_fit_errors = bootstrap(err_lst,
                memb_prob_avrg_sort, completeness, ip_list, bf_params,
                sc_params, ga_params, ps_params)

        # For plotting purposes.
        # Get list of stored isochrones and their parameters.
        isoch_list, isoch_ma = ip_list[0], ip_list[1]
        # Read best fit values for all parameters.
        m, a, e, d = isoch_fit_params[0]
        # Find indexes for metallixity and age.
        m_indx, a_indx = next(((i, j) for i, x in enumerate(isoch_ma) for j, y
        in enumerate(x) if y == [m, a]), None)
        # Generate shifted best fit isochrone.
        shift_isoch = move_isoch(sys_sel,
                                 isoch_list[m_indx][a_indx][:2], e, d)
        # Generate best fit synthetic cluster.
        synth_clst = s_c(err_lst, completeness, sc_params,
                         isoch_list[m_indx][a_indx], [-1., -1., e, d], sys_sel)
    else:
        # Pass empty lists to make_plots.
        print 'Skipping parameters fitting process.'
        isoch_fit_params, isoch_fit_errors, shift_isoch, synth_clst = \
        [[-1., -1., -1., -1.]], [-1., -1., -1., -1.], [], []

    bf_return = [isoch_fit_params, isoch_fit_errors, shift_isoch, synth_clst]

    return bf_return
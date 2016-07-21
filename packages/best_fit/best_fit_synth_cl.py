
import numpy as np
from ..inp import input_params as g
import obs_clust_prepare
import genetic_algorithm
import brute_force_algor
import bootstrap
import synth_cluster
from ..errors import error_round
import imf
import move_isochrone


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
    # to some difference in the significant figures, use the indexes
    # [0, 0] to prevent the code from halting.
    try:
        m_i, a_i = param_values[0].index(m), param_values[1].index(a)
    except:
        m_i, a_i = [0, 0]
    # Generate shifted best fit isochrone.
    shift_isoch = move_isochrone.main(isoch_list[m_i][a_i][:2], e, d)
    # Generate best fit synthetic cluster.
    synth_clst = synth_cluster.main(
        err_lst, completeness, st_dist_mass, isoch_list[m_i][a_i],
        [-1., -1., e, d, mass, binar_f])

    return shift_isoch, synth_clst


def params_errors(ip_list, err_lst, memb_prob_avrg_sort, completeness,
                  st_dist_mass, isoch_fit_params):
    '''
    Obtain errors for the fitted parameters.
    '''

    best_fit_algor, N_b = g.bf_params[1], g.bf_params[-1]

    if best_fit_algor == 'brute':
        isoch_fit_errors = []
        # Assign errors as the largest step in each parameter.
        par_vals = ip_list[1]
        for pv in par_vals:
            # If any parameter has a single valued range, assign an error
            # of -1.
            if len(pv) > 1:
                # Find largest delta in this parameter used values.
                largest_delta = np.diff(pv).max()
                # Store the maximum value.
                isoch_fit_errors.append(largest_delta)
            else:
                isoch_fit_errors.append(-1.)

    elif best_fit_algor == 'genet':
        if N_b >= 2:
            # Call bootstrap function with resampling to get the uncertainty
            # in each parameter.
            isoch_fit_errors = bootstrap.main(
                err_lst, memb_prob_avrg_sort, completeness, ip_list,
                st_dist_mass)
        else:
            print('Skipping bootstrap process.')
            # No error assignment.
            isoch_fit_errors = [-1.] * len(isoch_fit_params[0])

    return isoch_fit_errors


def main(err_lst, memb_prob_avrg_sort, completeness, ip_list):
    '''
    Perform a best fitting process to find the cluster's fundamental
    parameters.
    '''

    bf_flag, best_fit_algor, lkl_method, bin_method = g.bf_params[:-1]

    # Check if algorithm should run.
    if bf_flag:

        print('Searching for optimal parameters.')

        obs_clust = obs_clust_prepare.main(memb_prob_avrg_sort)
        # Store for plotting purposes.
        syn_b_edges = obs_clust[1]

        # Obtain mass distribution using the selected IMF. We run it once
        # because the array only depends on the IMF selected.
        st_dist_mass = imf.main()

        # Call algorithm to calculate the likelihoods for the set of
        # isochrones and return the best fitting parameters.
        if best_fit_algor == 'brute':

            print('Using Brute Force algorithm ({}).'.format(
                lkl_method + '; ' + bin_method if lkl_method == 'dolphin'
                else lkl_method))
            # Brute force algorithm.
            isoch_fit_params = brute_force_algor.main(
                err_lst, obs_clust, completeness, ip_list, st_dist_mass)

        elif best_fit_algor == 'genet':

            print('Using Genetic Algorithm ({}).'.format(
                lkl_method + '; ' + bin_method if lkl_method == 'dolphin'
                else lkl_method))
            # Genetic algorithm.
            # Let the GA algor know this call comes from the main function
            # so it will print percentages to screen.
            flag_print_perc = True
            isoch_fit_params = genetic_algorithm.main(
                flag_print_perc, err_lst, obs_clust, completeness, ip_list,
                st_dist_mass)

        print("Best fit parameters obtained.")

        # Assign errors for each parameter.
        isoch_fit_errors = params_errors(ip_list, err_lst, memb_prob_avrg_sort,
                                         completeness, st_dist_mass,
                                         isoch_fit_params)

        # Generate shifted isochrone and synthetic cluster for plotting.
        # Do this BEFORE rounding the parameter values.
        shift_isoch, synth_clst = synth_cl_plot(ip_list, isoch_fit_params,
                                                err_lst, completeness,
                                                st_dist_mass)

        if not synth_clst.any():
            print("  WARNING: best fit synthetic cluster found is empty.")

        # Round errors to 1 significant digit and round params values
        # to the corresponding number of significant digits given by
        # the errors.
        isoch_fit_params[0], isoch_fit_errors = error_round.round_sig_fig(
            isoch_fit_params[0], isoch_fit_errors)

    else:
        # Pass empty lists to make_plots.
        print('Skipping parameters fitting process.')
        isoch_fit_params, isoch_fit_errors, shift_isoch, synth_clst, \
            syn_b_edges = [[-1., -1., -1., -1., -1., -1.]], \
            [-1., -1., -1., -1., -1., -1.], [], [], []

    bf_return = [isoch_fit_params, isoch_fit_errors, shift_isoch, synth_clst,
                 syn_b_edges]

    return bf_return

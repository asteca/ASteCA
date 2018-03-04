
import numpy as np
import max_mag_cut
import obs_clust_prepare
import genetic_algorithm
import brute_force_algor
import bootstrap
from ..synth_clust import extin_coefs
from ..synth_clust import imf


def params_errors(
    lkl_method, e_max, err_lst, completeness, fundam_params, cl_max_mag,
        max_mag_syn, theor_tracks, R_V, ext_coefs, st_dist_mass, N_fc,
        err_rnd, N_pop, N_gen, fit_diff, cross_prob, cross_sel, mut_prob, N_el,
        N_ei, N_es, lkl_binning, lkl_weight, best_fit_algor, isoch_fit_params,
        N_b):
    '''
    Obtain errors for the fitted parameters.
    '''
    if best_fit_algor == 'brute':
        isoch_fit_errors = []
        # Assign errors as the largest step in each parameter.
        for pv in fundam_params:
            # If any parameter has a single valued range, assign an error
            # of -1.
            if len(pv) > 1:
                # Find largest delta in this parameter used values.
                largest_delta = np.diff(pv).max()
                # Store the maximum value.
                isoch_fit_errors.append(largest_delta)
            else:
                isoch_fit_errors.append(np.nan)

    elif best_fit_algor == 'genet':
        if N_b >= 2:
            # Call bootstrap function with resampling to get the uncertainty
            # in each parameter.
            isoch_fit_errors = bootstrap.main(
                lkl_method, e_max, err_lst, completeness, fundam_params,
                cl_max_mag, max_mag_syn, theor_tracks, R_V, ext_coefs,
                st_dist_mass, N_fc, err_rnd, N_pop, N_gen, fit_diff,
                cross_prob, cross_sel, mut_prob, N_el, N_ei, N_es, lkl_binning,
                lkl_weight, best_fit_algor, N_b)
        else:
            print('Skip bootstrap process.')
            # No error assignment.
            isoch_fit_errors = [np.nan] * len(isoch_fit_params[0])

    return isoch_fit_errors


def main(clp, bf_flag, best_fit_algor, lkl_method, lkl_binning, lkl_weight,
         N_bootstrap, max_mag, IMF_name, m_high, m_sample_flag, R_V, bin_mr,
         fundam_params, N_pop, N_gen, fit_diff, cross_prob, cross_sel,
         mut_prob, N_el, N_ei, N_es, cmd_systs, filters, colors,
         theor_tracks, **kwargs):
    '''
    Perform a best fitting process to find the cluster's fundamental
    parameters.
    '''
    # Check if algorithm should run.
    if bf_flag:
        err_lst, cl_reg_fit, completeness, e_max = clp['err_lst'],\
            clp['cl_reg_fit'], clp['completeness'], clp['err_max']

        # Remove stars beyond the maximum magnitude limit, if it was set.
        cl_max_mag, max_mag_syn = max_mag_cut.main(cl_reg_fit, max_mag)

        # Process observed cluster. This list contains data used by the
        # likelihoods, and for plotting.
        obs_clust = obs_clust_prepare.main(
            cl_max_mag, lkl_method, lkl_binning, lkl_weight)

        # DELETE
        # print(filters)
        # print(colors)
        # cl_histo_f = obs_clust[2]
        # N, B = len(cl_max_mag), len(cl_histo_f)
        # print('N:', N, 'B:', B)
        # B_p = np.count_nonzero(cl_histo_f)
        # print('B != 0:', B_p)
        # # Mighell likelihood testing
        # print('Best chi approx:', B_p / (1. + float(N) / B_p))
        # mig_chi = np.sum(np.clip(cl_histo_f, 0, 1) / (cl_histo_f + 1.))
        # print("Best chi:", mig_chi)
        # DELETE

        # Obtain extinction coefficients.
        # This parameter determines the total number of sub-arrays for each
        # isochrone stored.
        ext_shape = len(theor_tracks[0][0])
        ext_coefs = extin_coefs.main(cmd_systs, filters, colors, ext_shape)

        # Obtain mass distribution using the selected IMF. We run it once
        # because the array only depends on the IMF selected.
        st_dist_mass = imf.main(
            IMF_name, m_high, m_sample_flag, fundam_params[4])

        # Store the number of defined filters and colors.
        N_fc = [len(filters), len(colors)]

        # HARDCODED: generate 1e6 random floats to use in the synthetic cluster
        # error adding.
        err_rnd = np.random.normal(0, 1, 1000000)

        print('Searching for optimal parameters.')

        # Call algorithm to calculate the likelihoods for the set of
        # isochrones and return the best fitting parameters.
        if best_fit_algor == 'brute':

            print('Using Brute Force algorithm ({}).'.format(
                lkl_method + '; ' + lkl_binning if lkl_method == 'dolphin'
                else lkl_method))
            # Brute force algorithm.
            isoch_fit_params = brute_force_algor.main(
                lkl_method, e_max, err_lst, completeness, max_mag_syn,
                fundam_params, obs_clust, theor_tracks, R_V, ext_coefs,
                st_dist_mass, N_fc, err_rnd)

        elif best_fit_algor == 'genet':

            print('Using Genetic Algorithm ({}).'.format(
                lkl_method + '; ' + lkl_binning if lkl_method == 'dolphin'
                else lkl_method))
            # Genetic algorithm.
            # Let the GA algor know this call comes from the main function
            # so it will print percentages to screen.
            flag_print_perc = True
            isoch_fit_params = genetic_algorithm.main(
                lkl_method, e_max, err_lst, completeness, max_mag_syn,
                fundam_params, obs_clust, theor_tracks, R_V, ext_coefs,
                st_dist_mass, N_fc, err_rnd, N_pop, N_gen, fit_diff,
                cross_prob, cross_sel, mut_prob, N_el, N_ei, N_es,
                flag_print_perc)

        print("Best fit parameters obtained.")

        # Assign uncertainties for each parameter.
        isoch_fit_errors = params_errors(
            lkl_method, e_max, err_lst, completeness, fundam_params,
            cl_max_mag, max_mag_syn, theor_tracks, R_V, ext_coefs,
            st_dist_mass, N_fc, err_rnd, N_pop, N_gen, fit_diff, cross_prob,
            cross_sel, mut_prob, N_el, N_ei, N_es, lkl_binning, lkl_weight,
            best_fit_algor, isoch_fit_params, N_bootstrap)
    else:
        # Pass empty lists to make_plots.
        print('Skip parameters fitting process.')
        cl_max_mag, max_mag_syn, ext_coefs, st_dist_mass, N_fc, err_rnd,\
            isoch_fit_params, isoch_fit_errors = [], -1., [], {}, [], [],\
            [[np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]],\
            [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]

    clp['cl_max_mag'], clp['max_mag_syn'], clp['ext_coefs'],\
        clp['st_dist_mass'], clp['N_fc'], clp['err_rnd'],\
        clp['isoch_fit_params'], clp['isoch_fit_errors'] =\
        cl_max_mag, max_mag_syn, ext_coefs, st_dist_mass, N_fc, err_rnd,\
        isoch_fit_params, isoch_fit_errors
    return clp

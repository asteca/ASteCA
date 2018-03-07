
import numpy as np
import max_mag_cut
import emcee_algor  # TODO add to needed packages
import obs_clust_prepare
import genetic_algorithm
import brute_force_algor
import bootstrap
from ..synth_clust import extin_coefs
from ..synth_clust import imf


def params_errors(best_fit_algor, args):
    '''
    Obtain uncertainties for the fitted parameters.
    '''
    if best_fit_algor == 'brute':
        fundam_params = args
        isoch_fit_errors = []
        # Assign errors as the largest step in each parameter.
        for pv in fundam_params:
            # If any parameter has a single valued range, assign 'nan'.
            if len(pv) > 1:
                # Find largest delta in this parameter used values.
                largest_delta = np.diff(pv).max()
                # Store the maximum value.
                isoch_fit_errors.append(largest_delta)
            else:
                isoch_fit_errors.append(np.nan)

    elif best_fit_algor == 'genet':
        isoch_fit_errors = bootstrap.main(*args)

    elif best_fit_algor == 'emcee':
        varIdxs, emcee_trace = args
        isoch_fit_errors = []
        # TODO hard-coded for 6 parameters
        j = 0
        print("Best sol (median +- (84, 16) perc)")
        for i in range(6):
            if i in varIdxs:
                pm = np.percentile(emcee_trace[i - j], 50)  # Median
                #  16th and 84th percentiles (1 sigma)
                ph = np.percentile(emcee_trace[i - j], 84) - pm
                pl = pm - np.percentile(emcee_trace[i - j], 16)
                print("  {:.4f} +- ({:.4f}, {:.4f})".format(pm, ph, pl))
                # TODO fix this
                err = .5 * (ph + pl) if max(ph, pl) > 0. else np.nan
                isoch_fit_errors.append(err)
            else:
                isoch_fit_errors.append(np.nan)
                j += 1

    return isoch_fit_errors


def main(clp, bf_flag, best_fit_algor, lkl_method, lkl_binning, lkl_weight,
         N_bootstrap, max_mag, IMF_name, m_high, m_sample_flag, R_V,
         fundam_params, N_pop, N_gen, fit_diff, cross_prob, cross_sel,
         mut_prob, N_el, N_ei, N_es, cmd_systs, filters, colors, theor_tracks,
         nwalkers=0, nsteps=0., nburn=0., **kwargs):
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

        # HARDCODED: generate random floats to use in the synthetic cluster
        # completeness removal and error adding.
        cmpl_rnd = np.random.uniform(0., 1., 1000000)
        err_rnd = np.random.normal(0., 1., 1000000)

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
                st_dist_mass, N_fc, cmpl_rnd, err_rnd)
            # Assign uncertainties for each parameter.
            isoch_fit_errors = params_errors(best_fit_algor, fundam_params)

        elif best_fit_algor == 'genet':

            print('Using Genetic Algorithm ({}).'.format(
                lkl_method + '; ' + lkl_binning if lkl_method == 'dolphin'
                else lkl_method))
            # Genetic algorithm.
            # Let the GA know this call comes from the main function
            # so it will print percentages to screen.
            flag_print_perc = True
            isoch_fit_params = genetic_algorithm.main(
                lkl_method, e_max, err_lst, completeness, max_mag_syn,
                fundam_params, obs_clust, theor_tracks, R_V, ext_coefs,
                st_dist_mass, N_fc, cmpl_rnd, err_rnd, N_pop, N_gen, fit_diff,
                cross_prob, cross_sel, mut_prob, N_el, N_ei, N_es,
                flag_print_perc)
            # Assign uncertainties.
            isoch_fit_errors = params_errors(
                best_fit_algor,
                [lkl_method, e_max, err_lst, completeness, fundam_params,
                 cl_max_mag, max_mag_syn, theor_tracks, R_V, ext_coefs,
                 st_dist_mass, N_fc, cmpl_rnd, err_rnd, N_pop, N_gen, fit_diff,
                 cross_prob, cross_sel, mut_prob, N_el, N_ei, N_es,
                 lkl_binning, lkl_weight, N_bootstrap, False,
                 isoch_fit_params])

        elif best_fit_algor == 'emcee':
            print('Using emcee algorithm ({}).'.format(
                lkl_method + '; ' + lkl_binning if lkl_method == 'dolphin'
                else lkl_method))
            isoch_fit_params = emcee_algor.main(
                lkl_method, e_max, err_lst, completeness, max_mag_syn,
                fundam_params, obs_clust, theor_tracks, R_V, ext_coefs,
                st_dist_mass, N_fc, cmpl_rnd, err_rnd, nwalkers, nsteps, nburn)
            # Assign uncertainties.
            isoch_fit_errors = params_errors(
                best_fit_algor, isoch_fit_params[2:])

        print("Best fit parameters obtained.")

    else:
        # Pass empty lists to make_plots.
        print('Skip parameters fitting process.')
        cl_max_mag, max_mag_syn, ext_coefs, st_dist_mass, N_fc, cmpl_rnd,\
            err_rnd, isoch_fit_params, isoch_fit_errors = [], -1., [], {}, [],\
            [], [], [[np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]],\
            [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]

    clp['cl_max_mag'], clp['max_mag_syn'], clp['ext_coefs'],\
        clp['st_dist_mass'], clp['N_fc'], clp['cmpl_rnd'], clp['err_rnd'],\
        clp['isoch_fit_params'], clp['isoch_fit_errors'] =\
        cl_max_mag, max_mag_syn, ext_coefs, st_dist_mass, N_fc, cmpl_rnd,\
        err_rnd, isoch_fit_params, isoch_fit_errors
    return clp

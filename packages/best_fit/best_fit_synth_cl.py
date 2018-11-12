
from ..core_imp import np
from . import max_mag_cut, obs_clust_prepare, brute_force_algor,\
    genetic_algorithm, bootstrap, emcee_algor, abcpmc_algor, ptemcee_algor
from ..synth_clust import extin_coefs
from ..synth_clust import imf


def params_errors(best_fit_algor, args):
    '''
    Obtain uncertainties for the fitted parameters.
    '''
    mean_boot_sol = .0
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
        # TODO fix this with #64
        isoch_fit_errors, mean_boot_sol = bootstrap.main(*args)

    elif best_fit_algor in ['ptemcee', 'emcee', 'abc']:
        isoch_fit_params, isoch_fit_errors, j = args, [], 0
        for i, _ in enumerate(isoch_fit_params['mean_sol']):
            if i in isoch_fit_params['varIdxs']:
                #  16th and 84th percentiles (1 sigma)
                ph = np.percentile(isoch_fit_params['mcmc_trace'][i - j], 84)
                pl = np.percentile(isoch_fit_params['mcmc_trace'][i - j], 16)
                std = np.std(isoch_fit_params['mcmc_trace'][i - j])
                # Use maximum error value. TODO is this a sensible approach?
                err = max(.5 * (ph - pl), std)
                isoch_fit_errors.append(err)
            else:
                isoch_fit_errors.append(np.nan)
                j += 1

    return isoch_fit_errors, mean_boot_sol


def main(clp, bf_flag, best_fit_algor, lkl_method, lkl_binning,
         lkl_weight, N_bootstrap, max_mag, IMF_name, m_high, m_sample_flag,
         R_V, fundam_params, N_pop, N_gen, fit_diff, cross_prob, cross_sel,
         mut_prob, N_el, N_ei, N_es, cmd_systs, filters, colors, theor_tracks,
         nwalkers_emc, nsteps_emc, N_burn_emc, nburn_emc, emcee_a, priors_emc,
         ntemps, nwalkers_ptm, nsteps_ptm, nburn_ptm, pt_adapt, tmax_ptm,
         priors_ptm, nwalkers_abc, nsteps_abc, nburn_abc, priors_abc,
         **kwargs):
    '''
    Perform a best fitting process to find the cluster's fundamental
    parameters.
    '''
    # Check if algorithm should run.
    if bf_flag:
        err_lst, cl_reg_fit, completeness = clp['err_lst'],\
            clp['cl_reg_fit'], clp['completeness']

        # Remove stars beyond the maximum magnitude limit, if it was set.
        cl_max_mag, max_mag_syn = max_mag_cut.main(cl_reg_fit, max_mag)

        # Process observed cluster. This list contains data used by the
        # likelihoods, and for plotting.
        obs_clust = obs_clust_prepare.main(
            cl_max_mag, lkl_method, lkl_binning, lkl_weight)

        # TODO DELETE
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
                lkl_method, clp['em_float'], err_lst, completeness,
                max_mag_syn, fundam_params, obs_clust, theor_tracks, R_V,
                ext_coefs, st_dist_mass, N_fc, cmpl_rnd, err_rnd)
            # Assign uncertainties for each parameter.
            isoch_fit_errors, _ = params_errors(best_fit_algor, fundam_params)

        elif best_fit_algor == 'genet':

            print('Using Genetic Algorithm ({}).'.format(
                lkl_method + '; ' + lkl_binning if lkl_method == 'dolphin'
                else lkl_method))
            # Genetic algorithm.
            # Let the GA know this call comes from the main function
            # so it will print percentages to screen.
            flag_print_perc = True
            isoch_fit_params = genetic_algorithm.main(
                lkl_method, clp['em_float'], err_lst, completeness,
                max_mag_syn, fundam_params, obs_clust, theor_tracks, R_V,
                ext_coefs, st_dist_mass, N_fc, cmpl_rnd, err_rnd, N_pop, N_gen,
                fit_diff, cross_prob, cross_sel, mut_prob, N_el, N_ei, N_es,
                flag_print_perc)
            # Assign uncertainties.
            isoch_fit_errors, mean_boot_sol = params_errors(
                best_fit_algor,
                [lkl_method, clp['em_float'], err_lst, completeness,
                 fundam_params, cl_max_mag, max_mag_syn, theor_tracks, R_V,
                 ext_coefs, st_dist_mass, N_fc, cmpl_rnd, err_rnd, N_pop,
                 N_gen, fit_diff, cross_prob, cross_sel, mut_prob, N_el, N_ei,
                 N_es, lkl_binning, lkl_weight, N_bootstrap, False,
                 isoch_fit_params])
            # TODO fix this with #64
            isoch_fit_params['mean_sol'] = mean_boot_sol

        elif best_fit_algor == 'emcee':

            print('Using emcee algorithm ({}).'.format(
                lkl_method + '; ' + lkl_binning if lkl_method == 'dolphin'
                else lkl_method))
            isoch_fit_params = emcee_algor.main(
                lkl_method, clp['em_float'], err_lst, completeness,
                max_mag_syn, fundam_params, obs_clust, theor_tracks, R_V,
                ext_coefs, st_dist_mass, N_fc, cmpl_rnd, err_rnd, nwalkers_emc,
                nsteps_emc, nburn_emc, N_burn_emc, emcee_a, priors_emc)
            # Assign uncertainties.
            isoch_fit_errors, _ = params_errors(
                best_fit_algor, isoch_fit_params)

        elif best_fit_algor == 'abc':

            print('Using abcpmc algorithm ({}).'.format(
                lkl_method + '; ' + lkl_binning if lkl_method == 'dolphin'
                else lkl_method))
            isoch_fit_params = abcpmc_algor.main(
                lkl_method, clp['em_float'], err_lst, completeness,
                max_mag_syn, fundam_params, obs_clust, theor_tracks, R_V,
                ext_coefs, st_dist_mass, N_fc, cmpl_rnd, err_rnd, nwalkers_abc,
                nsteps_abc, nburn_abc, priors_abc)
            # Assign uncertainties.
            isoch_fit_errors, _ = params_errors(
                best_fit_algor, isoch_fit_params)

        elif best_fit_algor == 'ptemcee':
            print('Using ptemcee algorithm ({}).'.format(
                lkl_method + '; ' + lkl_binning if lkl_method == 'dolphin'
                else lkl_method))
            isoch_fit_params = ptemcee_algor.main(
                lkl_method, clp['em_float'], err_lst, completeness,
                max_mag_syn, fundam_params, obs_clust, theor_tracks, R_V,
                ext_coefs, st_dist_mass, N_fc, cmpl_rnd, err_rnd, ntemps,
                nwalkers_ptm, nsteps_ptm, nburn_ptm, pt_adapt, tmax_ptm,
                priors_ptm)
            # Assign uncertainties.
            isoch_fit_errors, _ = params_errors(
                best_fit_algor, isoch_fit_params)

        print("Best fit parameters obtained.")

    else:
        # Pass dummy data to make_plots.
        print('Skip parameters fitting process.')
        cl_max_mag, max_mag_syn, ext_coefs, st_dist_mass, N_fc, cmpl_rnd,\
            err_rnd, isoch_fit_params, isoch_fit_errors = [], -1., [], {}, [],\
            [], [],\
            {'mean_sol': [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
             'map_sol': [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]},\
            [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]

    clp['cl_max_mag'], clp['max_mag_syn'], clp['ext_coefs'],\
        clp['st_dist_mass'], clp['N_fc'], clp['cmpl_rnd'], clp['err_rnd'],\
        clp['isoch_fit_params'], clp['isoch_fit_errors'] =\
        cl_max_mag, max_mag_syn, ext_coefs, st_dist_mass, N_fc, cmpl_rnd,\
        err_rnd, isoch_fit_params, isoch_fit_errors
    return clp

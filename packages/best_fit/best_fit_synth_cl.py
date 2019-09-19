
from ..core_imp import np
from . import max_mag_cut, obs_clust_prepare, bootstrap, ptemcee_algor
# brute_force_algor
# emcee_algor, abcpmc_algor,
# TODO in place for #397: hopp_algor
from ..synth_clust import extin_coefs
from ..synth_clust import imf


def main(clp, pd):
    '''
    Perform a best fitting process to find the cluster's fundamental
    parameters.
    '''

    # Check if algorithm should run.
    if pd['bf_flag']:

        # Remove stars beyond the maximum magnitude limit, if it was set.
        cl_max_mag, max_mag_syn = max_mag_cut.main(
            clp['cl_reg_fit'], pd['max_mag'])

        # Processed observed cluster.
        obs_clust = obs_clust_prepare.main(
            cl_max_mag, pd['lkl_method'], pd['lkl_binning'], pd['lkl_weight'])

        # Obtain extinction coefficients.
        # This parameter determines the total number of sub-arrays for each
        # isochrone stored.
        ext_shape = len(pd['theor_tracks'][0][0])
        ext_coefs = extin_coefs.main(
            pd['cmd_systs'], pd['filters'], pd['colors'], ext_shape)

        # Obtain mass distribution using the selected IMF. We run it once
        # because the array only depends on the IMF selected.
        st_dist_mass = imf.main(
            pd['IMF_name'], pd['m_high'], pd['fundam_params'][4])

        # Store the number of defined filters and colors.
        N_fc = [len(pd['filters']), len(pd['colors'])]

        # HARDCODED: generate random floats to use in the synthetic cluster
        # completeness removal and error adding.
        cmpl_rnd = np.random.uniform(0., 1., 1000000)
        err_rnd = np.random.normal(0., 1., 1000000)

        print("Searching for optimal parameters")

        # Calculate the best fitting parameters.
        if pd['best_fit_algor'] == 'boot+GA':

            print("Using bootstrap + Genetic Algorithm ({})".format(
                pd['lkl_method'] + '; ' + pd['lkl_binning'] if
                pd['lkl_method'] == 'dolphin' else pd['lkl_method']))
            isoch_fit_params = bootstrap.main(
                pd, clp, cl_max_mag, max_mag_syn, obs_clust, ext_coefs,
                st_dist_mass, N_fc, cmpl_rnd, err_rnd)
            # Assign uncertainties.
            isoch_fit_errors = params_errors(pd, isoch_fit_params)

        elif pd['best_fit_algor'] == 'ptemcee':
            print("Using ptemcee algorithm ({})".format(
                pd['lkl_method'] + '; ' + pd['lkl_binning'] if
                pd['lkl_method'] == 'dolphin' else pd['lkl_method']))

            isoch_fit_params = ptemcee_algor.main(
                clp['err_lst'], clp['completeness'], clp['em_float'],
                max_mag_syn, obs_clust, ext_coefs, st_dist_mass, N_fc,
                cmpl_rnd, err_rnd, **pd)
            # Assign uncertainties.
            isoch_fit_errors = params_errors(pd, isoch_fit_params)

        # TODO DEPRECATED May 2019
        # if best_fit_algor == 'brute':
        #     print('Using Brute Force algorithm ({}).'.format(
        #         lkl_method + '; ' + lkl_binning if lkl_method == 'dolphin'
        #         else lkl_method))
        #     # Brute force algorithm.
        #     isoch_fit_params = brute_force_algor.main()

        # TODO not working yet
        # elif best_fit_algor == 'emcee':
        #     print('Using emcee algorithm ({}).'.format(
        #         lkl_method + '; ' + lkl_binning if lkl_method == 'dolphin'
        #         else lkl_method))
        #     isoch_fit_params = emcee_algor.main()
        #     # Assign uncertainties.
        #     isoch_fit_errors = params_errors(
        #         best_fit_algor, isoch_fit_params)

        # TODO not working yet
        # elif best_fit_algor == 'abc':
        #     print('Using abcpmc algorithm ({}).'.format(
        #         lkl_method + '; ' + lkl_binning if lkl_method == 'dolphin'
        #         else lkl_method))
        #     isoch_fit_params = abcpmc_algor.main()
        #     # Assign uncertainties.
        #     isoch_fit_errors, _ = params_errors(
        #         best_fit_algor, isoch_fit_params)

        print("Best fit parameters obtained")

        clp['max_mag_syn'], clp['ext_coefs'], clp['st_dist_mass'], \
            clp['N_fc'], clp['cmpl_rnd'], clp['err_rnd'], =\
            max_mag_syn, ext_coefs, st_dist_mass, N_fc, cmpl_rnd, err_rnd

    else:
        print("Skip parameters fitting process")
        # Pass dummy data used by data output and some plot functions.
        cl_max_mag = []
        isoch_fit_params = {
            'mean_sol': [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
            'median_sol': [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
            'map_sol': [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
            'mode_sol': [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
            'param_r2': [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
            'N_total': np.nan}
        isoch_fit_errors = [[np.nan, np.nan, np.nan]] * 6

    clp['cl_max_mag'], clp['isoch_fit_params'], clp['isoch_fit_errors'] = \
        cl_max_mag, isoch_fit_params, isoch_fit_errors

    return clp


def params_errors(pd, isoch_fit_params):
    '''
    Obtain uncertainties for the fitted parameters.
    '''
    # DEPRECATED May 2019
    # if best_fit_algor == 'brute':
    #     fundam_params = args
    #     isoch_fit_errors = []
    #     # Assign errors as the largest step in each parameter.
    #     for pv in fundam_params:
    #         # If any parameter has a single valued range, assign 'nan'.
    #         if len(pv) > 1:
    #             # Find largest delta in this parameter used values.
    #             largest_delta = np.diff(pv).max()
    #             # Store the maximum value.
    #             isoch_fit_errors.append(largest_delta)
    #         else:
    #             isoch_fit_errors.append(np.nan)

    def assignUncertns(varIdxs, trace):
        isoch_fit_errors, j = [], 0
        for i in range(6):
            if i in varIdxs:
                # 16th and 84th percentiles (1 sigma), and STDDEV
                ph = np.percentile(trace[i - j], 84)
                pl = np.percentile(trace[i - j], 16)
                std = np.std(trace[i - j])
                isoch_fit_errors.append((pl, ph, std))
            else:
                isoch_fit_errors.append((np.nan, np.nan, np.nan))
                j += 1

        return isoch_fit_errors

    if pd['best_fit_algor'] in ('boot+GA'):

        pb = isoch_fit_params['params_boot']
        if pb.any():
            isoch_fit_errors = assignUncertns(isoch_fit_params['varIdxs'], pb)
        else:
            # No error assignment.
            isoch_fit_errors = [[np.nan] * 3] * 6

    elif pd['best_fit_algor'] in ('ptemcee'):  # , 'emcee', 'abc'
        isoch_fit_errors = assignUncertns(
            isoch_fit_params['varIdxs'], isoch_fit_params['mcmc_trace'])

    return isoch_fit_errors

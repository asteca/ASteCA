
import numpy as np
from . import max_mag_cut, obs_clust_prepare, bootstrap, ptemcee_algor
from . import emcee_algor
# brute_force_algor , abcpmc_algor,
from ..synth_clust import add_errors, imf, extin_coefs, synth_clust_gen


def main(clp, pd):
    """
    Perform a best fitting process to find the cluster's fundamental
    parameters.
    """

    # Check if algorithm should run.
    if pd['bf_flag']:
        print("Searching for optimal parameters")

        cl_max_mag, max_mag_syn, obs_clust, ext_coefs, st_dist_mass, N_fc,\
            err_pars = dataPrep(pd, clp)

        lkl_bin = pd['lkl_method'] + '; ' + pd['lkl_binning'] if\
            pd['lkl_method'] in ('dolphin', 'tremmel') else pd['lkl_method']
        print("Using {} algorithm ({})".format(pd['best_fit_algor'], lkl_bin))

        # Calculate the best fitting parameters.
        if pd['best_fit_algor'] == 'ptemcee':
            isoch_fit_params = ptemcee_algor.main(
                clp['completeness'], max_mag_syn, obs_clust, ext_coefs,
                st_dist_mass, N_fc, err_pars, **pd)

        elif pd['best_fit_algor'] == 'boot+GA':
            isoch_fit_params = bootstrap.main(
                pd, clp, cl_max_mag, max_mag_syn, obs_clust, ext_coefs,
                st_dist_mass, N_fc, err_pars)

        elif pd['best_fit_algor'] == 'emcee':
            isoch_fit_params = emcee_algor.main(
                clp['err_lst'], clp['completeness'], clp['em_float'],
                max_mag_syn, obs_clust, ext_coefs, st_dist_mass, N_fc,
                err_pars, **pd)

        # TODO DEPRECATED May 2019
        # if best_fit_algor == 'brute':
        #     # Brute force algorithm.
        #     isoch_fit_params = brute_force_algor.main()
        # TODO not working yet
        # elif best_fit_algor == 'abc':
        #     isoch_fit_params = abcpmc_algor.main()

        # Assign uncertainties.
        isoch_fit_errors = params_errors(pd, isoch_fit_params)
        print("Best fit parameters obtained")

    else:
        print("Skip parameters fitting process")
        # Pass dummy data used by data output and some plot functions.
        cl_max_mag, max_mag_syn, ext_coefs, st_dist_mass, N_fc, err_pars =\
            [[]] * 6
        obs_clust = [[]]
        isoch_fit_params = {
            'mean_sol': [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
            'median_sol': [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
            'map_sol': [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
            'mode_sol': [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
            'param_r2': [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
            'N_total': np.nan}
        isoch_fit_errors = [[np.nan, np.nan, np.nan]] * 6

    # In place for #239
    if pd['best_fit_algor'] == 'synth_gen':
        cl_max_mag, max_mag_syn, obs_clust, ext_coefs, st_dist_mass, N_fc,\
            err_pars = dataPrep(pd, clp)
        clp = synth_clust_gen.main(
            clp, max_mag_syn, obs_clust,
            ext_coefs, st_dist_mass, N_fc, err_pars, **pd)

    #
    clp['cl_max_mag'], clp['max_mag_syn'], clp['ext_coefs'],\
        clp['st_dist_mass'], clp['N_fc'], clp['err_pars'],\
        clp['bf_bin_edges'], clp['isoch_fit_params'],\
        clp['isoch_fit_errors'] = cl_max_mag, max_mag_syn, ext_coefs,\
        st_dist_mass, N_fc, err_pars, obs_clust[0], isoch_fit_params,\
        isoch_fit_errors

    return clp


def dataPrep(pd, clp):
    """
    Obtain several parameters needed for the best fit process / synthetic
    clusters generation.
    """

    # Remove stars beyond the maximum magnitude limit, if it was set.
    cl_max_mag, max_mag_syn = max_mag_cut.main(
        clp['cl_reg_fit'], pd['max_mag'])

    # Processed observed cluster.
    obs_clust = obs_clust_prepare.main(
        cl_max_mag, pd['lkl_method'], pd['lkl_binning'],
        pd['lkl_manual_bins'])

    # Obtain extinction coefficients.
    # This parameter determines the total number of sub-arrays for each
    # isochrone stored.
    ext_shape = len(pd['theor_tracks'][0][0])
    ext_coefs = extin_coefs.main(
        pd['cmd_systs'], pd['filters'], pd['colors'], ext_shape)

    if pd['synth_rand_seed'] is not None:
        print("Random seed set for synthetic clusters: {}".format(
            pd['synth_rand_seed']))

    # Obtain mass distribution using the selected IMF. We run it once
    # because the array only depends on the IMF selected.
    st_dist_mass = imf.main(pd['IMF_name'], pd['fundam_params'][4])

    # Store the number of defined filters and colors.
    N_fc = [len(pd['filters']), len(pd['colors'])]

    err_rand = add_errors.randIdxs(pd['lkl_method'])
    err_pars = clp['err_lst'], clp['em_float'], err_rand

    # # TEMPORARY
    # # Use for saving the necessary input for the 'perf_test.py' file.
    # import pickle
    # with open('perf_test.pickle', 'wb') as f:
    #     pickle.dump([
    #         obs_clust, cl_max_mag, pd['fundam_params'], pd['theor_tracks'],
    #         pd['lkl_method'], pd['R_V'], clp['completeness'],
    #         max_mag_syn, st_dist_mass, ext_coefs, N_fc, pd['m_ini_idx'],
    #         pd['binar_flag'], err_pars], f)
    # print("finished")
    # import sys
    # sys.exit()
    # # TEMPORARY

    return cl_max_mag, max_mag_syn, obs_clust, ext_coefs, st_dist_mass, N_fc,\
        err_pars


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

    elif pd['best_fit_algor'] in ('ptemcee', 'emcee'):  # , , 'abc'
        isoch_fit_errors = assignUncertns(
            isoch_fit_params['varIdxs'], isoch_fit_params['mcmc_trace'])

    return isoch_fit_errors

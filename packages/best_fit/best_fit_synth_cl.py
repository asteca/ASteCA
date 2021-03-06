
import numpy as np
import pickle
from . import max_mag_cut, obs_clust_prepare, ptemcee_algor
from ..synth_clust import add_errors, imf, extin_coefs, synth_clust_gen
from .mcmc_convergence import convergenceVals
from .bf_common import r2Dist, modeKDE, fillParams  # thinChain


def main(npd, pd, clp):
    """
    Perform a best fitting process to find the cluster's fundamental
    parameters.
    """

    # General parameters
    if pd['best_fit_algor'] != 'n':
        cl_max_mag, max_mag_syn, obs_clust, ext_coefs, st_dist_mass, N_fc,\
            err_pars = dataPrep(pd, clp)
    else:
        # Pass dummy data used by data output and some plot functions.
        cl_max_mag, max_mag_syn, ext_coefs, st_dist_mass, N_fc, err_pars =\
            [[]] * 6
        obs_clust = [[]]

    if pd['best_fit_algor'] == 'ptemcee':
        print("Searching for optimal parameters")
        lkl_bin = pd['lkl_method'] + '; ' + pd['lkl_binning'] if\
            pd['lkl_method'] in ('dolphin', 'tremmel') else pd['lkl_method']
        print("Using {} algorithm ({})".format(pd['best_fit_algor'], lkl_bin))

        # Calculate the best fitting parameters.
        isoch_fit_params = ptemcee_algor.main(
            clp['completeness'], max_mag_syn, obs_clust, ext_coefs,
            st_dist_mass, N_fc, err_pars, **pd)

        # TODO DEPRECATED May 2020
        # elif pd['best_fit_algor'] == 'emcee':
        #     isoch_fit_params = emcee_algor.main(
        #         clp['completeness'], max_mag_syn, obs_clust, ext_coefs,
        #         st_dist_mass, N_fc, err_pars, **pd)

        # elif pd['best_fit_algor'] == 'boot+GA':
        #     isoch_fit_params = bootstrap.main(
        #         pd, clp, cl_max_mag, max_mag_syn, obs_clust, ext_coefs,
        #         st_dist_mass, N_fc, err_pars)

        # TODO DEPRECATED May 2019
        # if best_fit_algor == 'brute':
        #     # Brute force algorithm.
        #     isoch_fit_params = brute_force_algor.main()
        # TODO not working yet
        # elif best_fit_algor == 'abc':
        #     isoch_fit_params = abcpmc_algor.main()

        # Save MCMC trace (and some other variables) to file
        mcmcSamplesSave(
            isoch_fit_params, pd['full_trace_flag'], npd['mcmc_file_out'])

        # Obtain convergence parameters
        isoch_fit_params = convergenceParams(isoch_fit_params, **pd)

        # Assign uncertainties.
        isoch_fit_errors = params_errors(pd, isoch_fit_params)
        print("Best fit parameters obtained")

    elif pd['best_fit_algor'] == 'read':
        isoch_fit_params = mcmcSamplesRead(
            npd['clust_name'], npd['mcmc_file_out'])
        isoch_fit_params = convergenceParams(isoch_fit_params, **pd)
        isoch_fit_errors = params_errors(pd, isoch_fit_params)
        print("Best fit parameters read from file")

    # In place for #239
    elif pd['best_fit_algor'] == 'synth_gen':
        synth_clust_gen.main(
            npd, clp, max_mag_syn, obs_clust, ext_coefs, st_dist_mass, N_fc,
            err_pars, **pd)

    else:
        print("Skip parameters fitting process")
        # Pass dummy data used by data output and some plot functions.
        isoch_fit_params = {
            'mean_sol': [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
            'median_sol': [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
            'map_sol': [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
            'mode_sol': [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
            'param_r2': [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
            'N_total': np.nan}
        isoch_fit_errors = [[np.nan, np.nan, np.nan]] * 6

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


def convergenceParams(isoch_fit_params, fundam_params, nburn_mcee, **kwargs):
    """
    """

    # Store burn-in chain phase.
    bi_steps = int(nburn_mcee * isoch_fit_params['cold_chain'].shape[0])
    # chains_nruns.shape: (bi_steps, nchains, ndim)
    chains_nruns = isoch_fit_params['cold_chain'][:bi_steps]
    # pars_chains_bi.shape: (ndim, nchains, bi_steps)
    pars_chains_bi = chains_nruns.T

    # After burn-in
    chains_nruns = isoch_fit_params['cold_chain'][bi_steps:]
    pars_chains = chains_nruns.T

    # Convergence parameters.
    tau_autocorr, acorr_t, med_at_c, all_taus, geweke_z, acorr_function,\
        mcmc_ess = convergenceVals(
            'ptemcee', isoch_fit_params['ndim'], isoch_fit_params['varIdxs'],
            chains_nruns, bi_steps)

    # Re-shape trace for all parameters (flat chain).
    # Shape: (ndim, runs * nchains)
    mcmc_trace = chains_nruns.reshape(-1, isoch_fit_params['ndim']).T

    # # Thin chains
    # mcmc_trace = thinChain(mcmc_trace, acorr_t)

    param_r2 = r2Dist(fundam_params, isoch_fit_params['varIdxs'], mcmc_trace)
    mode_sol, pardist_kde = modeKDE(
        fundam_params, isoch_fit_params['varIdxs'], mcmc_trace)

    # Mean and median.
    mean_sol = np.mean(mcmc_trace, axis=1)
    median_sol = np.median(mcmc_trace, axis=1)

    # Fill the spaces of the parameters not fitted with their fixed values.
    mean_sol = fillParams(fundam_params, isoch_fit_params['varIdxs'], mean_sol)
    mode_sol = fillParams(fundam_params, isoch_fit_params['varIdxs'], mode_sol)
    median_sol = fillParams(
        fundam_params, isoch_fit_params['varIdxs'], median_sol)

    # Total number of values used to estimate the parameter's distributions.
    N_total = mcmc_trace.shape[-1]

    isoch_fit_params.update({
        'mean_sol': mean_sol, 'median_sol': median_sol, 'mode_sol': mode_sol,
        'pardist_kde': pardist_kde, 'param_r2': param_r2,
        'mcmc_trace': mcmc_trace, 'pars_chains_bi': pars_chains_bi,
        'pars_chains': pars_chains,
        'tau_autocorr': tau_autocorr, 'acorr_t': acorr_t, 'med_at_c': med_at_c,
        'all_taus': all_taus, 'acorr_function': acorr_function,
        # 'max_at_c': max_at_c, 'min_at_c': min_at_c,
        # 'minESS': minESS, 'mESS': mESS, 'mESS_epsilon': mESS_epsilon,
        'geweke_z': geweke_z, 'mcmc_ess': mcmc_ess, 'N_total': N_total,
    })

    return isoch_fit_params


def params_errors(pd, isoch_fit_params):
    """
    Obtain uncertainties for the fitted parameters.
    """

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

    isoch_fit_errors = assignUncertns(
        isoch_fit_params['varIdxs'], isoch_fit_params['mcmc_trace'])

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

    # DEPRECATED May 2020
    # if pd['best_fit_algor'] in ('boot+GA'):

    #     pb = isoch_fit_params['params_boot']
    #     if pb.any():
    #         isoch_fit_errors = assignUncertns(isoch_fit_params['varIdxs'], pb)
    #     else:
    #         # No error assignment.
    #         isoch_fit_errors = [[np.nan] * 3] * 6

    # elif pd['best_fit_algor'] in ('ptemcee', 'emcee'):  # , , 'abc'

    return isoch_fit_errors


def mcmcSamplesSave(isoch_fit_params, full_trace_flag, mcmc_file_out):
    """
    Create output data file with the MCMC sampler values for each parameter,
    for all chains/walkers.
    """
    if full_trace_flag:
        with open(mcmc_file_out, 'wb') as f:
            pickle.dump(isoch_fit_params, f)

        print("OUtput of MCMC sampler saved to file")


def mcmcSamplesRead(clust_name, mcmc_file_out):
    """
    Create output data file with the MCMC sampler values for each parameter,
    for all chains/walkers.
    """
    mcmc_file_out = mcmc_file_out.replace('output' + '/' + clust_name, 'input')
    with open(mcmc_file_out, 'rb') as f:
        isoch_fit_params = pickle.load(f)

    return isoch_fit_params

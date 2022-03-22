
import numpy as np
from ..inp import data_IO
from . import max_mag_cut, obs_clust_prepare, ptemcee_algor
from ..synth_clust import synth_clust_gen
from .mcmc_convergence import convergenceVals
from .bf_common import r2Dist, modeKDE, fillParams  # thinChain


def main(npd, pd, td, clp):
    """
    Perform a best fitting process to find the cluster's fundamental
    parameters.
    """

    # General parameters
    if pd['best_fit_algor'] != 'n':
        clp = dataPrep(pd, clp)

    # # Un-comment to fit [Fe/H] instead of Z
    # td['fundam_params'][0] = list(np.log10(np.array(
    #     td['fundam_params'][0]) / 0.0152))

    if pd['best_fit_algor'] == 'ptemcee':
        print("Searching for optimal parameters")
        lkl_bin = pd['lkl_method'] + '; ' + pd['lkl_binning'] if\
            pd['lkl_method'] in ('dolphin', 'tremmel') else pd['lkl_method']
        print("Using {} algorithm ({})".format(pd['best_fit_algor'], lkl_bin))

        # from . import kombine_algor
        # isoch_fit_params = kombine_algor.main(
        #     clp['completeness'], clp['err_lst'],
        #     clp['max_mag_syn'], clp['obs_clust'], pd['lkl_method'],
        #     pd['pt_ntemps'], pd['pt_adapt'], pd['pt_tmax'],
        #     pd['nsteps_mcee'], pd['nwalkers_mcee'], pd['mins_max'], **td)

        # from . import dynesty_algor
        # isoch_fit_params = dynesty_algor.main(
        #     clp['completeness'], clp['err_lst'],
        #     clp['max_mag_syn'], clp['obs_clust'], pd['lkl_method'],
        #     pd['pt_ntemps'], pd['pt_adapt'], pd['pt_tmax'],
        #     pd['nsteps_mcee'], pd['nwalkers_mcee'], pd['mins_max'], **td)

        # from . import dynesty_algor_2
        # isoch_fit_params = dynesty_algor_2.main(
        #     clp['completeness'], clp['err_lst'],
        #     clp['max_mag_syn'], clp['obs_clust'], pd['lkl_method'],
        #     pd['pt_ntemps'], pd['pt_adapt'], pd['pt_tmax'],
        #     pd['nsteps_mcee'], pd['nwalkers_mcee'], pd['mins_max'], **td)

        # from . import pypesto_algor
        # isoch_fit_params = pypesto_algor.main(
        #     clp['completeness'], clp['err_lst'],
        #     clp['max_mag_syn'], clp['obs_clust'], pd['lkl_method'],
        #     pd['pt_ntemps'], pd['pt_adapt'], pd['pt_tmax'],
        #     pd['nsteps_mcee'], pd['nwalkers_mcee'], pd['mins_max'], **td)

        # from . import pyabc_algor
        # isoch_fit_params = pyabc_algor.main(
        #     clp['completeness'], clp['err_lst'],
        #     clp['max_mag_syn'], clp['obs_clust'], pd['lkl_method'],
        #     pd['pt_ntemps'], pd['pt_adapt'], pd['pt_tmax'],
        #     pd['nsteps_mcee'], pd['nwalkers_mcee'], pd['mins_max'], **td)

        # Calculate the best fitting parameters.
        isoch_fit_params = ptemcee_algor.main(
            clp['completeness'], clp['err_lst'],
            clp['max_mag_syn'], clp['obs_clust'], pd['lkl_method'],
            pd['pt_ntemps'], pd['pt_adapt'], pd['pt_tmax'],
            pd['nsteps_mcee'], pd['nwalkers_mcee'], pd['mins_max'], **td)

        # Save MCMC trace (and some other variables) to file
        if pd['save_trace_flag']:
            data_IO.dataSave(isoch_fit_params, npd['mcmc_file_out'])
            print("Output of MCMC sampler saved to file")

        # Obtain convergence parameters
        isoch_fit_params = convergenceParams(
            isoch_fit_params, pd['nburn_mcee'], td['fundam_params'])

        # Assign uncertainties.
        isoch_fit_errors = params_errors(isoch_fit_params)
        print("Best fit parameters obtained")

    elif pd['best_fit_algor'] == 'read':
        isoch_fit_params = data_IO.dataRead(
            npd['clust_name'], npd['mcmc_file_out'])
        isoch_fit_params = convergenceParams(
            isoch_fit_params, pd['nburn_mcee'], td['fundam_params'])
        isoch_fit_errors = params_errors(isoch_fit_params)
        print("Best fit parameters read from file")

    # In place for #239. NOT WORKING AS OF FEB 2021, after #488 #503 #506
    elif pd['best_fit_algor'] == 'synth_gen':
        synth_clust_gen.main(npd, clp, pd)
        import sys
        sys.exit()

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

    clp['isoch_fit_params'], clp['isoch_fit_errors'] = isoch_fit_params,\
        isoch_fit_errors

    return clp


def dataPrep(pd, clp):
    """
    Obtain several parameters needed for the best fit process / synthetic
    clusters generation.
    """

    # Remove stars beyond the maximum magnitude limit, if it was set.
    clp['cl_max_mag'], clp['max_mag_syn'] = max_mag_cut.main(
        clp['cl_reg_fit'], pd['max_mag'])

    # Processed observed cluster.
    clp['obs_clust'] = obs_clust_prepare.main(
        clp['cl_max_mag'], pd['lkl_method'], pd['lkl_binning'],
        pd['lkl_manual_bins'])

    return clp


def convergenceParams(isoch_fit_params, nburn_mcee, fundam_params):
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


def params_errors(isoch_fit_params):
    """
    Obtain uncertainties for the fitted parameters.
    """

    def assignUncertns(varIdxs, trace):
        isoch_fit_errors, j = [], 0
        for i in range(len(isoch_fit_params)):
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

    return isoch_fit_errors

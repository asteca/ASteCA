
import numpy as np
from ..inp import data_IO
from . import ptemcee_algor
from ..synth_clust import synth_clust_gen
from .mcmc_convergence import convergenceVals
from .bf_common import varPars, modeKDE, fillParams


def main(npd, pd, td, clp):
    """
    Perform a best fitting process to find the cluster's fundamental
    parameters.
    """
    # # Un-comment to fit [Fe/H] instead of Z
    # td['fundam_params'][0] = list(np.log10(np.array(
    #     td['fundam_params'][0]) / 0.0152))

    if pd['best_fit_algor'] == 'ptemcee':
        print("Searching for optimal parameters")
        lkl_bin = pd['lkl_method'] + '; ' + pd['lkl_binning'] if\
            pd['lkl_method'] in ('dolphin', 'tremmel') else pd['lkl_method']
        print("Using {} algorithm ({})".format(pd['best_fit_algor'], lkl_bin))

        clp['varIdxs'], clp['ndim'], clp['ranges'] = varPars(
            td['fundam_params'])

        model_proper = np.array([_[0] for _ in td['fundam_params']])

        # Pack common args for the 'synth_cluster.py' function.
        clp['syntClustArgs'] = (
            pd['DR_dist'], pd['alpha'], model_proper,
            clp['varIdxs'], clp['completeness'], clp['err_lst'],
            clp['max_mag_syn'], clp['N_obs_stars'], td['fundam_params'],
            td['ext_coefs'], td['N_fc'], td['m_ini_idx'], td['st_dist_mass'],
            td['theor_tracks'], td['rand_norm_vals'], td['rand_unif_vals'])

        # import matplotlib.pyplot as plt
        # mags_cols_cl = [[], []]
        # for mag in list(zip(*list(zip(*clp['cl_syn_fit']))[1:][2])):
        #     mags_cols_cl[0].append(mag)
        # for col in list(zip(*list(zip(*clp['cl_syn_fit']))[1:][4])):
        #     mags_cols_cl[1].append(col)
        # from .bf_common import getSynthClust
        # from ..synth_clust import move_isochrone

        # def plot(model):
        #     synth_clust = getSynthClust(
        #         model, True, clp['syntClustArgs'])[0].T

        #     # Selected solution values for all the parameters.
        #     am = 9.6
        #     zm, e, dr, d, R_V = model
        #     # Values in grid
        #     zg = np.argmin(abs(np.array(td['fundam_params'][0]) - zm))
        #     ag = np.argmin(abs(np.array(td['fundam_params'][1]) - am))
        #     # Move isochrone
        #     isochrone = np.array(td['theor_tracks'][zg][ag])
        #     # dr = 0
        #     rand_v = np.ones(isochrone.shape[-1]) * .5
        #     shift_isoch = move_isochrone.main(
        #         isochrone, e, dr, d, R_V, td['ext_coefs'], td['N_fc'],
        #         pd['DR_dist'], rand_v, rand_v, td['m_ini_idx'])
        #     shift_isoch = shift_isoch[:sum(td['N_fc'])]

        #     plt.scatter(mags_cols_cl[1][0], mags_cols_cl[0][0], alpha=.5)
        #     plt.scatter(synth_clust[1], synth_clust[0], marker='x', alpha=.5)
        #     plt.plot(shift_isoch[1], shift_isoch[0], c='k')
        #     # plt.gca().invert_yaxis()
        #     plt.grid(which='both')
        #     plt.xlim(.5, 3)
        #     plt.ylim(20.5, 10.5)
        #     plt.show()

        # plot(np.array([0.01534, 0.7, 0.25, 12.1, 3.1]))
        # breakpoint()

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
            pd['lkl_method'], pd['pt_ntemps'], pd['pt_adapt'], pd['pt_tmax'],
            pd['nsteps_mcee'], pd['nwalkers_mcee'], pd['mins_max'],
            td['priors_mcee'], td['fundam_params'], clp['obs_clust'],
            clp['varIdxs'], clp['ndim'], clp['ranges'], clp['syntClustArgs'])

        # Save MCMC trace (and some other variables) to file
        if pd['save_trace_flag']:
            data_IO.dataSave(isoch_fit_params, npd['mcmc_file_out'])
            print("Output of MCMC sampler saved to file")

        # Obtain convergence parameters
        isoch_fit_params = convergenceParams(
            pd['nburn_mcee'], td['fundam_params'], clp['ndim'],
            clp['varIdxs'], isoch_fit_params)

        # Assign uncertainties.
        isoch_fit_errors = params_errors(clp['varIdxs'], isoch_fit_params)
        print("Best fit parameters obtained")

    elif pd['best_fit_algor'] == 'read':
        isoch_fit_params = data_IO.dataRead(
            npd['clust_name'], npd['mcmc_file_out'])
        isoch_fit_params = convergenceParams(
            pd['nburn_mcee'], td['fundam_params'], clp['ndim'],
            clp['varIdxs'], isoch_fit_params)
        isoch_fit_errors = params_errors(clp['varIdxs'], isoch_fit_params)
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


def convergenceParams(
        nburn_mcee, fundam_params, ndim, varIdxs, isoch_fit_params):
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
            'ptemcee', ndim, varIdxs, chains_nruns, bi_steps)

    # Re-shape trace for all parameters (flat chain).
    # Shape: (ndim, runs * nchains)
    mcmc_trace = chains_nruns.reshape(-1, ndim).T

    # # Thin chains
    # mcmc_trace = thinChain(mcmc_trace, acorr_t)

    # DEPRECATED 02/04/22
    # param_r2 = r2Dist(fundam_params, isoch_fit_params['varIdxs'], mcmc_trace)
    mode_sol, pardist_kde = modeKDE(fundam_params, varIdxs, mcmc_trace)

    # Mean and median.
    mean_sol = np.mean(mcmc_trace, axis=1)
    median_sol = np.median(mcmc_trace, axis=1)

    # Fill the spaces of the parameters not fitted with their fixed values.
    mean_sol = fillParams(fundam_params, varIdxs, mean_sol)
    mode_sol = fillParams(fundam_params, varIdxs, mode_sol)
    median_sol = fillParams(fundam_params, varIdxs, median_sol)

    # Total number of values used to estimate the parameter's distributions.
    N_total = mcmc_trace.shape[-1]

    isoch_fit_params.update({
        'mean_sol': mean_sol, 'median_sol': median_sol, 'mode_sol': mode_sol,
        'pardist_kde': pardist_kde,
        'mcmc_trace': mcmc_trace, 'pars_chains_bi': pars_chains_bi,
        'pars_chains': pars_chains,
        'tau_autocorr': tau_autocorr, 'acorr_t': acorr_t, 'med_at_c': med_at_c,
        'all_taus': all_taus, 'acorr_function': acorr_function,
        # 'max_at_c': max_at_c, 'min_at_c': min_at_c,
        # 'minESS': minESS, 'mESS': mESS, 'mESS_epsilon': mESS_epsilon,
        'geweke_z': geweke_z, 'mcmc_ess': mcmc_ess, 'N_total': N_total,
    })

    return isoch_fit_params


def params_errors(varIdxs, isoch_fit_params):
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

    isoch_fit_errors = assignUncertns(varIdxs, isoch_fit_params['mcmc_trace'])

    return isoch_fit_errors

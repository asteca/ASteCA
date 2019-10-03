
import numpy as np
import warnings
import time as t
from . import likelihood
# from .emcee3rc2 import autocorr
from .mcmc_convergence import convergenceVals
from .bf_common import initPop, varPars, synthClust, rangeCheck, fillParams,\
    r2Dist, modeKDE  # , thinChain
from .ptemcee import sampler


def main(
    err_lst, completeness, e_max, max_mag_syn, obs_clust, ext_coefs,
    st_dist_mass, N_fc, m_ini, cmpl_rnd, err_rnd, lkl_method, fundam_params,
    theor_tracks, R_V, ntemps, nsteps_pt, nwalkers_pt, nburn_pt, pt_adapt,
        tmax_pt, priors_pt, hmax, **kwargs):
    """
    """

    varIdxs, ndim, ranges = varPars(fundam_params)
    # Pack synthetic cluster arguments.
    synthcl_args = [
        theor_tracks, e_max, err_lst, completeness, max_mag_syn, st_dist_mass,
        R_V, ext_coefs, N_fc, m_ini, cmpl_rnd, err_rnd]

    if tmax_pt in ('n', 'none', 'None'):
        Tmax = None
    elif tmax_pt == 'inf':
        Tmax = np.inf
    else:
        Tmax = float(tmax_pt)
    if ntemps in ('n', 'none', 'None'):
        ntemps = None
    else:
        ntemps = int(float(ntemps))

    # Start timing.
    max_secs = hmax * 60. * 60.
    available_secs = max(30, max_secs)
    # elapsed, start_t = 0., t.time()

    # Define Parallel tempered sampler
    ptsampler = sampler.Sampler(
        nwalkers_pt, ndim, loglkl, logp,
        loglargs=[fundam_params, synthcl_args, lkl_method, obs_clust, ranges,
                  varIdxs, priors_pt], Tmax=Tmax, ntemps=ntemps)

    ntemps = ptsampler.ntemps
    # Initial population.
    pos0 = initPop(
        ranges, varIdxs, lkl_method, obs_clust, fundam_params, synthcl_args,
        ntemps, nwalkers_pt, 'random', None, None)

    # Track how the acceptance fractions, and temperature swaps acceptance
    # fractions.
    afs, tswaps = [], []
    # Store for Lkl values for plotting.
    prob_mean, map_lkl, map_sol_old = [], [], [[], -np.inf]

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        N_steps_store, runs = 50, 0

        elapsed, start = 0., t.time()
        milestones = list(range(10, 101, 10))
        for i, (pos, lnprob, lnlike) in enumerate(ptsampler.sample(
                pos0, iterations=nsteps_pt, adapt=pt_adapt)):

            # # Update position
            # pos0 = ptsampler._p0

            # # Check for dodgy inputs.
            # if np.any(np.isinf(pos0)):
            #     print(runs, "infs")  # <-- REMOVE
            #     pos0[np.isinf(pos0)] = 0.
            # if np.any(np.isnan(pos0)):
            #     print(runs, "nans")  # <-- REMOVE
            #     pos0[np.isnan(pos0)] = 0.
            # print("1", runs)

            # Only check convergence every 'N_steps_store' steps
            if (i + 1) % N_steps_store:
                continue
            runs += 1

            # Temperature swap acceptance fractions.
            tswaps.append(ptsampler.tswap_acceptance_fraction)
            # Mean acceptance fractions for all temperatures.
            afs.append(np.mean(ptsampler.acceptance_fraction, axis=1))

            # Compute the autocorrelation time using 'emcee'. tol=0 means that
            # we'll always get an estimate even if it isn't trustworthy.
            # tau0 = autocorr.integrated_time(
            #     ptsampler.chain[0, :, :i + 1, :].transpose(1, 0, 2), tol=0)
            # tau_emcee[tau_index] = np.mean(tau0)

            # # Autocorrelation time for the non-tempered chain. Mean across
            # # chains.
            # # ptsampler.chain.shape: (ntemps, nwalkers, nsteps, ndim)
            # # x.shape: (nsteps, ndim)
            # x = np.mean(ptsampler.chain[0, :, :i, :], axis=0)
            # # x = np.mean(ptsampler.chain[0], axis=0)
            # # tau.shape: ndim
            # tau = util.autocorr_integrated_time(x)
            # # Autocorrelation time. Mean across dimensions.
            # tau_autocorr.append(np.mean(tau))

            maf = np.mean(ptsampler.acceptance_fraction[0])
            # Store MAP solution in this iteration.
            prob_mean.append(np.mean(lnprob[0]))
            idx_best = np.argmax(lnprob[0])
            # Update if a new optimal solution was found.
            if lnprob[0][idx_best] > map_sol_old[1]:
                map_sol_old = [
                    fillParams(fundam_params, varIdxs, pos[0][idx_best]),
                    lnprob[0][idx_best]]
            map_lkl.append(map_sol_old[1])

            # Time used to check how fast the sampler is advancing.
            elapsed += t.time() - start
            start = t.time()
            # Print progress.
            percentage_complete = (100. * (i + 1) / nsteps_pt)
            if len(milestones) > 0 and percentage_complete >= milestones[0]:
                map_sol, logprob = map_sol_old
                m, s = divmod(nsteps_pt / (i / elapsed) - elapsed, 60)
                h, m = divmod(m, 60)
                print("{:>3}% ({:.3f}) LP={:.1f} ({:.5f}, {:.3f}, {:.3f}, "
                      "{:.2f}, {:.0f}, {:.2f})".format(
                          milestones[0], maf, logprob, *map_sol) +
                      " [{:.0f} m/s | {:.0f}h{:.0f}m]".format(
                          (ntemps * nwalkers_pt * i) / elapsed, h, m))
                milestones = milestones[1:]

            # elapsed += t.time() - start_t
            # Stop when available time is consumed.
            if elapsed >= available_secs:
                print("  Time consumed")
                break
            # start_t = t.time()

            # # Time used to check how fast the sampler is advancing.
            # percentage_complete = (100. * elapsed / available_secs)
            # if len(milestones) > 0 and percentage_complete >= milestones[0]:
            #     map_sol, logprob = map_sol_old
            #     m, s = divmod(max(1., available_secs - elapsed), 60)
            #     h, m = divmod(m, 60)
            #     print("{:>3}% ({:.3f}) LP={:.1f} ({:.5f}, {:.3f}, {:.3f}, "
            #           "{:.2f}, {:.0f}, {:.2f})".format(
            #               milestones[0], maf, logprob, *map_sol) +
            #           " [{:.0f} m/s | {:.0f}h{:.0f}m]".format(
            #               (ntemps * nwalkers_ptm * runs * steps) /
            #               elapsed, h, m))
            #     milestones = milestones[1:]

    # Total number of steps
    N_steps = N_steps_store * np.arange(1, runs + 1)

    # Mean acceptance fractions for all replicas.
    maf_allT = np.asarray(afs).T
    # Temperature swaps acceptance fractions.
    tswaps_afs = np.asarray(tswaps).T
    # Betas history
    betas_pt = ptsampler.beta_history[:, (N_steps_store - 1)::N_steps_store]

    # Final MAP fit.
    map_sol, map_lkl_final = map_sol_old

    # This number should be between approximately 0.25 and 0.5 if everything
    # went as planned.
    m_accpt_fr = np.mean(ptsampler.acceptance_fraction[0])
    if m_accpt_fr > .5 or m_accpt_fr < .25:
        print("  WARNING: mean acceptance fraction is outside of the\n"
              "  recommended range.")

    # ptsampler.chain.shape: (ntemps, nchains, nsteps, ndim)
    # cold_chain.shape: (i, nchains, ndim)
    cold_chain = ptsampler.chain[0, :, :i, :].transpose(1, 0, 2)

    # Store burn-in chain phase.
    bi_steps = int(nburn_pt * cold_chain.shape[0])
    # chains_nruns.shape: (bi_steps, nchains, ndim)
    chains_nruns = cold_chain[:bi_steps]
    # pars_chains_bi.shape: (ndim, nchains, bi_steps)
    pars_chains_bi = chains_nruns.T

    # After burn-in
    chains_nruns = cold_chain[bi_steps:]
    pars_chains = chains_nruns.T

    # Convergence parameters.
    tau_autocorr, acorr_t, med_at_c, all_taus, geweke_z, acorr_function,\
        mcmc_ess = convergenceVals(
            'ptemcee', ndim, varIdxs, None, chains_nruns, bi_steps)

    # Re-shape trace for all parameters (flat chain).
    # Shape: (ndim, runs * nchains)
    mcmc_trace = chains_nruns.reshape(-1, ndim).T

    # # Thin chains
    # mcmc_trace = thinChain(mcmc_trace, acorr_t)

    param_r2 = r2Dist(fundam_params, varIdxs, mcmc_trace)
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

    isoch_fit_params = {
        'varIdxs': varIdxs, 'mean_sol': mean_sol, 'Tmax': str(Tmax),
        'median_sol': median_sol, 'map_sol': map_sol, 'map_lkl': map_lkl,
        'mode_sol': mode_sol, 'pardist_kde': pardist_kde, 'param_r2': param_r2,
        'map_lkl_final': map_lkl_final, 'prob_mean': prob_mean,
        'bf_elapsed': elapsed, 'mcmc_trace': mcmc_trace,
        'pars_chains_bi': pars_chains_bi, 'pars_chains': pars_chains,
        'maf_allT': maf_allT, 'tswaps_afs': tswaps_afs, 'betas_pt': betas_pt,
        #
        'tau_autocorr': tau_autocorr, 'acorr_t': acorr_t, 'med_at_c': med_at_c,
        'all_taus': all_taus, 'acorr_function': acorr_function,
        # 'max_at_c': max_at_c, 'min_at_c': min_at_c,
        # 'minESS': minESS, 'mESS': mESS, 'mESS_epsilon': mESS_epsilon,
        'geweke_z': geweke_z, 'mcmc_ess': mcmc_ess, 'N_total': N_total,
        'N_steps': N_steps
    }

    return isoch_fit_params


def loglkl(
    model, fundam_params, synthcl_args, lkl_method, obs_clust, ranges,
        varIdxs, priors_ptm):
    """
    """
    rangeFlag = rangeCheck(model, ranges, varIdxs)

    logpost = -1e9
    if rangeFlag:
        # Generate synthetic cluster.
        synth_clust = synthClust(fundam_params, varIdxs, synthcl_args, model)
        # Call likelihood function for this model.
        lkl = likelihood.main(lkl_method, synth_clust, obs_clust)
        log_p = 0.
        for i, pr in enumerate(priors_ptm):
            # Gaussian prior. If the prior is uniform, simply pass.
            if pr[0] == 'g':
                log_p += np.log(1 / pr[2]) - .5 * np.square(
                    (model[i] - pr[1]) / pr[2])

        # The negative likelihood is returned since Dolphin requires a
        # minimization of the PLR. Here we are maximizing, hence the minus.
        logpost = log_p + (-lkl)

    return logpost


def logp(_):
    """
    Just here as a place holder for 'ptemcee'.

    The prior is moved inside the log-likelihood to save computation time
    because 'ptemcee' calls both functions and then adds the result. But if
    the model is outside of the allowed range, there is no point in checking
    twice (ie: calling rangeCheck() twice). As I can not "send" information
    from here to the likelihood (as can be done in 'emcee'), it's better to
    just put everything inside the loglkl() function.
    """
    return 0.

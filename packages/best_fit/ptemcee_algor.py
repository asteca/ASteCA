
import numpy as np
import warnings
import time as t
from .. import update_progress
from . import likelihood
# from .emcee3rc2 import autocorr
from .mcmc_convergence import convergenceVals
from .bf_common import initPop, varPars, synthClust, rangeCheck, fillParams,\
    closeSol, discreteParams, r2Dist, modeKDE  # , thinChain
from .ptemcee import sampler, util


def main(
    err_lst, completeness, e_max, max_mag_syn, obs_clust, ext_coefs,
    st_dist_mass, N_fc, cmpl_rnd, err_rnd, lkl_method, fundam_params,
    theor_tracks, R_V, init_mode_ptm, popsize_ptm, maxiter_ptm, ntemps,
    nwalkers_ptm, nsteps_ptm, nburn_ptm, pt_adapt, tmax_ptm, priors_ptm, hmax,
        N_conv, tol_conv, **kwargs):
    """
    """

    varIdxs, ndim, ranges = varPars(fundam_params)
    # Pack synthetic cluster arguments.
    synthcl_args = [
        theor_tracks, e_max, err_lst, completeness, max_mag_syn, st_dist_mass,
        R_V, ext_coefs, N_fc, cmpl_rnd, err_rnd]

    if tmax_ptm in ('n', 'none', 'None'):
        Tmax = None
    elif tmax_ptm == 'inf':
        Tmax = np.inf
    else:
        Tmax = float(tmax_ptm)
    if ntemps in ('n', 'none', 'None'):
        ntemps = None
    else:
        ntemps = int(float(ntemps))

    # Start timing.
    max_secs = hmax * 60. * 60.
    available_secs = max(30, max_secs)
    elapsed, start_t = 0., t.time()

    ptsampler = sampler.Sampler(
        nwalkers_ptm, ndim, loglkl, logp,
        loglargs=[fundam_params, synthcl_args, lkl_method, obs_clust, ranges,
                  varIdxs, priors_ptm], Tmax=Tmax, ntemps=ntemps)

    ntemps = ptsampler.ntemps
    # Initial population.
    p0 = initPop(
        ranges, varIdxs, lkl_method, obs_clust, fundam_params, synthcl_args,
        ntemps, nwalkers_ptm, init_mode_ptm, popsize_ptm, maxiter_ptm)

    print("     Burn-in stage")
    N_steps_check = max(1, int(nburn_ptm * .1))
    for i, (pos0, lnprob, lnlike) in enumerate(ptsampler.sample(
            p0, iterations=nburn_ptm, adapt=pt_adapt)):

        elapsed += t.time() - start_t
        if elapsed >= available_secs:
            print("  Time consumed during burn-in (runs={})".format(i + 1))
            break
        start_t = t.time()

        if (i + 1) % N_steps_check:
            continue
        maf = np.mean(ptsampler.acceptance_fraction[0])
        update_progress.updt(nburn_ptm, i + 1, "MAF={:.3f}".format(maf))

    # ptsampler.chain.shape: (ntemps, nwalkers, nsteps, ndim)
    # Store burn-in chain phase.
    # Shape: (runs, nwalkers, ndim)
    chains_nruns = ptsampler.chain[0].transpose(1, 0, 2)
    # The Mass parameter is not interpolated, use its grid values.
    chains_nruns = discreteParams(fundam_params, varIdxs, chains_nruns, [4])
    # print(chains_nruns.shape, pars_chains_bi.shape)
    pars_chains_bi = chains_nruns.T

    # Store MAP solution.
    idx_best = np.argmax(lnprob[0])
    map_sol_old = [
        fillParams(fundam_params, varIdxs, pos0[0][idx_best]),
        lnprob[0][idx_best]]

    ptsampler.reset()

    # We'll track how the average autocorrelation time estimate changes.
    # This will be useful to testing convergence.
    tau_index, autocorr_vals, old_tau = 0, np.empty(nsteps_ptm), np.inf

    # Check for convergence every 5% of steps or 100, whichever value
    # is lower.
    N_steps_conv = max(min(int(nsteps_ptm * .1), 100), 10)

    afs, tswaps = [], []
    # actimes = []
    # tau_emcee = np.empty(nsteps_ptm)

    # Check for dodgy inputs.
    if np.any(np.isinf(pos0)):
        print("  WARNING: At least one parameter value was infinite")
        # print(pos0)
        pos0[np.isinf(pos0)] = 0.
    if np.any(np.isnan(pos0)):
        print("  WARNING: At least one parameter value was NaN")
        pos0[np.isnan(pos0)] = 0.
        # print(pos0)

    maf_steps, prob_mean, map_lkl = [], [], []
    elapsed_in, start_in = 0., t.time()
    milestones = list(range(10, 101, 10))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        for i, result in enumerate(ptsampler.sample(
                pos0, iterations=nsteps_ptm, adapt=pt_adapt)):

            elapsed += t.time() - start_t
            if elapsed >= available_secs and i > 2:
                print("  Time consumed (runs={})".format(i + 1))
                break
            start_t = t.time()

            # Only check convergence every 'N_steps_conv' steps
            if (i + 1) % N_steps_conv:
                continue

            # Temperature swap acceptance fractions.
            tswaps.append(ptsampler.tswap_acceptance_fraction)
            # Mean acceptance fractions for all temperatures.
            afs.append(np.mean(ptsampler.acceptance_fraction, axis=1))

            # Compute the autocorrelation time using 'emcee'. tol=0 means that
            # we'll always get an estimate even if it isn't trustworthy.
            # tau0 = autocorr.integrated_time(
            #     ptsampler.chain[0, :, :i + 1, :].transpose(1, 0, 2), tol=0)
            # tau_emcee[tau_index] = np.mean(tau0)

            # Autocorrelation time for the non-tempered chain. Mean across
            # chains.
            # ptsampler.chain.shape: (ntemps, nwalkers, nsteps, ndim)
            # x.shape: (nsteps, ndim)
            x = np.mean(ptsampler.chain[0, :, :i + 1, :], axis=0)
            # tau.shape: ndim
            tau = util.autocorr_integrated_time(x)

            # acors = np.zeros(ntemps)
            # for temp in range(ntemps):
            #     x = np.mean(ptsampler.chain[temp, :, :i + 1, :], axis=0)
            #     acors[temp] = np.mean(util.autocorr_integrated_time(x))
            # actimes.append(acors)

            # Autocorrelation time. Mean across dimensions.
            autocorr_vals[tau_index] = np.mean(tau)
            tau_index += 1

            # Check convergence
            converged = np.all(tau * N_conv < (i + 1))
            converged &= np.all(np.abs(old_tau - tau) / tau < tol_conv)
            if converged:
                print("  Convergence achieved (runs={})".format(i + 1))
                break
            old_tau = tau

            pos, lnprob, lnlike = result

            maf = np.mean(ptsampler.acceptance_fraction[0])
            # Store MAP solution in this iteration.
            prob_mean.append([i, np.mean(lnprob[0])])
            idx_best = np.argmax(lnprob[0])
            # Update if a new optimal solution was found.
            if lnprob[0][idx_best] > map_sol_old[1]:
                map_sol_old = [
                    fillParams(fundam_params, varIdxs, pos[0][idx_best]),
                    lnprob[0][idx_best]]
            map_lkl.append([i, map_sol_old[1]])

            # Time used to check how fast the sampler is advancing.
            elapsed_in += t.time() - start_in
            start_in = t.time()
            # Print progress.
            percentage_complete = (100. * (i + 1) / nsteps_ptm)
            if len(milestones) > 0 and percentage_complete >= milestones[0]:
                map_sol, logprob = map_sol_old
                m, s = divmod(nsteps_ptm / (i / elapsed_in) - elapsed_in, 60)
                h, m = divmod(m, 60)
                print("{:>3}% ({:.3f}) LP={:.1f} ({:.5f}, {:.3f}, {:.3f}, "
                      "{:.2f}, {:.0f}, {:.2f})".format(
                          milestones[0], maf, logprob, *map_sol) +
                      " [{:.0f} m/s | {:.0f}h{:.0f}m]".format(
                          (ntemps * nwalkers_ptm * i) / elapsed_in, h, m))
                milestones = milestones[1:]

    runs = i + 1

    # Evolution of the mean autocorrelation time.
    tau_autocorr = autocorr_vals[:tau_index]
    # Mean acceptance fractions for all replicas.
    maf_steps = (
        N_steps_conv * np.arange(1, tau_index + 1), np.asarray(afs).T)
    # Temperature swaps acceptance fractions.
    tswaps_afs = (
        N_steps_conv * np.arange(1, tau_index + 1), np.asarray(tswaps).T)
    # Betas history
    betas_pt = (np.arange(runs), ptsampler.beta_history)

    # Final MAP fit.
    map_sol, map_lkl_final = map_sol_old

    # This number should be between approximately 0.25 and 0.5 if everything
    # went as planned.
    m_accpt_fr = np.mean(ptsampler.acceptance_fraction[0])
    if m_accpt_fr > .5 or m_accpt_fr < .25:
        print("  WARNING: mean acceptance fraction is outside of the\n"
              "  recommended range.")

    # ptsampler.chain.shape: (ntemps, nchains, nsteps, ndim)
    # chains_nruns.shape: (runs, nchains, ndim)
    chains_nruns = ptsampler.chain[0, :, :runs, :].transpose(1, 0, 2)
    chains_nruns = discreteParams(fundam_params, varIdxs, chains_nruns, [4])

    # Convergence parameters.
    acorr_t, med_at_c, all_taus, geweke_z, lag_zero, acorr_function,\
        mcmc_ess = convergenceVals(
            'ptemcee', ndim, varIdxs, N_conv, chains_nruns)

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

    # Push Mass value to grid value for mean, map, and mode solutions.
    mean_sol = closeSol(fundam_params, mean_sol, [4])
    mode_sol = closeSol(fundam_params, mode_sol, [4])
    map_sol = closeSol(fundam_params, map_sol, [4])

    # Total number of values used to estimate the parameter's distributions.
    N_total = mcmc_trace.shape[-1]

    isoch_fit_params = {
        'varIdxs': varIdxs, 'nsteps_ptm': runs, 'mean_sol': mean_sol,
        'median_sol': median_sol, 'map_sol': map_sol, 'map_lkl': map_lkl,
        'mode_sol': mode_sol, 'pardist_kde': pardist_kde, 'param_r2': param_r2,
        'map_lkl_final': map_lkl_final, 'prob_mean': prob_mean,
        'bf_elapsed': elapsed, 'mcmc_trace': mcmc_trace,
        'pars_chains_bi': pars_chains_bi, 'pars_chains': chains_nruns.T,
        'maf_steps': maf_steps, 'tswaps_afs': tswaps_afs, 'betas_pt': betas_pt,
        'acorr_t': acorr_t, 'med_at_c': med_at_c,
        'all_taus': all_taus,
        # 'max_at_c': max_at_c, 'min_at_c': min_at_c,
        # 'minESS': minESS, 'mESS': mESS, 'mESS_epsilon': mESS_epsilon,
        'lag_zero': lag_zero, 'acorr_function': acorr_function,
        'geweke_z': geweke_z, 'mcmc_ess': mcmc_ess, 'N_total': N_total,
        'N_steps_conv': N_steps_conv, 'N_conv': N_conv, 'tol_conv': tol_conv,
        'tau_index': tau_index, 'tau_autocorr': tau_autocorr
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
        synth_clust = synthClust(fundam_params, varIdxs, model, synthcl_args)
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


import numpy as np
import time as t
import warnings

from ..synth_clust import synth_cluster
from . import likelihood
from .mcmc_convergence import convergenceVals
from .bf_common import initPop, varPars, rangeCheck, fillParams,\
    r2Dist, modeKDE  # , thinChain


def main(
    err_lst, completeness, e_max, max_mag_syn, obs_clust, ext_coefs,
    st_dist_mass, N_fc, m_ini, err_rnd, lkl_method, fundam_params,
    theor_tracks, R_V, nsteps_mcee, nwalkers_mcee, nburn_mcee, priors_mcee,
        emcee_moves, hmax, **kwargs):
    """
    """
    from emcee import ensemble
    # This is used when the moves are defined below by eval()
    from emcee import moves

    varIdxs, ndim, ranges = varPars(fundam_params)
    # Pack synthetic cluster arguments.
    synthcl_args = [
        theor_tracks, e_max, err_lst, completeness, max_mag_syn, st_dist_mass,
        R_V, ext_coefs, N_fc, m_ini, err_rnd]

    # Start timing.
    max_secs = hmax * 60. * 60.
    available_secs = max(30, max_secs)

    # Define moves.
    mv = [(eval("moves." + _)) for _ in emcee_moves]
    # Define sampler.
    sampler = ensemble.EnsembleSampler(
        nwalkers_mcee, ndim, log_posterior,
        args=[priors_mcee, varIdxs, ranges, fundam_params, synthcl_args,
              lkl_method, obs_clust], moves=mv)

    ntemps = 1
    pos0 = initPop(
        ranges, varIdxs, lkl_method, obs_clust, fundam_params, synthcl_args,
        ntemps, nwalkers_mcee, 'random', None, None)
    pos0 = pos0[0]

    maf_steps, prob_mean, tau_autocorr = [], [], []
    map_lkl, map_sol_old = [], [[], -np.inf]

    # We'll track how the average autocorrelation time estimate changes.
    # This will be useful to testing convergence.
    # tau_index, autocorr_vals, old_tau = 0, np.empty(nsteps), np.inf

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        N_steps_store, runs = 50, 0

        elapsed, start = 0., t.time()
        milestones = list(range(10, 101, 10))
        for i, (pos, prob, stat) in enumerate(
                sampler.sample(pos0, iterations=nsteps_mcee)):

            # Only check convergence every 'N_steps_store' steps
            if (i + 1) % N_steps_store:
                continue
            runs += 1

            # Compute the autocorrelation time so far. Using tol=0 means that
            # we'll always get an estimate even if it isn't trustworthy.
            tau = sampler.get_autocorr_time(tol=0)
            tau_autocorr.append([i, np.mean(tau)])

            # # Check convergence
            # converged = np.all(tau * N_conv < (i + 1))
            # converged &= np.all(np.abs(old_tau - tau) / tau < tol_conv)
            # autocorr_vals[tau_index] = np.nanmean(tau)
            # tau_index += 1

            maf = np.mean(sampler.acceptance_fraction)
            maf_steps.append(maf)

            # Store MAP solution in this iteration.
            prob_mean.append(np.mean(prob))
            idx_best = np.argmax(prob)
            # Update if a new optimal solution was found.
            if prob[idx_best] > map_sol_old[1]:
                map_sol_old = [
                    fillParams(fundam_params, varIdxs, pos[idx_best]),
                    prob[idx_best]]
            map_lkl.append(map_sol_old[1])

            # Time used to check how fast the sampler is advancing.
            elapsed += t.time() - start
            start = t.time()
            # Print progress.
            percentage_complete = (100. * (i + 1) / nsteps_mcee)
            if len(milestones) > 0 and percentage_complete >= milestones[0]:
                map_sol, logprob = map_sol_old
                m, s = divmod(nsteps_mcee / (i / elapsed) - elapsed, 60)
                h, m = divmod(m, 60)
                print("{:>3}% ({:.3f}) LP={:.1f} ({:.5f}, {:.3f}, {:.3f}, "
                      "{:.2f}, {:.0f}, {:.2f})".format(
                          milestones[0], maf, logprob, *map_sol) +
                      " [{:.0f} m/s | {:.0f}h{:.0f}m]".format(
                          (ntemps * nwalkers_mcee * i) / elapsed, h, m))
                milestones = milestones[1:]

            # Stop when available time is consumed.
            if elapsed >= available_secs:
                print("  Time consumed")
                break

    # Total number of steps
    N_steps = N_steps_store * np.arange(1, runs + 1)

    # Final MAP fit.
    map_sol, map_lkl_final = map_sol_old

    # This number should be between approximately 0.25 and 0.5 if everything
    # went as planned.
    m_accpt_fr = np.mean(sampler.acceptance_fraction)
    if m_accpt_fr > .5 or m_accpt_fr < .25:
        print("  WARNING: mean acceptance fraction is outside of the\n"
              "  recommended range.")

    # all_chains.shape = (N_steps_store * runs, nchains, ndims)
    all_chains = sampler.get_chain()
    # Store burn-in chain phase.
    bi_steps = int(nburn_mcee * all_chains.shape[0])
    # chains_nruns.shape: (bi_steps, nchains, ndim)
    chains_nruns = all_chains[:bi_steps]
    # pars_chains_bi.shape: (ndim, nchains, bi_steps)
    pars_chains_bi = chains_nruns.T

    # After burn-in
    chains_nruns = all_chains[bi_steps:]
    pars_chains = chains_nruns.T

    # Convergence parameters.
    tau_autocorr = np.array(tau_autocorr).T
    _, acorr_t, med_at_c, all_taus, geweke_z, acorr_function,\
        mcmc_ess = convergenceVals(
            'emcee', ndim, varIdxs, chains_nruns, bi_steps)

    # Re-shape trace for all parameters (flat chain).
    # Shape: (ndim, runs * nchains)
    mcmc_trace = chains_nruns.reshape(-1, ndim).T

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
        'varIdxs': varIdxs, 'mean_sol': mean_sol,
        'median_sol': median_sol, 'map_sol': map_sol, 'map_lkl': map_lkl,
        'mode_sol': mode_sol, 'pardist_kde': pardist_kde, 'param_r2': param_r2,
        'map_lkl_final': map_lkl_final, 'prob_mean': prob_mean,
        'bf_elapsed': elapsed, 'mcmc_trace': mcmc_trace,
        'pars_chains_bi': pars_chains_bi, 'pars_chains': pars_chains,
        'maf_steps': maf_steps,
        #
        'tau_autocorr': tau_autocorr, 'acorr_t': acorr_t, 'med_at_c': med_at_c,
        'all_taus': all_taus, 'acorr_function': acorr_function,
        # 'max_at_c': max_at_c, 'min_at_c': min_at_c,
        # 'minESS': minESS, 'mESS': mESS, 'mESS_epsilon': mESS_epsilon,
        'geweke_z': geweke_z, 'mcmc_ess': mcmc_ess, 'N_total': N_total,
        'N_steps': N_steps
    }

    return isoch_fit_params


def log_posterior(
    model, priors, varIdxs, ranges, fundam_params, synthcl_args, lkl_method,
        obs_clust):
    """
    Log posterior function.
    """
    lp = log_prior(model, priors, varIdxs, ranges)
    if not np.isfinite(lp):
        return -np.inf
    lkl = log_likelihood(
        model, varIdxs, fundam_params, synthcl_args, lkl_method, obs_clust)
    return lp + lkl


def log_prior(model, priors, varIdxs, ranges):
    """
    met, age, ext, dm, mass, bf = model
    """
    rangeFlag = rangeCheck(model, ranges, varIdxs)

    lp = -np.inf
    # If some parameter is outside of the given ranges, don't bother obtaining
    # the proper model and just pass -inf
    if rangeFlag:
        lp = 0.
        for i, pr in enumerate(priors):
            # Gaussian prior. If the prior is uniform, simply pass.
            if pr[0] == 'g':
                lp += np.log(1 / pr[2]) - .5 * np.square(
                    (model[i] - pr[1]) / pr[2])

    return lp


def log_likelihood(
        model, varIdxs, fundam_params, synthcl_args, lkl_method, obs_clust):
    """
    The Dolphin likelihood needs to be *minimized*. Be careful with the priors.
    """

    # Generate synthetic cluster.
    synth_clust = synth_cluster.main(
        fundam_params, varIdxs, model, *synthcl_args)

    # Call likelihood function for this model. RETURNS THE INVERSE lkl.
    lkl = likelihood.main(lkl_method, synth_clust, obs_clust)

    # TODO, is this a general approach?
    # The negative likelihood is returned since Dolphin requires a minimization
    # of the PLR, and here we are maximizing
    return -lkl

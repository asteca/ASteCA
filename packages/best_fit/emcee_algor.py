
import numpy as np
import random
import logging
import time as t

# TODO emcee3
# import emcee
from .emcee3rc1 import ensemble
from .emcee3rc1 import moves
from .emcee3rc1 import utils
from .emcee3rc1 import autocorr

from ..synth_clust import synth_cluster
from . import likelihood
from .mcmc_convergence import multiESS, fminESS, geweke, effective_n,\
    pdfHalfves
from .. import update_progress


def main(
    lkl_method, e_max, err_lst, completeness, max_mag_syn, fundam_params,
        obs_clust, theor_tracks, R_V, ext_coefs, st_dist_mass, N_fc, cmpl_rnd,
        err_rnd, nwalkers, nsteps, nburn, N_burn, emcee_a, priors, *args):
    """

    nwalkers: number of MCMC walkers
    nwalkers: number of MCMC steps to take
    nburn: "burn-in" period to let chains stabilize

    """

    varIdxs, ndim, ranges = varPars(fundam_params)

    synthcl_args = [
        theor_tracks, e_max, err_lst, completeness, max_mag_syn, st_dist_mass,
        R_V, ext_coefs, N_fc, cmpl_rnd, err_rnd]

    # TODO make this a proper parameter
    if emcee_a <= 0.:
        mv = moves.KDEMove()
    elif 0 <= emcee_a <= 2.:
        mv = moves.DEMove(.1, nsplits=2)  # , gamma0=emcee_a)
    else:
        mv = moves.DESnookerMove()

    # Define sampler.
    model_done = {}  # TODO remove?
    # TODO emcee3: pass selected moves and its parameters
    sampler = ensemble.EnsembleSampler(
        nwalkers, ndim, log_posterior,
        args=[priors, varIdxs, ranges, fundam_params, synthcl_args, lkl_method,
              obs_clust, model_done], moves=mv)

    # Burn-in period.
    t0 = t.time()
    pos, best_sol_old, pars_chains_bi = burnIn(
        nwalkers, nburn, N_burn, fundam_params, varIdxs, sampler,
        ranges, priors, synthcl_args, lkl_method, obs_clust)
    # Reset the chain to remove the burn-in samples.
    sampler.reset()

    # TODO add this parameter to the input params file
    max_secs = 22. * 60. * 60.
    elapsed = t.time() - t0
    available_secs = max(30, max_secs - elapsed)

    s = t.time()
    # We'll track how the average autocorrelation time estimate changes
    tau_index, autocorr_vals = 0, np.empty(nsteps)
    # This will be useful to testing convergence
    old_tau = np.inf

    # Check for convergence every 10% of steps or 100, whichever value
    # is lower.
    N_steps_conv = min(int(nsteps / 10.), 100)
    # TODO emcee3 input as params
    N_conv, tol_conv = 100., 0.01

    maf_steps, map_lkl = [], []
    milestones = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    for i, result in enumerate(sampler.sample(pos, iterations=nsteps)):

        # Only check convergence every 'N_steps_conv' steps
        if (i + 1) % N_steps_conv:
            continue

        # Compute the autocorrelation time so far. Using tol=0 means that
        # we'll always get an estimate even if it isn't trustworthy.
        try:
            tau = sampler.get_autocorr_time(tol=0)
            autocorr_vals[tau_index] = np.nanmean(tau)
            tau_index += 1

            # Check convergence
            converged = np.all(tau * N_conv < (i + 1))
            converged &= np.all(np.abs(old_tau - tau) / tau < tol_conv)
            if converged:
                print("  Convergence achieved.")
                break
            old_tau = tau
        except FloatingPointError:
            pass

        elapsed += t.time() - s
        if elapsed >= available_secs:
            print("  Time consumed.")
            break

        pos, prob, state = result

        maf = np.mean(sampler.acceptance_fraction)
        maf_steps.append([i, maf])

        # Discard -np.inf posterior values.
        idx_keep = prob > -np.inf
        # Store MAP solution in this iteration.
        try:
            idx_best = np.argmin(-prob[idx_keep])
        except ValueError:
            idx_best = 0

        # Update if a new optimal solution was found.
        if -prob[idx_best] < best_sol_old[1]:
            best_sol_old = [
                closeSol(fundam_params, varIdxs, pos[idx_best]),
                -prob[idx_best]]
        map_lkl.append([i, best_sol_old[1]])

        # Print progress.
        percentage_complete = (100. * (i + 1) / nsteps)
        if len(milestones) > 0 and percentage_complete >= milestones[0]:
            map_sol, logprob = best_sol_old
            m, s = divmod(nsteps / (i / elapsed) - elapsed, 60)
            h, m = divmod(m, 60)
            print(" {:>3}% ({:.3f}) LP={:.1f} ({:g}, {:g}, {:.3f}, {:.2f}"
                  ", {:g}, {:.2f})".format(
                      milestones[0], maf, logprob, *map_sol) +
                  " [{:.0f} m/s | {:.0f}h{:.0f}m]".format(
                      (nwalkers * i) / elapsed, h, m))
            milestones = milestones[1:]

        s = t.time()
    runs = i + 1

    # Evolution of the mean autocorrelation time.
    tau_autocorr = autocorr_vals[:tau_index]

    # Final MAP fit.
    idx_keep = prob > -np.inf
    try:
        idx_best = np.argmin(-prob[idx_keep])
    except ValueError:
        idx_best = 0
        logging.warning(
            " No valid solution could be found. Run a longer chain\n"
            "or modify the sampler parameters.")
    map_sol = closeSol(fundam_params, varIdxs, pos[idx_best])
    map_lkl_final = -prob[idx_best]

    # This number should be between approximately 0.25 and 0.5 if everything
    # went as planned.
    m_accpt_fr = np.mean(sampler.acceptance_fraction)
    if m_accpt_fr > .5 or m_accpt_fr < .25:
        print("  WARNING: mean acceptance fraction is outside of the\n"
              "  recommended range.")

    # Use only the sampled runs (necessary if there was a time break out)
    # Shape: (runs, nwalkers, ndim)
    chains_nruns = sampler.get_chain()[:runs, :, :]
    # Change values for the discrete parameters with the closest valid values.
    chains_nruns = discreteParams(fundam_params, varIdxs, chains_nruns)
    # Re-shape trace for all parameters (flat chain).
    # Shape: (ndim, runs * nwalkers)
    emcee_trace = chains_nruns.reshape(-1, ndim).T

    # Convergence parameters.
    acorr_t, max_at_5c, min_at_5c, geweke_z, emcee_acorf, pymc3_ess, minESS,\
        mESS, mESS_epsilon, mcmc_halves = convergenceVals(
            ndim, varIdxs, N_conv, chains_nruns, emcee_trace)

    # Pass the mean as the best model fit found.
    best_sol = closeSol(fundam_params, varIdxs, np.mean(emcee_trace, axis=1))

    isoch_fit_params = {
        'varIdxs': varIdxs, 'nsteps': runs, 'best_sol': best_sol,
        'map_sol': map_sol, 'map_lkl': map_lkl, 'map_lkl_final': map_lkl_final,
        'mcmc_elapsed': elapsed, 'mcmc_trace': emcee_trace,
        'pars_chains_bi': pars_chains_bi, 'pars_chains': chains_nruns.T,
        'maf_steps': maf_steps, 'autocorr_time': acorr_t,
        'max_at_5c': max_at_5c, 'min_at_5c': min_at_5c,
        'minESS': minESS, 'mESS': mESS, 'mESS_epsilon': mESS_epsilon,
        'emcee_acorf': emcee_acorf, 'geweke_z': geweke_z,
        'pymc3_ess': pymc3_ess, 'mcmc_halves': mcmc_halves,
        'N_steps_conv': N_steps_conv, 'N_conv': N_conv, 'tol_conv': tol_conv,
        'tau_index': tau_index, 'tau_autocorr': tau_autocorr
    }

    return isoch_fit_params


def varPars(fundam_params):
    """
    Check which parameters are fixed and which have a dynamic range. This also
    dictates the number of free parameters.
    """
    ranges = [[min(_), max(_)] for _ in fundam_params]

    varIdxs = []
    for i, r in enumerate(ranges):
        if r[0] != r[1]:
            varIdxs.append(i)

    # Number of free parameters.
    ndim = len(varIdxs)

    return varIdxs, ndim, np.array(ranges)


def burnIn(
    nwalkers, nburn, N_burn, fundam_params, varIdxs, sampler,
        ranges, priors, synthcl_args, lkl_method, obs_clust):
    """
    N_burn: number of times the burn-in phase will be applied, focusing on a
    small ball around the MAP.
    """
    # Random initial models.
    starting_guesses = random_population(fundam_params, varIdxs, nwalkers)

    # starting_guesses = np.asarray(
    #     isoch_fit_params[-1][0][:nwalkers]).T[varIdxs].T

    # starting_guesses[0] = [0.0535, 9.85, .27, 12.6, 5000., .5]

    N_total, N_done = nburn * N_burn, 0
    print("     Burn-in stage")
    for _ in range(N_burn):

        for i, result in enumerate(
                sampler.sample(starting_guesses, iterations=nburn)):
            update_progress.updt(N_total, N_done + i + 1)
        N_done += nburn

        pos, prob, state = result

        # Discard -np.inf posterior values.
        idx_keep = prob > -np.inf
        # Best solution (MAP) in this iteration.
        try:
            idx_best = np.argmin(-prob[idx_keep])
        except ValueError:
            # All values were -inf, so there are no left in prob[idx_keep].
            idx_best = 0

        # TODO change when emcee3 is properly imported
        starting_guesses = utils.sample_ball(
            pos[idx_best], pos.std(axis=0), size=nwalkers)

        # best_sol = [pos[idx_best], -prob[idx_best]]
        # # Small ball (N sigma) around the best solution.
        # mean_p, std_p = pos.mean(axis=0), pos.std(axis=0)
        # low, high = mean_p - 1. * std_p, mean_p + 1. * std_p
        # starting_guesses = []
        # for p in pos:
        #     # Check that all elements in 'p' are between the (low,high)
        #     # values.
        #     flag = True
        #     for i, v in enumerate(p):
        #         # If any element is not inside its range, break out.
        #         if v <= low[i] or v >= high[i]:
        #             flag = False
        #             break

        #     if flag is True:
        #         starting_guesses.append(p)
        #     else:
        #         # If at least one element was not inside its range, store the
        #         # fill_value
        #         starting_guesses.append(best_sol[0])
        # starting_guesses = np.asarray(starting_guesses)

    # Store burn-in chain phase.
    chains_nruns = sampler.get_chain()[-nburn:, :, :]
    pars_chains_bi = discreteParams(fundam_params, varIdxs, chains_nruns).T

    # Discard -np.inf posterior values.
    idx_keep = prob > -np.inf
    # Store MAP solution.
    try:
        idx_best = np.argmin(-prob[idx_keep])
    except ValueError:
        # All values were -inf, so there are no left in prob[idx_keep].
        idx_best = 0
    best_sol = [pos[idx_best], -prob[idx_best]]

    return pos, best_sol, pars_chains_bi


def random_population(fundam_params, varIdxs, n_ran):
    """
    Generate a random set of parameter values to use as a random population.

    Pick n_ran initial random solutions from each list storing all the
    possible parameters values.
    """
    p_lst = []
    for i in varIdxs:
        p_lst.append([random.choice(fundam_params[i]) for _ in range(n_ran)])

    return np.array(p_lst).T


def log_posterior(
    model, priors, varIdxs, ranges, fundam_params, synthcl_args, lkl_method,
        obs_clust, model_done):
    """
    Log posterior function.
    """
    lp, model_proper = log_prior(model, priors, fundam_params, varIdxs, ranges)
    if not np.isfinite(lp):
        return -np.inf
    lkl, model_done = log_likelihood(
        model_proper, fundam_params, synthcl_args, lkl_method, obs_clust,
        model_done)
    return lp + lkl


def log_prior(model, priors, fundam_params, varIdxs, ranges):
    """
    met, age, ext, dm, mass, bf = model
    """
    check_ranges = [
        r[0] <= p <= r[1] for p, r in zip(*[model, ranges[varIdxs]])]

    lp, model_proper = -np.inf, []
    # If some parameter is outside of the given ranges, don't bother obtaining
    # the proper model and just pass -inf
    if all(check_ranges):
        model_proper = closeSol(fundam_params, varIdxs, model)
        # if model_proper:

        if priors == 'unif':
            # Flat prior
            lp = 0.

        elif priors == 'gauss':
            # TODO finish this
            model_mean = closeSol(
                fundam_params, varIdxs,
                [np.mean(fundam_params[_]) for _ in varIdxs])
            # Gaussian prior.
            model_std = np.array([0.005, .3, .1, .2, 500., 0.3])
            lp = np.sum(-np.square(
                (np.asarray(model_proper) - model_mean) / model_std))
        # else:
        #     return lp, model_proper

    return lp, model_proper


def closeSol(fundam_params, varIdxs, model):
    """
    Find the closest value in the parameters list for the discrete parameters
    metallicity, age, and mass.
    """
    model_proper, j = [], 0
    for i, par in enumerate(fundam_params):
        # If this parameter is one of the 'free' parameters.
        if i in varIdxs:
            # If it is the parameter metallicity, age or mass.
            if i in [0, 1, 4]:
                # Select the closest value in the array of allowed values.
                # pp = min(par, key=lambda x: abs(x - model[i - j]))
                # if abs(pp - model[i - j]) < model[i - j] * .01:
                #     model_proper.append(pp)
                model_proper.append(min(
                    par, key=lambda x: abs(x - model[i - j])))
                # else:
                #     return []

            else:
                model_proper.append(model[i - j])
        else:
            model_proper.append(par[0])
            j += 1

    # Store rounded parameter values to make the 'model_done'
    # key assignment more stable (else there could be a large
    # number of decimals)

    # return np.round(model_proper, 5)
    # TODO required for Autograd to work when using sampyl
    # return np.array([np.round(_, 5) for _ in model_proper])
    return model_proper


def log_likelihood(
    model_proper, fundam_params, synthcl_args, lkl_method, obs_clust,
        model_done):
    """
    The Dolphin likelihood needs to be *minimized*. Be careful with the priors.
    """

    # TODO if 'm_sample' is True, always process the synthetic cluster and
    # the likelihood

    # # If the model was already processed, use the generated likelihood value.
    # try:
    #     lkl = model_done[''.join(map(str, model_proper))]
    # except KeyError:

    theor_tracks, e_max, err_lst, completeness, max_mag_syn, st_dist_mass,\
        R_V, ext_coefs, N_fc, cmpl_rnd, err_rnd = synthcl_args

    # Metallicity and age indexes to identify isochrone.
    m_i = fundam_params[0].index(model_proper[0])
    a_i = fundam_params[1].index(model_proper[1])
    isochrone = theor_tracks[m_i][a_i]

    # Generate synthetic cluster.
    synth_clust = synth_cluster.main(
        e_max, err_lst, completeness, max_mag_syn, st_dist_mass, isochrone,
        R_V, ext_coefs, N_fc, cmpl_rnd, err_rnd, model_proper)

    # import matplotlib.pyplot as plt

    # Call likelihood function for this model.
    lkl = likelihood.main(lkl_method, synth_clust, obs_clust)

        # # Store processed model and its likelihood value.
        # model_done[''.join(map(str, model_proper))] = lkl

    # TODO, is this a general approach?
    # The negative likelihood is returned since Dolphin requires a minimization
    # of the PLR, and here we are maximizing
    return -lkl, model_done


def discreteParams(fundam_params, varIdxs, chains_nruns):
    """
    Push values in each chain for each discrete parameter to the closest
    accepted value.

    chains_nruns.shape: (runs, nwalkers, ndim)
    """
    params = []
    for i, par in enumerate(fundam_params):
        p = np.array(par)
        # If this parameter is one of the 'free' parameters.
        if i in varIdxs:
            # If it is the parameter metallicity, age or mass.
            if i in [0, 1, 4]:
                pc = chains_nruns.T[i]
                chains = []
                for c in pc:
                    chains.append(
                        p[abs(c[None, :] - p[:, None]).argmin(axis=0)])
                params.append(chains)
            else:
                params.append(chains_nruns.T[i])

    return np.array(params).T


def convergenceVals(ndim, varIdxs, N_conv, chains_nruns, emcee_trace):
    """
    Convergence statistics.
    """

    # Autocorrelation time for each parameter.
    acorr_t = autocorr.integrated_time(chains_nruns, tol=N_conv, quiet=True)

    # Autocorrelation time for each chain for each parameter.
    logger = logging.getLogger()
    logger.disabled = True
    at = []
    for p in chains_nruns.T:
        at_p = []
        for c in p:
            at_p.append(autocorr.integrated_time(c, quiet=True)[0])
        at.append(at_p)
    logger.disabled = False

    # Select the indexes of the 5 chains with the largest acorr times, and
    # the 5 chains with the smallest acorr times, for each parameter.
    if len(at[0]) >= 5:
        max_at_5c = [np.argpartition(a, -5)[-5:] for a in at]
        min_at_5c = [np.argpartition(a, 5)[:5] for a in at]
    else:
        max_at_5c, min_at_5c = [np.array([0])] * ndim, [np.array([0])] * ndim

    # Worst chain: chain with the largest acorr time.
    # max_at_c = [np.argmax(a) for a in at]

    # Mean Geweke z-scores and autocorrelation functions for all chains.
    geweke_z, emcee_acorf = [[] for _ in range(ndim)],\
        [[] for _ in range(ndim)]
    for i, p in enumerate(chains_nruns.T):
        for c in p:
            try:
                geweke_z[i].append(geweke(c))  # p[max_at_c[i]]
            except ZeroDivisionError:
                geweke_z[i].append([np.nan, np.nan])
            try:
                emcee_acorf[i].append(autocorr.function_1d(c))
            except FloatingPointError:
                emcee_acorf[i].append([np.nan])
    geweke_z = np.nanmean(geweke_z, axis=1)
    emcee_acorf = np.nanmean(emcee_acorf, axis=1)

    # PyMC3 effective sample size.
    try:
        # Change shape to (nchains, nstesp, ndim)
        pymc3_ess = effective_n(chains_nruns.transpose(1, 0, 2))
    except FloatingPointError:
        pymc3_ess = np.array([np.nan] * ndim)

    # TODO fix this function
    # Minimum effective sample size (ESS), and multi-variable ESS.
    minESS, mESS = fminESS(ndim), multiESS(emcee_trace.T)
    mESS_epsilon = [[], [], []]
    for alpha in [.01, .05, .1, .2, .3, .4, .5, .6, .7, .8, .9, .95]:
        mESS_epsilon[0].append(alpha)
        mESS_epsilon[1].append(fminESS(ndim, alpha=alpha, ess=minESS))
        mESS_epsilon[2].append(fminESS(ndim, alpha=alpha, ess=mESS))

    mcmc_halves = pdfHalfves(varIdxs, emcee_trace)

    return acorr_t, max_at_5c, min_at_5c, geweke_z, emcee_acorf, pymc3_ess,\
        minESS, mESS, mESS_epsilon, mcmc_halves

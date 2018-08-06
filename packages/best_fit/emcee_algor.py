
import numpy as np
import random
import emcee
from ..synth_clust import synth_cluster
from . import likelihood
from .mcmc_convergence import multiESS, fminESS, geweke, effective_n
from .emcee_new_AT import integrated_time
from .. import update_progress

import time as t


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

    # Generate sampler. TODO do something with 'a' parameter
    model_done = {}
    sampler = emcee.EnsembleSampler(
        nwalkers, ndim, log_posterior,
        args=[priors, varIdxs, ranges, fundam_params, synthcl_args, lkl_method,
              obs_clust, model_done], a=emcee_a)

    # Burn-in period.
    t0 = t.time()
    pos, prob, state, best_sol_old, pars_chains_bi = burnIn(
        nwalkers, nburn, N_burn, fundam_params, varIdxs, sampler)

    # Reset the chain to remove the burn-in samples.
    sampler.reset()

    # # Obtain the number of steps given the time constrain.
    # try:
    #     nsteps = int(float(nsteps))
    #     flag_time_run, total = False, nsteps
    # except ValueError:
    #     flag_time_run = True
    #     max_secs = float(nsteps[:-1]) * 60. * 60.
    #     available_secs = max_secs - (t.time() - t0)
    #     print(available_secs)
    #     total = available_secs + 30.
    #     if available_secs < 10:
    #         print("  WARNING: the maximum allowed time has been consumed.\n"
    #               "  Running with nsteps=10.")
    #         nsteps = 10
    #     else:
    #         import os
    #         free_ram_byte = float(os.popen(
    #             'free -t -m').readlines()[-1].split()[1:][-1]) * 1024. * 1024.
    #         # Use 75% of available mem
    #         print(free_ram_byte/1024./1024.)
    #         nsteps = int((.75 * free_ram_byte) / (nwalkers * ndim * 8.))
    #         print(nsteps)
    #         nsteps = 100000

    nsteps = int(float(nsteps))
    total = nsteps
    max_secs = 7. * 60. * 60.
    available_secs = max_secs - (t.time() - t0)
    print(nsteps, available_secs)
    # total = available_secs + 30.
    if available_secs < 10:
        print("  WARNING: the maximum allowed time has been consumed.\n"
              "  Running with nsteps=10.")
        nsteps = 10

    elapsed = t.time() - t0

    s = t.time()
    milestones = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    # Run sampler
    for i, result in enumerate(sampler.sample(
            pos, lnprob0=prob, rstate0=state, iterations=nsteps)):

        # if flag_time_run:
        elapsed += t.time() - s
        if elapsed >= available_secs:
            print("time consumed")
            break

        # Print progress.
        percentage_complete = (100. * (i + 1) / total)
        if len(milestones) > 0 and percentage_complete >= milestones[0]:

            # Discard -np.inf posterior values.
            idx_keep = result[1] > -np.inf
            # Store MAP solution in this iteration.
            idx_best = np.argmin(-result[1][idx_keep])

            best_sol_run = [result[0][idx_best], -result[1][idx_best]]
            # Check if a new optimal solution was found.
            if best_sol_run[1] < best_sol_old[1]:
                map_sol = closeSol(fundam_params, varIdxs, best_sol_run[0])
                logprob = best_sol_run[1]
                best_sol_old = [best_sol_run[0], best_sol_run[1]]
            else:
                map_sol = closeSol(fundam_params, varIdxs, best_sol_old[0])
                logprob = best_sol_old[1]
            maf = np.mean(sampler.acceptance_fraction)
            print(" {:>3}% ({:.3f}) LP={:.1f} ({:g}, {:g}, {:.3f}, {:.2f}"
                  ", {:g}, {:.2f})".format(
                      milestones[0], maf, logprob, *map_sol))
            milestones = milestones[1:]

        s = t.time()
    runs = i + 1

    print(elapsed, runs)
    # Final MAP fit.
    idx_keep = result[1] > -np.inf
    idx_best = np.argmin(-result[1][idx_keep])
    map_sol = closeSol(fundam_params, varIdxs, result[0][idx_best])

    # This number should be between approximately 0.25 and 0.5 if everything
    # went as planned.
    m_accpt_fr = np.mean(sampler.acceptance_fraction)
    # print("  Mean acceptance fraction: {:.3f}".format(m_accpt_fr))
    if m_accpt_fr > .5 or m_accpt_fr < .25:
        print("  WARNING: mean acceptance fraction is outside of the\n"
              "  recommended range.")

    chains_nruns = sampler.chain[:, :runs, :]
    # Re-shape trace for all parameters (flat chain).
    emcee_trace = chains_nruns.reshape(-1, ndim).T

    # Convergence parameters.
    acorr_t, max_at_10c, geweke_z, emcee_acorf, pymc3_ess, minESS,\
        mESS, mESS_epsilon = convergenceVals(ndim, chains_nruns, emcee_trace)

    # Pass the mean as the best model fit found.
    best_sol = closeSol(fundam_params, varIdxs, np.mean(emcee_trace, axis=1))

    isoch_fit_params = {
        'nsteps': runs, 'best_sol': best_sol, 'map_sol': map_sol,
        'pars_chains_bi': pars_chains_bi, 'pars_chains': chains_nruns.T,
        'm_accpt_fr': m_accpt_fr, 'varIdxs': varIdxs,
        'mcmc_trace': emcee_trace, 'autocorr_time': acorr_t,
        'max_at_10c': max_at_10c, 'minESS': minESS, 'mESS': mESS,
        'mESS_epsilon': mESS_epsilon, 'emcee_acorf': emcee_acorf,
        'geweke_z': geweke_z, 'pymc3_ess': pymc3_ess
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


def burnIn(nwalkers, nburn, N_burn, fundam_params, varIdxs, sampler):
    """
    N_burn: number of times the burn-in phase will be applied, focusing on a
    small ball around the MAP.
    """
    # Random initial models.
    starting_guesses = random_population(fundam_params, varIdxs, nwalkers)

    # # GA initial models.
    # N_pop, N_gen, fit_diff, cross_prob, cross_sel, mut_prob, N_el,\
    #     N_ei, N_es = 10, 100, 1., 0.85, '2P', 0.05, 5, 50, 20
    # # total sols: (N_pop * N_gen) + N_pop
    # if (N_pop * N_gen) + N_pop < nwalkers:
    #     N_pop = int(nwalkers / (1 + N_gen)) + 1

    # import genetic_algorithm
    # isoch_fit_params = genetic_algorithm.main(
    #     lkl_method, e_max, err_lst, completeness, max_mag_syn,
    #     fundam_params, obs_clust, theor_tracks, R_V, ext_coefs,
    #     st_dist_mass, N_fc, N_pop, N_gen, fit_diff, cross_prob,
    #     cross_sel, mut_prob, N_el, N_ei, N_es, False)

    # starting_guesses = np.asarray(
    #     isoch_fit_params[-1][0][:nwalkers]).T[varIdxs].T

    np.seterr(all='raise')

    N_total, N_done = nburn * N_burn, 0
    # Run burn-in.
    print("     Burn-in stage")
    for _ in range(N_burn):

        try:
            for i, result in enumerate(
                    sampler.sample(starting_guesses, iterations=nburn)):
                update_progress.updt(N_total, N_done + i + 1)
            N_done += nburn
        except FloatingPointError:
            print("Fp error")
            import pdb; pdb.set_trace()  # breakpoint c3f0d703 //

        pos, prob, state = result

        # Discard -np.inf posterior values.
        idx_keep = prob > -np.inf
        # Store best solution in this iteration.
        idx_best = np.argmin(-prob[idx_keep])

        # starting_guesses = emcee.utils.sample_ball(
        #     pos[idx_best], pos.std(axis=0), size=nwalkers)

        best_sol = [pos[idx_best], -prob[idx_best]]
        # Small ball (N sigma) around the best solution.
        mean_p, std_p = pos.mean(axis=0), pos.std(axis=0)
        low, high = mean_p - 1. * std_p, mean_p + 1. * std_p

        starting_guesses = []
        for p in pos:

            # Check that all elements in 'p' are between the (low,high)
            # values.
            flag = True
            for i, v in enumerate(p):
                # If any element is not inside its range, break out.
                if v <= low[i] or v >= high[i]:
                    flag = False
                    break

            if flag is True:
                starting_guesses.append(p)
            else:
                # If at least one element was not inside its range, store the
                # fill_value
                starting_guesses.append(best_sol[0])

        starting_guesses = np.asarray(starting_guesses)

    # Store burn-in chain phase.
    pars_chains_bi = sampler.chain.T[:, -nburn:, :]

    # Discard -np.inf posterior values.
    idx_keep = prob > -np.inf
    # Store MAP solution.
    idx_best = np.argmin(-prob[idx_keep])
    best_sol = [pos[idx_best], -prob[idx_best]]

    return pos, prob, state, best_sol, pars_chains_bi


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

        if priors == 'unif':
            # Flat prior
            lp = 0.

        elif priors == 'gauss':
            # TODO fix this
            model_mean = closeSol(
                fundam_params, varIdxs,
                [np.mean(fundam_params[_]) for _ in varIdxs])
            # Gaussian prior.
            model_std = np.array([0.005, .3, .1, .2, 500., 0.3])
            lp = np.sum(-np.square(
                (np.asarray(model_proper) - model_mean) / model_std))

    return lp, model_proper


def closeSol(fundam_params, varIdxs, model):
    """
    Find the closest value in the parameters list for metallicity, age, and
    mass.
    """
    model_proper, j = [], 0
    for i, par in enumerate(fundam_params):
        # If this parameter is one of the 'free' parameters.
        if i in varIdxs:
            # If it is the parameter metallicity, age or mass.
            if i in [0, 1, 4]:
                # Select the closest value in the array of allowed values.
                model_proper.append(
                    min(par, key=lambda x: abs(x - model[i - j])))
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


def convergenceVals(ndim, chains_nruns, emcee_trace):
    """
    Several convergence statistics.
    """
    # Autocorrelation time for each parameter.
    try:
        acorr_t = integrated_time(chains_nruns)
    except FloatingPointError:
        acorr_t = np.array([np.nan] * ndim)

    # Autocorrelation time for each chain for each parameter.
    at = []
    for p in chains_nruns.T:
        at_p = []
        for c in p.T:
            try:
                at_p.append(integrated_time(c, quiet=True)[0])
            except FloatingPointError:
                at_p.append(np.nan)
        at.append(at_p)
    # Select the indexes of the 10 chains with the largest acorr times, for
    # each parameter.
    max_at_10c = [np.argpartition(a, -10)[-10:] for a in at]

    # Worst chain: chain with the largest acor time.
    max_at_c = [np.argmax(a) for a in at]

    # Geweke z-scores and autocorrelation functions for the *worst* chain.
    geweke_z, emcee_acorf = [], []
    for i, p in enumerate(chains_nruns.T):
        geweke_z.append(geweke(p.T[max_at_c[i]]))
        try:
            emcee_acorf.append(emcee.autocorr.function(p.T[max_at_c[i]]))
        except FloatingPointError:
            emcee_acorf.append([np.nan])

    # PyMC3 effective sample size.
    try:
        pymc3_ess = effective_n(chains_nruns)
    except FloatingPointError:
        pymc3_ess = np.array([np.nan] * ndim)

    # Minimum effective sample size (ESS), and multi-variable ESS.
    minESS, mESS = fminESS(ndim), multiESS(emcee_trace.T)
    mESS_epsilon = [[], [], []]
    for alpha in [.01, .05, .1, .2, .3, .4, .5, .6, .7, .8, .9, .95]:
        mESS_epsilon[0].append(alpha)
        mESS_epsilon[1].append(fminESS(ndim, alpha=alpha, ess=minESS))
        mESS_epsilon[2].append(fminESS(ndim, alpha=alpha, ess=mESS))

    return acorr_t, max_at_10c, geweke_z, emcee_acorf, pymc3_ess, minESS,\
        mESS, mESS_epsilon

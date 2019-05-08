
import numpy as np
import time as t
import warnings

from .emcee3rc2 import ensemble
from .emcee3rc2 import moves
from .emcee3rc2 import utils
from .mcmc_common import interpSol

from ..synth_clust import synth_cluster
from . import likelihood
from .mcmc_convergence import convergenceVals
from .mcmc_common import random_population, varPars, closeSol, rangeCheck
from .. import update_progress


def main(
    lkl_method, e_max, err_lst, completeness, max_mag_syn, fundam_params,
        obs_clust, theor_tracks, R_V, ext_coefs, st_dist_mass, N_fc, cmpl_rnd,
        err_rnd, nwalkers, nsteps, nburn, N_burn, emcee_a, priors, hmax,
        *args):
    """

    nwalkers: number of MCMC walkers
    nwalkers: number of MCMC steps to take
    nburn: "burn-in" period to let chains stabilize

    """

    varIdxs, ndim, ranges = varPars(fundam_params)
    # Pack synthetic cluster arguments.
    synthcl_args = [
        theor_tracks, e_max, err_lst, completeness, max_mag_syn, st_dist_mass,
        R_V, ext_coefs, N_fc, cmpl_rnd, err_rnd]

    # TODO make this a proper parameter
    if emcee_a <= 0.:
        print("KDE")
        mv = moves.KDEMove()
    elif 0 < emcee_a <= 1.:
        print("StretchMove")
        mv = moves.StretchMove()
    elif 1 < emcee_a <= 2.:
        sigma = .1
        print("DE", sigma)
        mv = moves.DEMove(sigma)  # , gamma0=emcee_a)
    elif 2 < emcee_a <= 3.:
        gammas = 1.7
        print("DESnooker", gammas)
        mv = moves.DESnookerMove(gammas=gammas)
    elif 3 < emcee_a <= 4.:
        cov = .1
        print("Metropolis-Hastings", cov)
        mv = moves.GaussianMove(cov)
    elif 4 < emcee_a <= 5.:
        print("KDE + DESnooker")
        mv = [(moves.KDEMove(), 0.5), (moves.DESnookerMove(), 0.5)]
    elif 5 < emcee_a <= 6.:
        sigma = .05
        print("KDE + DE", sigma)
        mv = [(moves.KDEMove(), 0.5), (moves.DEMove(sigma), 0.5)]

    # TODO add this parameter to the input params file
    N_conv, tol_conv = 500., 0.01

    # Define sampler.
    sampler = ensemble.EnsembleSampler(
        nwalkers, ndim, log_posterior,
        args=[priors, varIdxs, ranges, fundam_params, synthcl_args, lkl_method,
              obs_clust], moves=mv)

    # Burn-in period.
    t0 = t.time()
    pos, map_sol_old, pars_chains_bi = burnIn(
        nwalkers, nburn, N_burn, fundam_params, varIdxs, sampler,
        ranges, priors, synthcl_args, lkl_method, obs_clust)

    elapsed = t.time() - t0
    max_secs = hmax * 60. * 60.
    available_secs = max(30, max_secs - elapsed)

    start_t = t.time()
    # We'll track how the average autocorrelation time estimate changes.
    # This will be useful to testing convergence.
    tau_index, autocorr_vals, old_tau = 0, np.empty(nsteps), np.inf

    # Check for convergence every 10% of steps or 100, whichever value
    # is lower.
    N_steps_conv = min(int(nsteps / 10.), 100)

    maf_steps, prob_mean, map_lkl = [], [], []
    milestones = list(range(10, 101, 10))
    for i, result in enumerate(sampler.sample(pos, iterations=nsteps)):

        elapsed += t.time() - start_t
        if elapsed >= available_secs:
            print("  Time consumed.")
            break
        start_t = t.time()

        # Only check convergence every 'N_steps_conv' steps
        if (i + 1) % N_steps_conv:
            continue

        # Compute the autocorrelation time so far. Using tol=0 means that
        # we'll always get an estimate even if it isn't trustworthy.
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                tau = sampler.get_autocorr_time(tol=0)
                # Check convergence
                converged = np.all(tau * N_conv < (i + 1))
                converged &= np.all(np.abs(old_tau - tau) / tau < tol_conv)
                autocorr_vals[tau_index] = np.nanmean(tau)
            tau_index += 1
            if converged:
                print("  Convergence achieved.")
                break
            old_tau = tau
        except FloatingPointError:
            pass

        pos, prob, state = result
        maf = np.mean(sampler.acceptance_fraction)
        maf_steps.append([i, maf])

        # Store MAP solution in this iteration.
        prob_mean.append([i, np.mean(prob)])
        idx_best = np.argmax(prob)
        # Update if a new optimal solution was found.
        if prob[idx_best] > map_sol_old[1]:
            map_sol_old = [
                closeSol(fundam_params, varIdxs, pos[idx_best]),
                prob[idx_best]]
        map_lkl.append([i, map_sol_old[1]])

        # Print progress.
        percentage_complete = (100. * (i + 1) / nsteps)
        if len(milestones) > 0 and percentage_complete >= milestones[0]:
            map_sol, logprob = map_sol_old
            m, s = divmod(nsteps / (i / elapsed) - elapsed, 60)
            h, m = divmod(m, 60)
            print("{:>3}% ({:.3f}) LP={:.1f} ({:g}, {:g}, {:.3f}, {:.2f}"
                  ", {:g}, {:.2f})".format(
                      milestones[0], maf, logprob, *map_sol) +
                  " [{:.0f} m/s | {:.0f}h{:.0f}m]".format(
                      (nwalkers * i) / elapsed, h, m))
            milestones = milestones[1:]

    runs = i + 1

    # Evolution of the mean autocorrelation time.
    tau_autocorr = autocorr_vals[:tau_index]

    # Final MAP fit.
    map_sol, map_lkl_final = map_sol_old

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
    # chains_nruns = discreteParams(fundam_params, varIdxs, chains_nruns)
    # Re-shape trace for all parameters (flat chain).
    # Shape: (ndim, runs * nwalkers)
    emcee_trace = chains_nruns.reshape(-1, ndim).T

    # Convergence parameters.
    acorr_t, max_at_c, min_at_c, geweke_z, emcee_acorf, mcmc_ess, minESS,\
        mESS, mESS_epsilon = convergenceVals(
            'emcee', ndim, varIdxs, N_conv, chains_nruns, emcee_trace)

    # Pass the mean as the best model fit found.
    mean_sol = closeSol(fundam_params, varIdxs, np.mean(emcee_trace, axis=1))

    isoch_fit_params = {
        'varIdxs': varIdxs, 'nsteps_emc': runs, 'mean_sol': mean_sol,
        'map_sol': map_sol, 'map_lkl': map_lkl, 'map_lkl_final': map_lkl_final,
        'prob_mean': prob_mean, 'mcmc_elapsed': elapsed,
        'mcmc_trace': emcee_trace,
        'pars_chains_bi': pars_chains_bi, 'pars_chains': chains_nruns.T,
        'maf_steps': maf_steps, 'autocorr_time': acorr_t,
        'max_at_c': max_at_c, 'min_at_c': min_at_c,
        'minESS': minESS, 'mESS': mESS, 'mESS_epsilon': mESS_epsilon,
        'emcee_acorf': emcee_acorf, 'geweke_z': geweke_z,
        'mcmc_ess': mcmc_ess,
        'N_steps_conv': N_steps_conv, 'N_conv': N_conv, 'tol_conv': tol_conv,
        'tau_index': tau_index, 'tau_autocorr': tau_autocorr
    }

    return isoch_fit_params


def burnIn(
    nwalkers, nburn, N_burn, fundam_params, varIdxs, sampler,
        ranges, priors, synthcl_args, lkl_method, obs_clust):
    """
    N_burn: number of times the burn-in phase will be applied, focusing on a
    small ball around the MAP.
    """
    # Random initial models.
    p0 = random_population(fundam_params, varIdxs, nwalkers)

    # # Initial population.
    # ntemps = 1
    # p0 = initPop(
    #     ranges, varIdxs, lkl_method, obs_clust, fundam_params, synthcl_args,
    #     ntemps, nwalkers, init_mode_ptm, popsize_ptm, maxiter_ptm)

    N_total, N_done, maf = nburn * N_burn, 0, np.nan
    print("     Burn-in stage")
    for _ in range(N_burn):

        while True:
            try:
                for i, result in enumerate(
                        sampler.sample(p0, iterations=nburn)):

                    if (i + 1) % int(nburn * .1):
                        maf = np.mean(sampler.acceptance_fraction)
                    update_progress.updt(
                        N_total, N_done + i + 1, "MAF={:.3f}".format(maf))
                pos, prob, state = result
                break
            except ValueError:
                print("  NaN/inf value found. Trying again.")
                pass

        N_done += nburn

        # Best solution (MAP) in this iteration.
        idx_best = np.argmax(prob)
        # TODO change when emcee3 is properly imported
        p0 = utils.sample_ball(
            pos[idx_best], pos.std(axis=0), size=nwalkers)

    # Store burn-in chain phase.
    chains_nruns = sampler.get_chain()[-nburn:, :, :]
    # pars_chains_bi = discreteParams(fundam_params, varIdxs, chains_nruns).T
    pars_chains_bi = chains_nruns.T

    # Store MAP solution.
    idx_best = np.argmax(prob)
    map_sol = [
        closeSol(fundam_params, varIdxs, pos[idx_best]), prob[idx_best]]

    # Reset the chain to remove the burn-in samples.
    sampler.reset()

    return pos, map_sol, pars_chains_bi


def log_posterior(
    model, priors, varIdxs, ranges, fundam_params, synthcl_args, lkl_method,
        obs_clust):
    """
    Log posterior function.
    """
    lp, isochrone, model_proper = log_prior(
        model, priors, fundam_params, varIdxs, ranges, synthcl_args[0])
    if not np.isfinite(lp):
        return -np.inf
    lkl = log_likelihood(
        isochrone, model_proper, fundam_params, synthcl_args, lkl_method,
        obs_clust)
    return lp + lkl


def log_prior(model, priors, fundam_params, varIdxs, ranges, theor_tracks):
    """
    met, age, ext, dm, mass, bf = model
    """
    rangeFlag = rangeCheck(model, ranges, varIdxs)

    lp, isochrone, model_proper = -np.inf, [], []
    # If some parameter is outside of the given ranges, don't bother obtaining
    # the proper model and just pass -inf
    if rangeFlag:
        # model_proper = closeSol(fundam_params, varIdxs, model)
        isochrone, model_proper = interpSol(
            theor_tracks, fundam_params, varIdxs, model)

        if priors == 'unif':
            # Flat prior
            lp = 0.

        elif priors == 'gauss':
            # TODO finish this
            model_mean = closeSol(
                fundam_params, varIdxs,
                [np.mean(fundam_params[_]) for _ in varIdxs])
            # Gaussian prior.
            model_std = np.array([0.005, .3, .001, .2, 500., 0.3])
            lp = np.sum(-np.square(
                (np.asarray(model_proper) - model_mean) / model_std))

    return lp, isochrone, model_proper


def log_likelihood(
        isochrone, model_proper, fundam_params, synthcl_args, lkl_method,
        obs_clust):
    """
    The Dolphin likelihood needs to be *minimized*. Be careful with the priors.
    """

    # TODO if 'm_sample' is True, always process the synthetic cluster and
    # the likelihood

    theor_tracks, e_max, err_lst, completeness, max_mag_syn, st_dist_mass,\
        R_V, ext_coefs, N_fc, cmpl_rnd, err_rnd = synthcl_args

    # # Metallicity and age indexes to identify isochrone.
    # m_i = fundam_params[0].index(model_proper[0])
    # a_i = fundam_params[1].index(model_proper[1])
    # isochrone = theor_tracks[m_i][a_i]

    # Generate synthetic cluster.
    synth_clust = synth_cluster.main(
        e_max, err_lst, completeness, max_mag_syn, st_dist_mass, isochrone,
        R_V, ext_coefs, N_fc, cmpl_rnd, err_rnd, model_proper)

    # Call likelihood function for this model. RETURNS THE INVERSE lkl.
    lkl = likelihood.main(lkl_method, synth_clust, obs_clust)

    # TODO, is this a general approach?
    # The negative likelihood is returned since Dolphin requires a minimization
    # of the PLR, and here we are maximizing
    return -lkl


import numpy as np
from scipy.optimize import differential_evolution as DE
import warnings
import time as t
from .. import update_progress
from ..synth_clust import synth_cluster
from . import likelihood
from .emcee_algor import varPars, random_population, closeSol, discreteParams
from .emcee3rc2 import autocorr
from .mcmc_convergence import convergenceVals
from .ptemcee import sampler
from .ptemcee import util


def main(
    lkl_method, e_max, err_lst, completeness, max_mag_syn,
    fundam_params, obs_clust, theor_tracks, R_V,
    ext_coefs, st_dist_mass, N_fc, cmpl_rnd, err_rnd, ntemps, nwalkers_ptm,
        nsteps_ptm, nburn_ptm, pt_adapt, tmax_ptm, priors_ptm):
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

    # TODO add these parameters to the input params file
    h_max = 20.
    max_secs = h_max * 60. * 60.

    ptsampler = sampler.Sampler(
        nwalkers_ptm, ndim, loglkl, logp,
        loglargs=[fundam_params, synthcl_args, lkl_method, obs_clust, ranges,
                  varIdxs, priors_ptm], Tmax=Tmax, ntemps=ntemps)

    # Start timing.
    available_secs = max(30, max_secs)
    elapsed, start_t = 0., t.time()

    # Initial population.
    p0 = initPop(
        ranges, varIdxs, lkl_method, obs_clust, fundam_params, synthcl_args,
        ntemps, nwalkers_ptm)

    print("     Burn-in stage")
    N_steps_check = max(1, int(nburn_ptm * .1))
    for i, (pos0, lnprob, lnlike) in enumerate(ptsampler.sample(
            p0, iterations=nburn_ptm, adapt=pt_adapt)):

        if (i + 1) % N_steps_check:
            continue
        maf = np.mean(ptsampler.acceptance_fraction[0])
        update_progress.updt(nburn_ptm, i + 1, "MAF={:.3f}".format(maf))

    # ptsampler.chain.shape: (ntemps, nwalkers, nsteps, ndim)
    # Store burn-in chain phase.
    # Shape: (runs, nwalkers, ndim)
    chains_nruns = ptsampler.chain[0].transpose(1, 0, 2)
    pars_chains_bi = discreteParams(fundam_params, varIdxs, chains_nruns).T

    # Store MAP solution.
    idx_best = np.argmax(lnprob[0])
    map_sol_old = [
        closeSol(fundam_params, varIdxs, pos0[0][idx_best]),
        lnprob[0][idx_best]]

    ptsampler.reset()

    # We'll track how the average autocorrelation time estimate changes.
    # This will be useful to testing convergence.
    tau_index, autocorr_vals, old_tau = 0, np.empty(nsteps_ptm), np.inf

    # Check for convergence every 5% of steps or 100, whichever value
    # is lower.
    N_steps_conv = min(int(nsteps_ptm * .1), 100)
    # TODO input as params
    N_conv, tol_conv = 1000., 0.01

    afs, tswaps, betas = [], [], []
    # actimes = []

    tau_emcee = np.empty(nsteps_ptm)

    # Check for dodgy inputs.
    if np.any(np.isinf(pos0)):
        print('At least one parameter value was infinite.')
        print(pos0)
    if np.any(np.isnan(pos0)):
        print('At least one parameter value was NaN.')
        print(pos0)

    maf_steps, prob_mean, map_lkl = [], [], []
    elapsed_in, start_in = 0., t.time()
    milestones = list(range(10, 101, 10))
    for i, result in enumerate(ptsampler.sample(
            pos0, iterations=nsteps_ptm, adapt=pt_adapt)):

        elapsed += t.time() - start_t
        if elapsed >= available_secs:
            print("  Time consumed (runs={})".format(i + 1))
            break
        start_t = t.time()

        # Only check convergence every 'N_steps_conv' steps
        if (i + 1) % N_steps_conv:
            continue

        # Compute the autocorrelation time so far. Using tol=0 means that
        # we'll always get an estimate even if it isn't trustworthy.
        tau0 = autocorr.integrated_time(
            ptsampler.chain[0, :, :i + 1, :].transpose(1, 0, 2), tol=0)
        tau_emcee[tau_index] = np.mean(tau0)

        # ptsampler.chain.shape: (ntemps, nwalkers, nsteps, ndim)
        x = np.mean(ptsampler.chain[0, :, :i + 1, :], axis=0)
        tau = util.autocorr_integrated_time(x)

        tswaps.append(ptsampler.tswap_acceptance_fraction)
        betas.append(ptsampler.betas)
        afs.append(np.mean(ptsampler.acceptance_fraction, axis=1))
        # acors = np.zeros(ntemps)
        # for temp in range(ntemps):
        #     x = np.mean(ptsampler.chain[temp, :, :i + 1, :], axis=0)
        #     acors[temp] = np.mean(util.autocorr_integrated_time(x))
        # actimes.append(acors)

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
        maf_steps.append([i, maf])

        # Store MAP solution in this iteration.
        prob_mean.append([i, np.mean(lnprob[0])])
        idx_best = np.argmax(lnprob[0])
        # Update if a new optimal solution was found.
        if lnprob[0][idx_best] > map_sol_old[1]:
            map_sol_old = [
                closeSol(fundam_params, varIdxs, pos[0][idx_best]),
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
            print("{:>3}% ({:.3f}) LP={:.1f} ({:g}, {:g}, {:.3f}, {:.2f}"
                  ", {:g}, {:.2f})".format(
                      milestones[0], maf, logprob, *map_sol) +
                  " [{:.0f} m/s | {:.0f}h{:.0f}m]".format(
                      (ntemps * nwalkers_ptm * i) / elapsed_in, h, m))
            milestones = milestones[1:]

    runs = i + 1

    # Evolution of the mean autocorrelation time.
    tau_autocorr = autocorr_vals[:tau_index]

    # Final MAP fit.
    map_sol, map_lkl_final = map_sol_old

    # This number should be between approximately 0.25 and 0.5 if everything
    # went as planned.
    m_accpt_fr = np.mean(ptsampler.acceptance_fraction[0])
    if m_accpt_fr > .5 or m_accpt_fr < .25:
        print("  WARNING: mean acceptance fraction is outside of the\n"
              "  recommended range.")

    # ptsampler.chain.shape: (ntemps, nwalkers, nsteps, ndim)
    # chains_nruns.shape: (runs, nwalkers, ndim)
    chains_nruns = ptsampler.chain[0, :, :runs, :].transpose(1, 0, 2)
    chains_nruns = discreteParams(fundam_params, varIdxs, chains_nruns)
    # Re-shape trace for all parameters (flat chain).
    # Shape: (ndim, runs * nwalkers)
    mcmc_trace = chains_nruns.reshape(-1, ndim).T

    # TODO delete
    # ptsampler.chain.shape: (ntemps, nwalkers, nsteps, ndim)
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    fig = plt.figure(figsize=(20, 20))
    gs = gridspec.GridSpec(10, 10)
    n = N_steps_conv * np.arange(1, tau_index + 1)

    ax = plt.subplot(gs[0:2, 0:6])
    ax.set_title(
        "ntemps: {}, Tmax: {}, adapt: {}".format(ntemps, Tmax, pt_adapt))
    for tsw in np.asarray(tswaps).T:
        plt.plot(n, tsw)
    plt.ylabel("Tswap AF")

    ax = plt.subplot(gs[2:4, 0:6])
    ax.set_title("Cold chain MAF: {:.5f}".format(
        np.mean(ptsampler.acceptance_fraction[0])))
    for af in np.asarray(afs).T:
        plt.plot(n, af)
    plt.ylabel("Mean AFs")

    ax = plt.subplot(gs[4:6, 0:6])
    ax.set_title("Beta_min: {:.5f}".format(ptsampler.betas[-1]))
    for bet in np.asarray(betas).T:
        plt.plot(n, bet)
    plt.ylabel("Betas")

    # plt.subplot(gs[6:8, 0:6])
    # for act in np.asarray(actimes).T:
    #     plt.plot(n, act)
    # plt.ylabel("ACTs")

    # plt.subplot(gs[8:10, 0:6])
    # y_min, y_max = np.inf, -np.inf
    # for i in range(nwalkers_ptm):
    #     chain = chains_nruns[:, i, 2]
    #     plt.plot(range(runs), chain)
    #     y_min, y_max = min(min(chain[-int(runs * .5):]), y_min),\
    #         max(max(chain[-int(runs * .5):]), y_max)
    # plt.ylim(y_min, y_max)

    plt.subplot(gs[6:8, 0:6])
    plt.plot(n, tau_emcee[:tau_index], ls=':', label="emcee")
    plt.plot(n, tau_autocorr, label="ptemcee")
    plt.plot(n, n / 100.0, "--k")
    plt.ylabel("ACT, cold chain")
    plt.legend()

    plt.xlabel("steps")
    fig.tight_layout()
    plt.savefig("ptemcee_{}b_{}s_{}w_{}t_{}_{}.png".format(
        nburn_ptm, runs, nwalkers_ptm, ntemps, Tmax, pt_adapt), dpi=300,
        bbox_inches='tight')
    plt.close()
    # TODO delete

    # Convergence parameters.
    acorr_t, max_at_c, min_at_c, geweke_z, emcee_acorf, mcmc_ess, minESS,\
        mESS, mESS_epsilon = convergenceVals(
            'ptemcee', ndim, varIdxs, N_conv, chains_nruns, mcmc_trace)

    # Pass the mean as the best model fit found.
    mean_sol = closeSol(fundam_params, varIdxs, np.mean(mcmc_trace, axis=1))

    isoch_fit_params = {
        'varIdxs': varIdxs, 'nsteps_ptm': runs, 'mean_sol': mean_sol,
        'map_sol': map_sol, 'map_lkl': map_lkl, 'map_lkl_final': map_lkl_final,
        'prob_mean': prob_mean,
        'mcmc_elapsed': elapsed, 'mcmc_trace': mcmc_trace,
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


def loglkl(
        model, fundam_params, synthcl_args, lkl_method, obs_clust, ranges,
        varIdxs, priors_ptm):
    """
    """
    rangeFlag = rangeCheck(model, ranges, varIdxs)

    logpost = -10000
    if rangeFlag:
        # Generate synthetic cluster.
        synth_clust = synthClust(fundam_params, varIdxs, model, synthcl_args)
        # Call likelihood function for this model.
        lkl = likelihood.main(lkl_method, synth_clust, obs_clust)

        logp = 0.
        # Logarithm of the prior.
        if priors_ptm == 'unif':
            # Flat prior
            logp = 0.

        # The negative likelihood is returned since Dolphin requires a
        # minimization of the PLR, and here we are maximizing
        logpost = logp + (-lkl)

    return logpost


def logp(_):
    # Just here as a place holder for 'ptemcee', the prior is inside the
    # log-likelihood.
    return 0.


def initPop(
    ranges, varIdxs, lkl_method, obs_clust, fundam_params, synthcl_args,
        ntemps, nwalkers_ptm, init_mode='diffevol', popsize=5, maxiter=10):
    """
    Obtain initial parameter values using either a random distribution, or
    the Differential Evolution algorithm to approximate reasonable solutions.
    """

    p0 = []
    if init_mode == 'random':
        print("     Rnd init pop")
        for _ in range(ntemps):
            p0.append(random_population(fundam_params, varIdxs, nwalkers_ptm))

    elif init_mode == 'diffevol':
        print("     DE init pop")

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # Estimate initial threshold value using DE.
            def dist(synth_clust, obs_clust):
                if synth_clust:
                    return likelihood.main(lkl_method, synth_clust, obs_clust)
                return np.inf

            def postfn(model):
                if rangeCheck(model, ranges, varIdxs):
                    return synthClust(
                        fundam_params, varIdxs, model, synthcl_args)
                return []

            def DEdist(model):
                synth_clust = postfn(model)
                return dist(synth_clust, obs_clust)

            runs = 1
            for _ in range(ntemps):
                walkers_sols = []
                for _ in range(nwalkers_ptm):
                    result = DE(
                        DEdist, ranges, popsize=popsize, maxiter=maxiter)
                    walkers_sols.append(result.x)
                    update_progress.updt(ntemps * nwalkers_ptm, runs)
                    runs += 1
                p0.append(walkers_sols)

        # Remove fixed parameters.
        p0arr = np.array(p0)
        p0 = []
        for ip in varIdxs:
            p0.append(p0arr[:, :, ip])
        p0 = np.array(p0).transpose(1, 2, 0)

    return p0


def rangeCheck(model, ranges, varIdxs):
    """
    Check that all the model values are within the given ranges.
    """
    check_ranges = [
        r[0] <= p <= r[1] for p, r in zip(*[model, ranges[varIdxs]])]
    if all(check_ranges):
        return True
    return False


def synthClust(fundam_params, varIdxs, model, synthcl_args):
    """
    Generate synthetic cluster.
    """
    model_proper = closeSol(fundam_params, varIdxs, model)

    theor_tracks, e_max, err_lst, completeness, max_mag_syn, st_dist_mass,\
        R_V, ext_coefs, N_fc, cmpl_rnd, err_rnd = synthcl_args

    # Metallicity and age indexes to identify isochrone.
    m_i = fundam_params[0].index(model_proper[0])
    a_i = fundam_params[1].index(model_proper[1])
    isochrone = theor_tracks[m_i][a_i]

    # Generate synthetic cluster.
    return synth_cluster.main(
        e_max, err_lst, completeness, max_mag_syn, st_dist_mass, isochrone,
        R_V, ext_coefs, N_fc, cmpl_rnd, err_rnd, model_proper)


import numpy as np
import time as t
from .. import update_progress
from ..synth_clust import synth_cluster
from . import likelihood
from .emcee_algor import varPars, random_population, convergenceVals,\
    closeSol, discreteParams, autocorr
from .ptemcee import sampler


def main(
    lkl_method, e_max, err_lst, completeness, max_mag_syn,
    fundam_params, obs_clust, theor_tracks, R_V,
    ext_coefs, st_dist_mass, N_fc, cmpl_rnd, err_rnd, ntemps, nwalkers_ptm,
        nsteps_ptm, nburn_ptm, ptemcee_a, priors_ptm):
    """
    """

    varIdxs, ndim, ranges = varPars(fundam_params)
    # Pack synthetic cluster arguments.
    synthcl_args = [
        theor_tracks, e_max, err_lst, completeness, max_mag_syn, st_dist_mass,
        R_V, ext_coefs, N_fc, cmpl_rnd, err_rnd]

    # TODO add these parameters to the input params file
    max_secs = 22. * 60. * 60.

    ptsampler = sampler.Sampler(
        nwalkers_ptm, ndim, logl, logp,
        loglargs=[fundam_params, synthcl_args, lkl_method, obs_clust, ranges,
                  varIdxs],
        logpargs=[priors_ptm, fundam_params, varIdxs, ranges], ntemps=ntemps)

    # Initial parameter values
    p0 = []
    for _ in range(ntemps):
        p0.append(random_population(fundam_params, varIdxs, nwalkers_ptm))
    # Shape p0: (ntemps, nwalkers_ptm, ndim)
    print("     Burn-in stage")
    t0 = t.time()
    N_steps_check = max(1, int(nburn_ptm * .1))
    for i, (pos0, lnprob, lnlike) in enumerate(ptsampler.sample(
            p0, iterations=nburn_ptm, adapt=True)):

        if (i + 1) % N_steps_check:
            continue
        maf = np.mean(ptsampler.acceptance_fraction[0])
        update_progress.updt(nburn_ptm, i + 1, "MAF={:.3f}".format(maf))

    # ptsampler.chain.shape: (ntemps, nwalkers, nsteps, ndim)

    # import matplotlib.pyplot as plt
    # for i, temp in enumerate(ptsampler.chain):
    #     print("Temp: {}".format(i))
    #     # First walker (ndim, nsteps)
    #     t_w = temp[0].T
    #     for j, t_w_dim in enumerate(t_w):
    #         ax = plt.subplot(int('33' + str(j + 1)))
    #         ax.hist(t_w_dim, bins=25)
    #     plt.show()

    # Store burn-in chain phase.
    chains_nruns = ptsampler.chain[0].reshape(nburn_ptm, nwalkers_ptm, ndim)
    pars_chains_bi = discreteParams(fundam_params, varIdxs, chains_nruns).T

    # Store MAP solution.
    idx_best = np.argmax(lnprob[0])
    best_sol_old = [
        closeSol(fundam_params, varIdxs, pos0[0][idx_best]),
        lnprob[0][idx_best]]

    ptsampler.reset()

    # Start timing.
    available_secs = max(30, max_secs)
    elapsed = t.time() - t0
    start_t = t.time()

    # We'll track how the average autocorrelation time estimate changes.
    # This will be useful to testing convergence.
    tau_index, autocorr_vals, old_tau = 0, np.empty(nsteps_ptm), np.inf

    # Check for convergence every 10% of steps or 100, whichever value
    # is lower.
    N_steps_conv = min(int(nsteps_ptm / 10.), 100)
    # TODO emcee3 input as params
    N_conv, tol_conv = 100., 0.01

    maf_steps, prob_mean, map_lkl = [], [], []
    milestones = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    for i, result in enumerate(ptsampler.sample(
            pos0, iterations=nsteps_ptm, adapt=True)):

        # Only check convergence every 'N_steps_conv' steps
        if (i + 1) % N_steps_conv:
            continue

        # Compute the autocorrelation time so far. Using tol=0 means that
        # we'll always get an estimate even if it isn't trustworthy.
        try:
            tau0 = ptsampler.get_autocorr_time()[0]
            tau = autocorr.integrated_time(
                ptsampler.chain[0, :, :i + 1, :].reshape(
                    i + 1, nwalkers_ptm, ndim), tol=0)
            print(tau0)
            print(tau)
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

        pos, lnprob, lnlike = result

        maf = np.mean(ptsampler.acceptance_fraction[0])
        maf_steps.append([i, maf])

        # Store MAP solution in this iteration.
        prob_mean.append([i, np.mean(lnprob[0])])
        idx_best = np.argmax(lnprob[0])
        # Update if a new optimal solution was found.
        if lnprob[0][idx_best] > best_sol_old[1]:
            best_sol_old = [
                closeSol(fundam_params, varIdxs, pos[0][idx_best]),
                lnprob[0][idx_best]]
        map_lkl.append([i, best_sol_old[1]])

        # Print progress.
        percentage_complete = (100. * (i + 1) / nsteps_ptm)
        if len(milestones) > 0 and percentage_complete >= milestones[0]:
            map_sol, logprob = best_sol_old
            m, s = divmod(nsteps_ptm / (i / elapsed) - elapsed, 60)
            h, m = divmod(m, 60)
            print("{:>3}% ({:.3f}) LP={:.1f} ({:g}, {:g}, {:.3f}, {:.2f}"
                  ", {:g}, {:.2f})".format(
                      milestones[0], maf, logprob, *map_sol) +
                  " [{:.0f} m/s | {:.0f}h{:.0f}m]".format(
                      (ntemps * nwalkers_ptm * i) / elapsed, h, m))
            milestones = milestones[1:]

        elapsed += t.time() - start_t
        if elapsed >= available_secs:
            print("  Time consumed.")
            break
        start_t = t.time()
    runs = i + 1

    # Evolution of the mean autocorrelation time.
    tau_autocorr = autocorr_vals[:tau_index]

    # Final MAP fit.
    map_sol, map_lkl_final = best_sol_old

    # This number should be between approximately 0.25 and 0.5 if everything
    # went as planned.
    m_accpt_fr = np.mean(ptsampler.acceptance_fraction[0])
    if m_accpt_fr > .5 or m_accpt_fr < .25:
        print("  WARNING: mean acceptance fraction is outside of the\n"
              "  recommended range.")

    # Shape: (runs, nwalkers, ndim)
    chains_nruns = ptsampler.chain[0].reshape(runs, nwalkers_ptm, ndim)
    chains_nruns = discreteParams(fundam_params, varIdxs, chains_nruns)
    # Re-shape trace for all parameters (flat chain).
    # Shape: (ndim, runs * nwalkers)
    mcmc_trace = chains_nruns.reshape(-1, ndim).T

    # Convergence parameters.
    acorr_t, max_at_5c, min_at_5c, geweke_z, emcee_acorf, pymc3_ess, minESS,\
        mESS, mESS_epsilon = convergenceVals(
            ndim, varIdxs, N_conv, chains_nruns, mcmc_trace)

    # Pass the mean as the best model fit found.
    best_sol = closeSol(fundam_params, varIdxs, np.mean(mcmc_trace, axis=1))

    isoch_fit_params = {
        'varIdxs': varIdxs, 'nsteps_ptm': runs, 'best_sol': best_sol,
        'map_sol': map_sol, 'map_lkl': map_lkl, 'map_lkl_final': map_lkl_final,
        'prob_mean': prob_mean,
        'mcmc_elapsed': elapsed, 'mcmc_trace': mcmc_trace,
        'pars_chains_bi': pars_chains_bi, 'pars_chains': chains_nruns.T,
        'maf_steps': maf_steps, 'autocorr_time': acorr_t,
        'max_at_5c': max_at_5c, 'min_at_5c': min_at_5c,
        'minESS': minESS, 'mESS': mESS, 'mESS_epsilon': mESS_epsilon,
        'emcee_acorf': emcee_acorf, 'geweke_z': geweke_z,
        'pymc3_ess': pymc3_ess,
        'N_steps_conv': N_steps_conv, 'N_conv': N_conv, 'tol_conv': tol_conv,
        'tau_index': tau_index, 'tau_autocorr': tau_autocorr
    }

    return isoch_fit_params


def logl(
        model, fundam_params, synthcl_args, lkl_method, obs_clust, ranges,
        varIdxs):
    """
    """

    check_ranges = [
        r[0] <= p <= r[1] for p, r in zip(*[model, ranges[varIdxs]])]

    lkl = -6000
    if all(check_ranges):
        model_proper = closeSol(fundam_params, varIdxs, model)

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

        # Call likelihood function for this model.
        lkl = likelihood.main(lkl_method, synth_clust, obs_clust)

        # Normalization  # TODO
        norm_lkl = 1.
        lkl = -lkl * norm_lkl

    # TODO, is this a general approach?
    # The negative likelihood is returned since Dolphin requires a minimization
    # of the PLR, and here we are maximizing
    return lkl


def logp(model, priors_emc, fundam_params, varIdxs, ranges):
    """
    """
    check_ranges = [
        r[0] <= p <= r[1] for p, r in zip(*[model, ranges[varIdxs]])]

    lp = -1e6
    # If some parameter is outside of the given ranges, don't bother obtaining
    # the proper model and just pass -inf
    if all(check_ranges):
        if priors_emc == 'unif':
            # Flat prior
            lp = 0.
    return lp

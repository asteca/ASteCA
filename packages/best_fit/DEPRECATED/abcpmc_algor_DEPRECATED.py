
import numpy as np
from scipy.optimize import differential_evolution as DE
import time as t
from .abcpmc import sampler, threshold
from ..synth_clust import synth_cluster
from . import likelihood
from .emcee_algor import varPars, closeSol, discreteParams, convergenceVals


def main(
    lkl_method, e_max, err_lst, completeness, max_mag_syn,
        fundam_params, obs_clust, theor_tracks, R_V, ext_coefs, st_dist_mass,
        N_fc, cmpl_rnd, err_rnd, nwalkers_abc, nsteps_abc, nburn_abc,
        priors_abc):

    varIdxs, ndim, ranges = varPars(fundam_params)

    def dist(synth_clust, obs_clust):
        lkl = np.inf
        if synth_clust:
            lkl = likelihood.main(lkl_method, synth_clust, obs_clust)
        return lkl

    def postfn(model):
        # Re-scale z and M
        model_scale = [
            model[0] / 100., model[1], model[2], model[3] * 10.,
            model[4] * 1000., model[5]]

        check_ranges = [
            r[0] <= p <= r[1] for p, r in zip(*[model_scale, ranges[varIdxs]])]

        synth_clust = []
        # If some parameter is outside of the given ranges, don't bother
        # obtaining the proper model.
        if all(check_ranges):
            model_proper = closeSol(fundam_params, varIdxs, model_scale)
            # Metallicity and age indexes to identify isochrone.
            m_i = fundam_params[0].index(model_proper[0])
            a_i = fundam_params[1].index(model_proper[1])
            isochrone = theor_tracks[m_i][a_i]

            # Generate synthetic cluster.
            synth_clust = synth_cluster.main(
                e_max, err_lst, completeness, max_mag_syn, st_dist_mass,
                isochrone, R_V, ext_coefs, N_fc, cmpl_rnd, err_rnd,
                model_proper)

        return synth_clust

    # TODO add these parameters to the input params file
    alpha, init_eps = 95, None
    N_conv, tol_conv = 50., 0.01
    max_secs = 22. * 60. * 60.
    # Break out when AF is low.
    # af_low, af_min_steps = 0.001, .1
    max_t_walker = 30.
    # eps_stuck_perc, N_eps_stuck_max = .005, 100

    # Start timing.
    elapsed = 0.
    available_secs = max(30, max_secs)
    start_t = t.time()

    abcsampler = sampler.Sampler(
        N=nwalkers_abc, Y=obs_clust, postfn=postfn, dist=dist)
    # Set proposal
    # sampler.particle_proposal_cls = sampler.OLCMParticleProposal

    if init_eps is None:
        # Estimate initial threshold value using DE.
        def lnprob(model):
            synth_clust = postfn(model)
            return dist(synth_clust, obs_clust)
        # Scale parameters bounds.
        bounds = [
            ranges[0] * 100., ranges[1], ranges[2], ranges[3] / 10.,
            ranges[4] / 1000., ranges[5]]
        result = DE(lnprob, bounds, maxiter=20)
        init_eps = 4. * result.fun
    print(" Initial threshold value: {:.2f}".format(init_eps))

    # old_eps = init_eps
    # TODO pass type of threshold from params file
    # eps = threshold.LinearEps(T, 5000, init_eps)
    eps = threshold.ConstEps(nsteps_abc, init_eps)

    # Stddev values as full range.
    std = np.eye(ndim) * (ranges.max(axis=1) - ranges.min(axis=1))
    # Means as middle points in ranges.
    means = (ranges.max(axis=1) + ranges.min(axis=1)) / 2.
    # Scale values.
    std[0], means[0] = std[0] * 100, means[0] * 100
    std[3], means[3] = std[3] / 10, means[3] / 10
    std[4], means[4] = std[4] / 1000., means[4] / 1000.
    # Gaussian prior.
    print(means)
    print(std)
    prior = sampler.GaussianPrior(mu=means, sigma=std)

    # # We'll track how the average autocorrelation time estimate changes
    # tau_index, autocorr_vals = 0, np.empty(nsteps_abc)
    # # This will be useful to testing convergence
    # old_tau = np.inf

    # Check for convergence every 2% of steps or 100, whichever value
    # is lower.
    # N_steps_conv = min(int(nsteps_abc * 0.02), 100)

    map_sol_old, N_models, prob_mean = [[], np.inf], 0, []
    # N_eps_stuck = 0
    chains_nruns, maf_steps, map_lkl = [], [], []
    milestones = list(range(5, 101, 5))
    for pool in abcsampler.sample(prior, eps):

        print(
            pool.t, pool.eps, pool.ratio, np.min(pool.dists),
            np.mean(pool.dists))

        chains_nruns.append(pool.thetas)
        maf = pool.ratio
        maf_steps.append([pool.t, maf])
        N_models += nwalkers_abc / maf

        # reduce eps value
        # old_eps = eps.eps
        eps.eps = np.percentile(pool.dists, alpha)

        # # Check if threshold is stuck.
        # if abs(eps.eps - old_eps) < eps_stuck_perc * eps.eps:
        #     N_eps_stuck += 1
        # else:
        #     N_eps_stuck = 0
        # if N_eps_stuck > N_eps_stuck_max:
        #     print("  Threshold is stuck (runs={}).".format(pool.t + 1))
        #     break

        # if maf < af_low and pool.t > int(af_min_steps * nsteps_abc):
        #     print("  AF<{} (runs={})".format(af_low, pool.t + 1))
        #     break

        if t.time() - start_t > (max_t_walker * nwalkers_abc):
            print("  Sampler is stuck (runs={})".format(pool.t + 1))
            break

        elapsed += t.time() - start_t
        if elapsed >= available_secs:
            print("  Time consumed (runs={})".format(pool.t + 1))
            break
        start_t = t.time()

        # # Only check convergence every 'N_steps_conv' steps
        # if (pool.t + 1) % N_steps_conv:
        #     continue

        # # Compute the autocorrelation time so far. Using tol=0 means that
        # # we'll always get an estimate even if it isn't trustworthy.
        # try:
        #     tau = autocorr.integrated_time(np.array(chains_nruns), tol=0)
        #     autocorr_vals[tau_index] = np.nanmean(tau)
        #     tau_index += 1

        #     # Check convergence
        #     converged = np.all(tau * N_conv < (pool.t + 1))
        #     converged &= np.all(np.abs(old_tau - tau) / tau < tol_conv)
        #     if converged:
        #         print("  Convergence achieved (runs={}).".format(pool.t + 1))
        #         break
        #     old_tau = tau
        # except FloatingPointError:
        #     pass

        # Store MAP solution in this iteration.
        prob_mean.append([pool.t, np.mean(pool.dists)])
        idx_best = np.argmin(pool.dists)
        # Update if a new optimal solution was found.
        if pool.dists[idx_best] < map_sol_old[1]:
            pars = pool.thetas[idx_best]
            # pars = scaleParams(model)
            pars = [pars[0] / 100., pars[1], pars[2], pars[3] * 10.,
                    pars[4] * 1000., pars[5]]
            map_sol_old = [
                closeSol(fundam_params, varIdxs, pars),
                pool.dists[idx_best]]
        map_lkl.append([pool.t, map_sol_old[1]])

        # Print progress.
        percentage_complete = (100. * (pool.t + 1) / nsteps_abc)
        if len(milestones) > 0 and percentage_complete >= milestones[0]:
            map_sol, logprob = map_sol_old
            print("{:>3}% ({:.3f}) LP={:.1f} ({:g}, {:g}, {:.3f}, {:.2f}"
                  ", {:g}, {:.2f})".format(
                      milestones[0], maf, logprob, *map_sol) +
                  " [{:.0f} m/s]".format(N_models / elapsed))
            milestones = milestones[1:]

    runs = pool.t + 1

    # Evolution of the mean autocorrelation time.
    tau_autocorr = np.array([np.nan] * 10)  # autocorr_vals[:tau_index]
    tau_index = np.nan
    N_steps_conv = runs

    # Final MAP fit.
    idx_best = np.argmin(pool.dists)
    pars = pool.thetas[idx_best]
    # pars = scaleParams(model)
    pars = [
        pars[0] / 100., pars[1], pars[2], pars[3] * 10., pars[4] * 1000.,
        pars[5]]
    map_sol = closeSol(fundam_params, varIdxs, pars)
    map_lkl_final = pool.dists[idx_best]

    abcsampler.close()

    # Shape: (runs, nwalkers, ndim)
    chains_nruns = np.array(chains_nruns)
    # De-scale parameters.
    chains_nruns[:, :, 0] = chains_nruns[:, :, 0] / 100.
    chains_nruns[:, :, 3] = chains_nruns[:, :, 3] * 10.
    chains_nruns[:, :, 4] = chains_nruns[:, :, 4] * 1000.

    # Burn-in range.
    Nb = int(runs * nburn_abc)

    # Burn-in. Shape: (ndim, nwalkers, runs)
    pars_chains_bi = discreteParams(
        fundam_params, varIdxs, chains_nruns[:Nb, :, :]).T

    # Change values for the discrete parameters with the closest valid values.
    chains_nruns = discreteParams(
        fundam_params, varIdxs, chains_nruns[Nb:, :, :])

    mcmc_trace = chains_nruns.reshape(-1, ndim).T

    # import matplotlib.pyplot as plt
    # import corner
    # corner.corner(
    #     mcmc_trace.T, quantiles=[0.16, 0.5, 0.84], show_titles=True)
    #     # levels=(1 - np.exp(-0.5),))
    # plt.savefig("corner.png", dpi=300)

    # Convergence parameters.
    acorr_t, max_at_c, min_at_c, geweke_z, emcee_acorf, mcmc_ess, minESS,\
        mESS, mESS_epsilon = convergenceVals(
            'abc', ndim, varIdxs, N_conv, chains_nruns, mcmc_trace)

    # Store mean solution.
    mean_sol = closeSol(fundam_params, varIdxs, np.mean(mcmc_trace, axis=1))

    isoch_fit_params = {
        'varIdxs': varIdxs, 'nsteps_abc': runs, 'mean_sol': mean_sol,
        'nburn_abc': Nb, 'map_sol': map_sol, 'map_lkl': map_lkl,
        'map_lkl_final': map_lkl_final, 'prob_mean': prob_mean,
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

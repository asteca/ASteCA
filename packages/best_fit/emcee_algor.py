
import numpy as np
import random
import emcee
from ..synth_clust import synth_cluster
import likelihood
from .. import update_progress


def varPars(fundam_params):
    """
    """
    ranges = [[min(_), max(_)] for _ in fundam_params]

    varIdxs = []
    for i, r in enumerate(ranges):
        if r[0] != r[1]:
            varIdxs.append(i)

    # Number of free parameters.
    ndim = len(varIdxs)

    return varIdxs, ndim, np.array(ranges)


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

    return model_proper


def log_prior(model, fundam_params, varIdxs, ranges):
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

        # # Gaussian prior.
        # mu_std = [[0.0155, 8., .3, 13., 2000., 0.],
        #           [0.005, .5, .5, .5, 1000., 0.01]]
        # pr_pars = -np.square(
        #     (np.asarray(model_proper) - mu_std[0]) / mu_std[1])
        # lp = np.sum(pr_pars)

        # Flat prior
        lp = 0.

    return lp, model_proper


def log_likelihood(model_proper, fundam_params, lkl_args):
    """
    """
    theor_tracks, lkl_method, obs_clust, e_max, err_lst, completeness,\
        max_mag_syn, st_dist_mass, R_V, ext_coefs, N_fc, cmpl_rnd, err_rnd =\
        lkl_args

    # Metallicity and age indexes to identify isochrone.
    m_i = fundam_params[0].index(model_proper[0])
    a_i = fundam_params[1].index(model_proper[1])
    isochrone = theor_tracks[m_i][a_i]

    # Generate synthetic cluster.
    synth_clust = synth_cluster.main(
        e_max, err_lst, completeness, max_mag_syn, st_dist_mass, isochrone,
        R_V, ext_coefs, N_fc, cmpl_rnd, err_rnd, model_proper)

    # Call likelihood function for this model.
    lkl = likelihood.main(
        lkl_method, synth_clust, obs_clust)

    return -lkl


def log_posterior(model, varIdxs, ranges, fundam_params, lkl_args):
    """
    """
    lp, model = log_prior(model, fundam_params, varIdxs, ranges)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(model, fundam_params, lkl_args)


def main(lkl_method, e_max, err_lst, completeness, max_mag_syn,
         fundam_params, obs_clust, theor_tracks, R_V, ext_coefs, st_dist_mass,
         N_fc, cmpl_rnd, err_rnd, nwalkers, nsteps, nburn, *args):
    """

    nwalkers: number of MCMC walkers
    nwalkers: number of MCMC steps to take
    nburn: "burn-in" period to let chains stabilize

    TODO: add error in check func:

    AssertionError: The number of walkers needs to be more than twice the
    dimension of your parameter space... unless you're crazy!

    AssertionError: The number of walkers must be even.
    """

    varIdxs, ndim, ranges = varPars(fundam_params)

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

    lkl_args = [
        theor_tracks, lkl_method, obs_clust, e_max, err_lst, completeness,
        max_mag_syn, st_dist_mass, R_V, ext_coefs, N_fc, cmpl_rnd, err_rnd]

    # from emcee import PTSampler
    # ntemps = 20

    # sampler = PTSampler(
    #     ntemps, nwalkers, ndim, logp=log_prior, logl=log_likelihood,
    #     logpargs=[fundam_params, varIdxs, ranges],
    #     loglargs=[fundam_params, varIdxs, lkl_args])

    # # Burn-in
    # p0 = []
    # for _ in range(ntemps):
    #     p0.append(random_population(fundam_params, varIdxs, nwalkers))
    # p0 = np.asarray(p0)
    # print(p0.shape, ntemps, nwalkers, ndim)
    # for p, lnprob, lnlike in sampler.sample(p0, iterations=nburn):
    #     pass
    # sampler.reset()

    # for p, lnprob, lnlike in sampler.sample(
    #         p, lnprob0=lnprob, lnlike0=lnlike, iterations=nsteps, thin=2):
    #     pass

    # Generate sampler.
    sampler = emcee.EnsembleSampler(
        nwalkers, ndim, log_posterior,
        args=[varIdxs, ranges, fundam_params, lkl_args])

    # Run burn-in.
    print("     Burn-in stage")
    for i, result in enumerate(
            sampler.sample(starting_guesses, iterations=nburn)):
        update_progress.updt(nburn, i + 1)
    pos, prob, state = result

    # pos, prob, state = sampler.run_mcmc(starting_guesses, nburn)

    # Store burn-in chain phase.
    pars_chains_bi = sampler.chain.T
    # Reset the chain to remove the burn-in samples.
    sampler.reset()
    # Discard -np.inf posterior values.
    idx_keep = prob > -np.inf
    # Store best solution in this iteration.
    idx_best = np.argmin(-prob[idx_keep])
    best_sol_old = [pos[idx_best], -prob[idx_best]]

    milestones = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    # Run 'nsteps'
    for i, result in enumerate(
            sampler.sample(
                pos, lnprob0=prob, rstate0=state, iterations=nsteps)):

        # Print progress.
        percentage_complete = (100.0 * (i + 1) / nsteps)
        if len(milestones) > 0 and percentage_complete >= milestones[0]:

            # Discard -np.inf posterior values.
            idx_keep = result[1] > -np.inf
            # Store MAP solution in this iteration.
            idx_best = np.argmin(-result[1][idx_keep])

            best_sol_run = [result[0][idx_best], -result[1][idx_best]]
            # Check if a new optimal solution was found.
            if best_sol_run[1] < best_sol_old[1]:
                best = closeSol(fundam_params, varIdxs, best_sol_run[0])
                lkl = best_sol_run[1]
                best_sol_old = [best_sol_run[0], best_sol_run[1]]
            else:
                best = closeSol(fundam_params, varIdxs, best_sol_old[0])
                lkl = best_sol_old[1]
            maf = np.mean(sampler.acceptance_fraction)
            print(" {:>3}% ({:.3f}) L={:.1f} ({:g}, {:g}, {:.3f}, {:.2f}"
                  ", {:g}, {:.2f})".format(milestones[0], maf, lkl, *best))
            milestones = milestones[1:]

    # TODO should I pass the MAP solution or the median?

    # Pass the median as the best model fit found.
    # best = closeSol(fundam_params, np.median(emcee_trace, axis=1))

    # Pass the best model fit found (MAP).
    idx = np.argmax(sampler.flatlnprobability)
    best = closeSol(fundam_params, varIdxs, sampler.flatchain[idx])
    print("MAP sol: {:g} {:g} {:.3f} {:.2f} {:g} {:.2f}".format(
        *best))

    # This number should be between approximately 0.25 and 0.5 if everything
    # went as planned.
    m_accpt_fr = np.mean(sampler.acceptance_fraction)
    print("Mean acceptance fraction: {:.3f}".format(m_accpt_fr))
    if m_accpt_fr > .5 or m_accpt_fr < .25:
        print("  WARNING: mean acceptance fraction is outside of the\n"
              "  recommended range.")

    # Estimate the integrated autocorrelation time for the time series in each
    # parameter.
    try:
        print("Autocorrelation time: {:.2f}".format(
            sampler.get_autocorr_time()))
    except Exception:
        print("  WARNING: the chain is too short to reliably estimate\n"
              "  the autocorrelation time.")

    # Throw-out the burn-in points and reshape.
    # emcee_trace = sampler.chain[:, nburn:, :].reshape(-1, ndim).T
    # Re-shape trace for all parameters.
    emcee_trace = sampler.chain.reshape(-1, ndim).T

    # Chains for all parameters, post burn-in.
    pars_chains = sampler.chain.T

    isoch_fit_params = [
        best, [pars_chains_bi, pars_chains], m_accpt_fr, varIdxs, emcee_trace]

    # # TODO un-comment to generate corner plot
    # import matplotlib.pyplot as plt
    # import corner
    # all_labels = [r"$z$", r"$\log (age)$", r"$E_{BV}$", r"$\mu_0$",
    #               r"$M_{\odot}$", r"$b_{frac}$"]
    # labs = []
    # # Only plot parameters with dynamical ranges.
    # for i, r in enumerate(ranges):
    #     if min(r) != max(r):
    #         labs.append(all_labels[i])
    #         # emcee_r.append(emcee_trace[i])
    # # plt.figure(figsize=(20, 20))
    # # Medians.
    # mdns = [np.percentile(p, 50) for p in emcee_trace]
    # corner.corner(
    #     np.array(emcee_trace).T, color='#4682b4', truths=mdns,
    #     truth_color='r', labels=labs, quantiles=[0.16, 0.84],
    #     # levels=(1. - np.exp(-0.5),),
    #     show_titles=True, title_kwargs={"fontsize": 11})
    # # plt.savefig(str(random.randint(0, 100000000)) + ".png", dpi=150)
    # # plt.clf()
    # # plt.close("all")
    # plt.show()

    # mx, my = 0, 1
    # import matplotlib.pyplot as plt
    # plt.figure(figsize=(10, 10))
    # ax = plt.subplot(111)
    # x = np.array(zip(*model_done[0])[mx])
    # y = np.array(zip(*model_done[0])[my])
    # h2d, xbins, ybins = plt.hist2d(
    #     x, y, bins=20, cmap=plt.get_cmap('Blues'), zorder=2)[:-1]
    # plt.contour(
    #     h2d.transpose(), levels=(1. - np.exp(-0.5),),
    #     extent=[xbins.min(), xbins.max(), ybins.min(), ybins.max()],
    #     colors='#551a8b', linewidths=2.5, zorder=3)
    # ax.set_aspect('auto')
    # plt.savefig("emcee2.png", dpi=150)

    return isoch_fit_params

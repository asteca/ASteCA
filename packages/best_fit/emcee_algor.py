
import numpy as np
import random
import emcee
from ..synth_clust import synth_cluster
import likelihood
from .. import update_progress

import time as t


def main(
    lkl_method, e_max, err_lst, completeness, max_mag_syn, fundam_params,
        obs_clust, theor_tracks, R_V, ext_coefs, st_dist_mass, N_fc, cmpl_rnd,
        err_rnd, nwalkers, nsteps, nburn, *args):
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

    synthcl_args = [
        theor_tracks, e_max, err_lst, completeness, max_mag_syn, st_dist_mass,
        R_V, ext_coefs, N_fc, cmpl_rnd, err_rnd]

    s = t.time()
    # Generate sampler.
    model_done = {}
    sampler = emcee.EnsembleSampler(
        nwalkers, ndim, log_posterior,
        args=[varIdxs, ranges, fundam_params, synthcl_args, lkl_method,
              obs_clust, model_done])

    # Burn-in period.
    pos, prob, state, best_sol_old, pars_chains_bi = burnIn(
        nwalkers, nburn, fundam_params, varIdxs, sampler)

    # Reset the chain to remove the burn-in samples.
    sampler.reset()

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
                logprob = best_sol_run[1]
                best_sol_old = [best_sol_run[0], best_sol_run[1]]
            else:
                best = closeSol(fundam_params, varIdxs, best_sol_old[0])
                logprob = best_sol_old[1]
            maf = np.mean(sampler.acceptance_fraction)
            print(" {:>3}% ({:.3f}) LP={:.1f} ({:g}, {:g}, {:.3f}, {:.2f}"
                  ", {:g}, {:.2f})".format(milestones[0], maf, logprob, *best))
            milestones = milestones[1:]

    # TODO should I pass the MAP solution or the median?

    # Pass the mean as the best model fit found.
    # best = closeSol(fundam_params, np.median(emcee_trace, axis=1))
    print("TIME: {:.3f}".format(t.time() - s))


    # Pass the best model fit found (MAP).
    idx = np.argmax(sampler.flatlnprobability)
    best = closeSol(fundam_params, varIdxs, sampler.flatchain[idx])
    print("  MAP: {:g} {:g} {:.3f} {:.2f} {:g} {:.2f}".format(
        *best))

    # This number should be between approximately 0.25 and 0.5 if everything
    # went as planned.
    m_accpt_fr = np.mean(sampler.acceptance_fraction)
    print("  Mean acceptance fraction: {:.3f}".format(m_accpt_fr))
    if m_accpt_fr > .5 or m_accpt_fr < .25:
        print("  WARNING: mean acceptance fraction is outside of the\n"
              "  recommended range.")

    # Re-shape trace for all parameters.
    emcee_trace = sampler.chain.reshape(-1, ndim).T

    # Estimate the integrated autocorrelation time for the time series in each
    # parameter.
    try:
        # print("  Autocorrelation time: {:.2f}".format(
        #     sampler.get_autocorr_time()))
        print("  Autocorrelation time: ", emcee.autocorr.integrated_time(
            emcee_trace))
    except Exception:
        print("  WARNING: the chain is too short to reliably estimate\n"
              "  the autocorrelation time.")

    print("  G-R: ", gelman_rubin(sampler.chain))
    import pdb; pdb.set_trace()  # breakpoint b512d290 //
    

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


def burnIn(nwalkers, nburn, fundam_params, varIdxs, sampler):
    """
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

    N_b = 5
    N_total, N_done = nburn * N_b, 0
    # Run burn-in.
    print("     Burn-in stage")
    for _ in range(N_b):

        for i, result in enumerate(
                sampler.sample(starting_guesses, iterations=nburn)):
            update_progress.updt(N_total, N_done + i + 1)
        N_done += nburn

        pos, prob, state = result

        # Discard -np.inf posterior values.
        idx_keep = prob > -np.inf
        # Store best solution in this iteration.
        idx_best = np.argmin(-prob[idx_keep])
        best_sol = [pos[idx_best], -prob[idx_best]]

        # Small ball (N sigma) around the best solution.
        mean_p, std_p = pos.mean(axis=0), pos.std(axis=0)
        low, high = mean_p - 3. * std_p, mean_p + 3. * std_p

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
    pars_chains_bi = sampler.chain.T[:, :nburn, :]
    # if _ == 0:
    #     pars_chains_bi = sampler.chain.T
    # else:
    #     pars_chains_bi = np.concatenate(
    #         (pars_chains_bi, sampler.chain.T), axis=1)

    # # Reset the chain to remove the burn-in samples.
    # sampler.reset()

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
    model, varIdxs, ranges, fundam_params, synthcl_args, lkl_method, obs_clust,
        model_done):
    """
    Log posterior function.
    """
    lp, model = log_prior(model, fundam_params, varIdxs, ranges)
    if not np.isfinite(lp):
        return -np.inf
    lkl, model_done = log_likelihood(
        model, fundam_params, synthcl_args, lkl_method, obs_clust, model_done)
    return lp + lkl


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
                # Store rounded parameter values to make the 'model_done'
                # key assignment more stable (else there could be a large
                # number of decimals)
                model_proper.append(model[i - j])
        else:
            model_proper.append(par[0])
            j += 1

    return np.round(model_proper, 5)


def log_likelihood(
    model_proper, fundam_params, synthcl_args, lkl_method, obs_clust,
        model_done):
    """
    """
    # TODO if 'm_sample' is True, always process the synthetic cluster and
    # the likelihood

    # If the model was already processed, use the generated likelihood value.
    try:
        lkl = model_done[''.join(map(str, model_proper))]
    except KeyError:

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

        # Store processed model and its likelihood value.
        model_done[''.join(map(str, model_proper))] = lkl

    return -lkl, model_done


def gelman_rubin(x, return_var=False):
    """
    Source PyMC: https://github.com/pymc-devs/pymc

    Returns estimate of R for a set of traces.
    The Gelman-Rubin diagnostic tests for lack of convergence by comparing
    the variance between multiple chains to the variance within each chain.
    If convergence has been achieved, the between-chain and within-chain
    variances should be identical. To be most effective in detecting evidence
    for nonconvergence, each chain should have been initialized to starting
    values that are dispersed relative to the target distribution.
    Parameters
    ----------
    x : array-like
      An array containing the 2 or more traces of a stochastic parameter.
      That is, an array of dimension m x n x k, where m is the number of
      traces, n the number of samples, and k the dimension of the stochastic.

    return_var : bool
      Flag for returning the marginal posterior variance instead of R-hat
      (defaults of False).
    Returns
    -------
    Rhat : float
      Return the potential scale reduction factor, :math:`\hat{R}`
    Notes
    -----
    The diagnostic is computed by:
      .. math:: \hat{R} = \sqrt{\frac{\hat{V}}{W}}
    where :math:`W` is the within-chain variance and :math:`\hat{V}` is
    the posterior variance estimate for the pooled traces. This is the
    potential scale reduction factor, which converges to unity when each
    of the traces is a sample from the target posterior. Values greater
    than one indicate that one or more chains have not yet converged.
    References
    ----------
    Brooks and Gelman (1998)
    Gelman and Rubin (1992)"""

    if np.shape(x) < (2,):
        print("Gelman-Rubin diagnostic requires multiple chains of the same"
              " length.")
        return np.nan

    try:
        m, n = np.shape(x)
    except ValueError:
        return [gelman_rubin(np.transpose(y)) for y in np.transpose(x)]

    # Calculate between-chain variance
    B_over_n = np.sum((np.mean(x, 1) - np.mean(x)) ** 2) / (m - 1)

    # Calculate within-chain variances
    W = np.sum(
        [(x[i] - xbar) ** 2 for i,
         xbar in enumerate(np.mean(x,
                                   1))]) / (m * (n - 1))

    # (over) estimate of variance
    s2 = W * (n - 1) / n + B_over_n

    if return_var:
        return s2

    # Pooled posterior variance estimate
    V = s2 + B_over_n / m

    # Calculate PSRF
    R = V / W

    return np.sqrt(R)

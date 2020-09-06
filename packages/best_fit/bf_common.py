
import numpy as np
from scipy.optimize import differential_evolution as DE
from scipy import stats
import warnings
from ..synth_clust import synth_cluster
from ..aux_funcs import kde1D, reject_outliers
from . import likelihood
from .. import update_progress


def varPars(fundam_params):
    """
    Check which parameters are fixed and which have a dynamic range. This also
    dictates the number of free parameters.
    """
    ranges = [np.array([min(_), max(_)]) for _ in fundam_params]

    varIdxs = []
    for i, r in enumerate(ranges):
        if r[0] != r[1]:
            varIdxs.append(i)

    # Number of free parameters.
    ndim = len(varIdxs)

    return varIdxs, ndim, np.array(ranges)


def fillParams(fundam_params, varIdxs, model):
    """
    Fills the places in 'model' of the parameters that were not fitted, with
    their fixed values.
    """

    model_filled, j = [], 0
    for i, par in enumerate(fundam_params):
        # If this parameter is one of the 'free' parameters, use its own value.
        if i in varIdxs:
            model_filled.append(model[i - j])
        else:
            # Fill with the fixed value for this parameter.
            model_filled.append(par[0])
            j += 1

    return model_filled


#  DEPRECATED 02-10-2019
# def closeSol(fundam_params, model, pushidxs):
#     """
#     Find the closest value in the parameters list for the discrete parameters
#     metallicity, age, and mass.
#     """
#     model_proper = []
#     for i, par in enumerate(fundam_params):
#         # If it is the parameter metallicity, age or mass.
#         if i in pushidxs:
#             # Select the closest value in the array of allowed values.
#             model_proper.append(min(par, key=lambda x: abs(x - model[i])))
#         else:
#             model_proper.append(model[i])

#     return model_proper


def initPop(
    ranges, varIdxs, lkl_method, obs_clust, fundam_params,
        synthcl_args, ntemps, nwalkers, init_mode, popsize, maxiter):
    """
    Obtain initial parameter values using either a random distribution, or
    the Differential Evolution algorithm to approximate reasonable solutions.
    """

    p0 = []
    if init_mode == 'random':
        # print("Random initial population")
        for _ in range(ntemps):
            p0.append(random_population(fundam_params, varIdxs, nwalkers))

    elif init_mode == 'diffevol':
        # DEPRECATED 05/09/19 when #425 was implemented
        print("     DE init pop")

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # Estimate initial threshold value using DE.
            def DEdist(model):
                synth_clust = synth_cluster.main(
                    fundam_params, varIdxs, model, *synthcl_args)
                if synth_clust:
                    lkl = likelihood.main(lkl_method, synth_clust, obs_clust)
                    return lkl
                return np.inf

            walkers_sols = []
            for _ in range(nwalkers):
                result = DE(
                    DEdist, ranges[varIdxs], popsize=popsize, maxiter=maxiter)
                walkers_sols.append(result.x)
                update_progress.updt(nwalkers, _ + 1)

        p0 = [walkers_sols for _ in range(ntemps)]

    return p0


def random_population(fundam_params, varIdxs, n_ran):
    """
    Generate a random set of parameter values to use as a random population.

    Pick n_ran initial random solutions from the valid ranges.
    """
    p_lst = []
    for i in varIdxs:
            p_lst.append(
                np.random.uniform(
                    fundam_params[i][0], fundam_params[i][-1], n_ran))

    return np.array(p_lst).T


#  DEPRECATED 24-09-2019
# def discreteParams(fundam_params, varIdxs, chains_nruns, pushidxs):
#     """
#     Push values in each chain for each discrete parameter in the 'pushidxs'
#     list to the closest grid value.

#     chains_nruns.shape: (runs, nwalkers, ndim)
#     """
#     params, j = [], 0
#     for i, par in enumerate(fundam_params):
#         p = np.array(par)
#         # If this parameter is one of the 'free' parameters.
#         if i in varIdxs:
#             # If it is the parameter metallicity, age or mass.
#             if i in pushidxs:
#                 pc = chains_nruns.T[j]
#                 chains = []
#                 for c in pc:
#                     chains.append(
#                         p[abs(c[None, :] - p[:, None]).argmin(axis=0)])
#                 params.append(chains)
#             else:
#                 params.append(chains_nruns.T[j])
#             j += 1

#     return np.array(params).T


def rangeCheck(model, ranges, varIdxs):
    """
    Check that all the model values are within the given ranges.
    """
    check_ranges = [
        r[0] <= p <= r[1] for p, r in zip(*[model, ranges[varIdxs]])]
    if all(check_ranges):
        return True
    return False


def r2Dist(fundam_params, varIdxs, params_trace):
    """
    R^2 for normal distribution.
    """
    param_r2 = []
    for i, _ in enumerate(fundam_params):
        if i in varIdxs:
            c_model = varIdxs.index(i)
            par = params_trace[c_model]
            param_r2.append(stats.probplot(par)[1][-1] ** 2)
        else:
            param_r2.append(np.nan)

    return param_r2


def modeKDE(fundam_params, varIdxs, mcmc_trace):
    """
    Estimate the mode for each fitted parameter from a KDE. Store the KDE
    for plotting.
    """
    mcmc_kde, mode_sol = [], []
    for i, mcmc_par in enumerate(mcmc_trace):
        mcmc_par = reject_outliers(mcmc_par)

        # Multiple elements are required
        if len(mcmc_par) <= 1:
            mcmc_kde.append([])
            mode_sol.append(np.nan)
            continue

        # Length of the last 10% of the chain.
        N = int(mcmc_par.shape[0] * .1)

        # Define KDE limits using only the last 10% of the chains.
        std = np.std(mcmc_par[-N:])
        pmin, pmax = np.min(mcmc_par[-N:]), np.max(mcmc_par[-N:])
        xp_min, xp_max = max(fundam_params[varIdxs[i]][0], pmin - std),\
            min(fundam_params[varIdxs[i]][-1], pmax + std)

        # KDE for plotting.
        try:
            # Use a slightly larger Scott bandwidth (looks better when plotted)
            bw = 1.25 * len(mcmc_par) ** (-1. / (len(varIdxs) + 4))
            x_kde, par_kde = kde1D(mcmc_par, xp_min, xp_max, bw)

            # Store for plotting
            mcmc_kde.append([x_kde, par_kde])
            # Mode (using KDE)
            mode_sol.append(x_kde[np.argmax(par_kde)])
        except (np.linalg.LinAlgError, FloatingPointError, UnboundLocalError):
            mcmc_kde.append([])
            mode_sol.append(np.nan)

    return mode_sol, mcmc_kde


def thinChain(mcmc_trace, acorr_t):
    """
    """
    return mcmc_trace[:, ::int(np.mean(acorr_t))]

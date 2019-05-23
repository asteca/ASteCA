
import numpy as np
import random
from scipy.optimize import differential_evolution as DE
from scipy import stats
import warnings
from ..synth_clust import synth_cluster
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


def closeSol(fundam_params, model, pushidxs):
    """
    Find the closest value in the parameters list for the discrete parameters
    metallicity, age, and mass.
    """
    model_proper = []
    for i, par in enumerate(fundam_params):
        # If it is the parameter metallicity, age or mass.
        if i in pushidxs:
            # Select the closest value in the array of allowed values.
            model_proper.append(min(
                par, key=lambda x: abs(x - model[i])))
        else:
            model_proper.append(model[i])

    return model_proper


def initPop(
    ranges, varIdxs, lkl_method, obs_clust, fundam_params,
        synthcl_args, ntemps, nwalkers_ptm, init_mode, popsize, maxiter):
    """
    Obtain initial parameter values using either a random distribution, or
    the Differential Evolution algorithm to approximate reasonable solutions.
    """

    p0 = []
    if init_mode == 'random':
        print("Random initial population")
        for _ in range(ntemps):
            p0.append(random_population(fundam_params, varIdxs, nwalkers_ptm))

    elif init_mode == 'diffevol':
        print("     DE init pop")

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # Estimate initial threshold value using DE.
            def DEdist(model):
                synth_clust = synthClust(
                    fundam_params, varIdxs, model, synthcl_args)
                if synth_clust:
                    lkl = likelihood.main(lkl_method, synth_clust, obs_clust)
                    return lkl
                return np.inf

            walkers_sols = []
            for _ in range(nwalkers_ptm):
                result = DE(
                    DEdist, ranges[varIdxs], popsize=popsize, maxiter=maxiter)
                walkers_sols.append(result.x)
                update_progress.updt(nwalkers_ptm, _ + 1)

        p0 = [walkers_sols for _ in range(ntemps)]

    return p0


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


def synthClust(fundam_params, varIdxs, model, synthcl_args):
    """
    Generate a synthetic cluster.
    """
    theor_tracks, e_max, err_lst, completeness, max_mag_syn, st_dist_mass,\
        R_V, ext_coefs, N_fc, cmpl_rnd, err_rnd = synthcl_args

    # Interpolate (z, a) data
    isochrone, model_proper = interpSol(
        theor_tracks, fundam_params, varIdxs, model)

    # Generate synthetic cluster.
    return synth_cluster.main(
        e_max, err_lst, completeness, max_mag_syn, st_dist_mass, isochrone,
        R_V, ext_coefs, N_fc, cmpl_rnd, err_rnd, model_proper)


def interpSol(theor_tracks, fundam_params, varIdxs, model):
    """
    Interpolate a new isochron from the four closest points in the (z, a) grid.

    The mass value is not interpolated.

    theor_tracks = [m1, m2, .., mN]
    mX = [age1, age2, ..., ageM]
    ageX = [f1,.., c1, c2,.., f1b,.., c1b, c2b,.., bp, mb, m_ini,.., m_bol]
    where:
    fX:  individual filters (mags)
    cX:  colors
    fXb: filters with binary data
    cXb: colors with the binary data
    bp:  binary probabilities
    mb:  binary masses
    m_ini,..., m_bol: six extra parameters.
    This means that '-6' is the index of the initial masses.

    """

    # Select the proper values for the model (ie: that exist in the grid), and
    # choose the (z, a) grid indexes to interpolate the isochrone.
    model_proper, j = [], 0
    for i, par in enumerate(fundam_params):
        # If this parameter is one of the 'free' parameters.
        if i in varIdxs:
            # If it is the parameter metallicity.
            if i == 0:
                # Select the closest value in the array of allowed values.
                mh = min(len(par) - 1, np.searchsorted(par, model[i - j]))
                ml = mh - 1
                model_proper.append(par[mh])
            elif i == 1:
                # Select the closest value in the array of allowed values.
                ah = min(len(par) - 1, np.searchsorted(par, model[i - j]))
                al = ah - 1
                model_proper.append(par[ah])
            elif i == 4:
                # Select the closest value in the array of allowed values for
                # the masses.
                model_proper.append(min(
                    par, key=lambda x: abs(x - model[i - j])))
            else:
                model_proper.append(model[i - j])
        else:
            if i == 0:
                ml = mh = 0
            elif i == 1:
                al = ah = 0
            model_proper.append(par[0])
            j += 1

    # # Minimum and maximum initial mass for each of the four isochrones.
    # mmin = np.min(isochs[:, -6, :], axis=1)
    # mmax = np.max(isochs[:, -6, :], axis=1)

    # Values of the four points in the (z, age) grid that contain the model
    # value (model[0], model[1])
    z1, z2 = fundam_params[0][ml], fundam_params[0][mh]
    a1, a2 = fundam_params[1][al], fundam_params[1][ah]
    pts = np.array([(z1, a1), (z1, a2), (z2, a1), (z2, a2)])

    # Define weights for the average, based on the inverse of the distance
    # to the grid (z, age) points.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # Fast euclidean distance: https://stackoverflow.com/a/47775357/1391441
        a_min_b = np.array([(model[0], model[1])]) - pts
        inv_d = 1. / np.sqrt(np.einsum('ij,ij->i', a_min_b, a_min_b))

        # Order: (z1, a1), (z1, a2), (z2, a1), (z2, a2)
        isochs = np.array([
            theor_tracks[ml][al], theor_tracks[ml][ah], theor_tracks[mh][al],
            theor_tracks[mh][ah]])

        # Sort according to smaller distance to the (model[0], model[1]) point.
        isoch_idx = np.argsort(inv_d)[::-1]

        # Masked weighted average. Source:
        # https://stackoverflow.com/a/35758345/1391441
        x1, x2, x3, x4 = isochs[isoch_idx]
        for x in (x2, x3, x4):
            # Maximum mass difference allowed
            msk = abs(x1[-6] - x[-6]) > .01
            # If the distance in this array is larger than the maximum allowed,
            # mask with the values taken from 'x1'.
            # x[:, msk] = x1[:, msk]
            np.copyto(x, x1, where=msk)

        # # Weighted average with mass "alignment".
        # isochrone = np.average(
        #     np.array([x1, x2, x3, x4]), weights=inv_d[isoch_idx], axis=0)
        # Scale weights so they add up to 1, then add based on them
        weights = inv_d[isoch_idx] / np.sum(inv_d[isoch_idx])
        isochrone = x1 * weights[0] + x2 * weights[1] + x3 * weights[2] +\
            x4 * weights[3]

    # This way is *marginally* faster
    # wgts = D * inv_d[isoch_idx]
    # x1, x2, x3, x4 = isochs[isoch_idx]
    # for x in (x2, x3, x4):
    #     # Maximum mass difference allowed
    #     msk = abs(x1[-6] - x[-6]) > .01
    #     x[:, msk] = x1[:, msk]
    # isochrone = np.sum(np.array([
    #     x1 * wgts[0], x2 * wgts[1], x3 * wgts[2], x4 * wgts[3]]), 0)

    # Weighted version, no mass alignment.
    # isochrone = np.average(np.array([
    #     theor_tracks[ml][al], theor_tracks[ml][ah], theor_tracks[mh][al],
    #     theor_tracks[mh][ah]]), weights=D * inv_d, axis=0)

    # Simpler mean version, no weights.
    # isochrone = np.mean([
    #     theor_tracks[ml][al], theor_tracks[ml][ah], theor_tracks[mh][al],
    #     theor_tracks[mh][ah]], axis=0)

    # nn = np.random.randint(0, 100)
    # if nn == 50:
    #     print(model)
    #     print(model_proper)
    #     import matplotlib.pyplot as plt
    #     plt.subplot(131)
    #     plt.scatter(*pts[isoch_idx][0], c='r')
    #     plt.scatter(*pts[isoch_idx][1:].T, c='g')
    #     plt.scatter(model[0], model[1], marker='x')
    #     # First color
    #     plt.subplot(132)
    #     plt.plot(theor_tracks[ml][al][1], theor_tracks[ml][al][0], c='b')
    #     plt.plot(theor_tracks[ml][ah][1], theor_tracks[ml][ah][0], c='r')
    #     plt.plot(theor_tracks[mh][al][1], theor_tracks[mh][al][0], c='cyan')
    #     plt.plot(theor_tracks[mh][ah][1], theor_tracks[mh][ah][0], c='orange')
    #     plt.plot(isochrone[1], isochrone[0], c='g', ls='--')
    #     plt.gca().invert_yaxis()
    #     # Second color
    #     plt.subplot(133)
    #     plt.plot(theor_tracks[ml][al][2], theor_tracks[ml][al][0], c='b')
    #     plt.plot(theor_tracks[ml][ah][2], theor_tracks[ml][ah][0], c='r')
    #     plt.plot(theor_tracks[mh][al][2], theor_tracks[mh][al][0], c='cyan')
    #     plt.plot(theor_tracks[mh][ah][2], theor_tracks[mh][ah][0], c='orange')
    #     plt.plot(isochrone[2], isochrone[0], c='g', ls='--')
    #     plt.gca().invert_yaxis()
    #     plt.show()

    return isochrone, model_proper


def discreteParams(fundam_params, varIdxs, chains_nruns, pushidxs):
    """
    Push values in each chain for each discrete parameter in the 'pushidxs'
    list to the closest grid value.

    chains_nruns.shape: (runs, nwalkers, ndim)
    """
    params, j = [], 0
    for i, par in enumerate(fundam_params):
        p = np.array(par)
        # If this parameter is one of the 'free' parameters.
        if i in varIdxs:
            # If it is the parameter metallicity, age or mass.
            if i in pushidxs:
                pc = chains_nruns.T[j]
                chains = []
                for c in pc:
                    chains.append(
                        p[abs(c[None, :] - p[:, None]).argmin(axis=0)])
                params.append(chains)
            else:
                params.append(chains_nruns.T[j])
            j += 1

    return np.array(params).T


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

        # Length of the last 10% of the chain.
        N = int(mcmc_par.shape[0] * .1)

        # Define KDE limits using only the last 10% of the chains.
        std = np.std(mcmc_par[-N:])
        pmin, pmax = np.min(mcmc_par[-N:]), np.max(mcmc_par[-N:])
        xp_min, xp_max = max(fundam_params[varIdxs[i]][0], pmin - std),\
            min(fundam_params[varIdxs[i]][-1], pmax + std)

        x_rang = .1 * (xp_max - xp_min)
        x_kde = np.mgrid[xp_min - x_rang:xp_max + x_rang:100j]
        # Use a slightly larger Scott bandwidth (looks better when plotted)
        bw = 1.25 * len(mcmc_par) ** (-1. / (len(varIdxs) + 4))
        # KDE for plotting.
        try:
            kernel_cl = stats.gaussian_kde(mcmc_par, bw_method=bw)
            # Parameter's KDE evaluated and reshaped.
            par_kde = np.reshape(kernel_cl(x_kde).T, x_kde.shape)

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

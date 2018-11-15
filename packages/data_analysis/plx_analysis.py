
import numpy as np
# from scipy import optimize
from scipy.optimize import differential_evolution as DE
from scipy.special import exp1
from ..best_fit.emcee3rc2 import ensemble
import warnings
from .. import update_progress


def main(clp):
    """
    """
    plx_flag = False
    mmag_clp, mp_clp, plx_clp, e_plx_clp, plx_Bys, plx_std_Bys, plx_wa =\
        [], [], [], [], np.nan, np.nan, np.nan

    # Extract parallax data.
    plx = np.array(list(zip(*list(zip(*clp['cl_reg_fit']))[7]))[0])
    # Array with no nan values
    plx_clrg = plx[~np.isnan(plx)]

    # Check that a range of parallaxes is possible.
    if plx_clrg.any() and np.min(plx_clrg) < np.max(plx_clrg):
        plx_flag = True
        print("  Bayesian Plx model")

        # Reject 2\sigma outliers.
        max_plx, min_plx = np.nanmedian(plx) + 2. * np.nanstd(plx),\
            np.nanmedian(plx) - 2. * np.nanstd(plx)

        # Suppress Runtimewarning issued when 'plx' contains 'nan' values.
        with np.warnings.catch_warnings():
            np.warnings.filterwarnings('ignore')
            plx_2s_msk = (plx < max_plx) & (plx > min_plx)

        # Prepare masked data.
        mmag_clp = np.array(
            list(zip(*list(zip(*clp['cl_reg_fit']))[3]))[0])[plx_2s_msk]
        mp_clp = np.array(list(zip(*clp['cl_reg_fit']))[9])[plx_2s_msk]
        plx_clp = plx[plx_2s_msk]
        e_plx_clp = np.array(
            list(zip(*list(zip(*clp['cl_reg_fit']))[8]))[0])[plx_2s_msk]
        # Take care of possible zero values that can produce issues since
        # errors are in the denominator.
        e_plx_clp[e_plx_clp == 0.] = 10.

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # Use DE to estimate the ML
            def DEdist(model):
                return -lnlike(model, plx_clp, e_plx_clp, mp_clp)
            # bounds = [[0., plx_clp.max()], [0., 1.]]
            bounds = [[0., 20.]]
            result = DE(DEdist, bounds, popsize=20, maxiter=50)
            # print(result)

        # Prior parameters.
        mu_p = result.x
        # Sampler parameters.
        ndim, nwalkers, nruns = 1, 50, 1000
        sampler = ensemble.EnsembleSampler(
            nwalkers, ndim, lnprob, args=(plx_clp, e_plx_clp, mp_clp, mu_p))

        # Random initial guesses.
        # pos0 = [np.random.uniform(0., 1., ndim) for i in range(nwalkers)]
        # Ball of initial guesses around 'mu_p'
        pos0 = [mu_p + .5 * np.random.normal() for i in range(nwalkers)]

        old_tau = np.inf
        for i, _ in enumerate(sampler.sample(pos0, iterations=nruns)):
            # Only check convergence every X steps
            if i % 50 and i < (nruns - 1):
                continue

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                tau = sampler.get_autocorr_time(tol=0)
                # Check convergence
                converged = np.all(tau * 100 < i)
                converged &= np.all(np.abs(old_tau - tau) / tau < 0.01)
                if converged:
                    print("")
                    break
                old_tau = tau
                tau_mean = np.nanmean(sampler.get_autocorr_time(tol=0))
            update_progress.updt(nruns, i + 1, "tau={:.0f}".format(tau_mean))

        # Remove burn-in (half of chain)
        nburn = int(i / 2.)
        samples = sampler.get_chain(discard=nburn, flat=True)

        # import matplotlib.pyplot as plt
        # import corner
        # corner.corner(samples)
        # plt.show()

        # plt.plot(samples.T[0])
        # plt.show()

        # 16th, median, 84th
        plx_Bys = np.percentile(samples, (16, 50, 84))
        tau_mean = np.mean(sampler.get_autocorr_time(tol=0))
        print("Bayesian plx estimated: {:.3f} (ESS={:.0f}, tau={:.0f})".format(
            1. / plx_Bys[1], samples.size / tau_mean, tau_mean))

        # Weighted average and its error.
        # Source: https://physics.stackexchange.com/a/329412/8514
        plx_w = mp_clp / np.square(e_plx_clp)
        # e_plx_w = np.sqrt(np.sum(np.square(e_plx * plx_w))) / np.sum(plx_w)
        plx_wa = np.average(plx_clp, weights=plx_w)

    clp.update({
        'plx_flag': plx_flag, 'plx_clrg': plx_clrg, 'mmag_clp': mmag_clp,
        'mp_clp': mp_clp, 'plx_clp': plx_clp, 'e_plx_clp': e_plx_clp,
        'plx_Bys': plx_Bys, 'plx_wa': plx_wa})
    return clp


def lnprob(mu_std, plx, plx_std, mp, mu_std_p):
    lp = lnprior(mu_std, mu_std_p)
    if np.isinf(lp):
        return lp
    return lp + lnlike(mu_std, plx, plx_std, mp)


def lnprior(mu_std, mu_std_p, std_p=(1., .1)):
    """
    Log prior, Gaussian > 0.
    """
    # if mu_std[0] < 0. or mu_std[1] < 0.:
    #     return -np.inf
    # return -0.5 * np.sum(
    #     ((mu_std[0] - mu_std_p[0]) / std_p)**2 +
    #     ((mu_std[1] - mu_std_p[1]) / std_p)**2)

    if mu_std < 0.: #[0] < 0. or mu_std[1] < 0.:
        return -np.inf
    return -0.5 * np.sum(
        ((mu_std - mu_std_p) / std_p[0])**2)# +
        # ((mu_std[1] - mu_std_p[1]) / std_p[1])**2)


# def lnlike(mu_std, plx, plx_std, mp):
#     """
#     Log likelihood, product of Gaussian functions.
#     """
#     # sigma = np.sqrt(plx_std**2 + mu_std[1]**2)
#     # A = -np.sum(np.log(np.sqrt(2. * np.pi) * sigma))
#     # B = -0.5 * (np.sum(mp * (plx - mu_std[0])**2 / sigma**2))

#     from scipy.integrate import quad

#     # Bailer-jones
#     A = -np.sum(np.log(2. * np.pi * plx_std * mu_std[1]))

#     int_exp = np.zeros(plx.size)
#     for i in range(plx.size):
#         int_exp[i] = quad(
#             distFunc, 0.1, 20., args=(plx[i], plx_std[i], mp[i], mu_std))[0]
#     B = np.sum(np.log(int_exp))

#     return A + B


# def distFunc(r_i, plx, plx_std, mp, mu_std):
#     B1 = ((plx - (1. / r_i)) / plx_std)**2
#     B2 = ((r_i - mu_std[0]) / mu_std[1])**2
#     C = np.exp(-.5 * 1. * (B1 + B2))
#     return C

# def lnlike(mu_std, plx, plx_std, mp):
#     """
#     Log likelihood, product of Gaussian functions.
#     """
#     # Bailer-jones
#     A = -np.sum(np.log(2. * np.pi * plx_std * mu_std[1]))

#     int_max = mu_std[0] + 3 * mu_std[1]
#     N = int(int_max / 0.01)

#     x = np.linspace(.1, int_max, N).reshape(-1, 1)
#     B1 = ((plx - (1. / x)) / plx_std)**2

#     def distFunc0(r_i):
#         B2 = ((r_i - mu_std[0]) / mu_std[1])**2
#         C = np.exp(-.5 * 1. * (B1 + B2))
#         return C

#     int_exp = np.trapz(distFunc0(x), x, axis=0)
#     B = np.sum(np.log(int_exp))

#     return A + B


def lnlike(mu, plx, plx_std, mp):

    # Define the 'r_i' values used to evaluate the integral.
    int_max = mu + 5.
    N = int(int_max / 0.01)
    x = np.linspace(.1, int_max, N).reshape(-1, 1)

    # Marginalization of the scale parameter 's_c'. We integrate over it
    # using the incomplete gamma function as per Wolfram:
    #
    # https://www.wolframalpha.com/input/
    # ?i=integral+exp(-.5*(a%5E2%2Fx%5E2))+%2F+x,+x%3D0+to+x%3Db
    #
    # This function is equivalent to scipy's 'exp1()', as stated in:
    #
    # https://stackoverflow.com/a/53148269/1391441
    #
    # so we use this function to marginalize the 's_c' parameter up to a
    # 5 kpc limit.
    lim_u = 5.

    # This is the first term in Eq (20) of Bailer-Jones (2015). We calculate
    # it outside of the integral as it only depends on the 'r_i' values (ie:
    # the 'x' array)

    def distFunc(r_i):
        sc_int = .5 * exp1(.5 * ((r_i - mu) / lim_u)**2)
        B1 = ((plx - (1. / r_i)) / plx_std)**2
        C = (np.exp(-.5 * B1) / plx_std) * sc_int
        return C

    # Double integral
    int_exp = np.trapz(distFunc(x), x, axis=0)

    return np.sum(np.log(mp * int_exp))

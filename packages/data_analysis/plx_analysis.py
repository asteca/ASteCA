
import numpy as np
from scipy.optimize import differential_evolution as DE
from scipy.special import exp1
from ..best_fit.emcee3rc2 import ensemble
import warnings
from .. import update_progress


def main(
    clp, plx_bayes_flag, plx_offset, plx_chains, plx_runs, flag_plx_mp,
        flag_make_plot, **kwargs):
    """
    Bayesian parallax distance using the Bailer-Jones (2015) model with the
    'shape parameter' marginalized.

    Hardcoded choices:

    * 2 sigma outliers are rejected
    * MPs are used
    * Bayesian prior is a Gaussian with a fixed standard deviation
    """

    plx_clrg, mmag_clp, mp_clp, plx_clp, e_plx_clp, plx_samples,\
        plx_tau_autocorr, mean_afs = [[] for _ in range(8)]
    plx_flag_clp, plx_bayes_flag_clp, plx_wa, plx_Bys, plx_ess = False, False,\
        np.nan, np.array([]), np.nan

    if ('C2' in flag_make_plot) or plx_bayes_flag:

        # Extract parallax data.
        plx = np.array(list(zip(*list(zip(*clp['cl_reg_fit']))[7]))[0])
        # Array with no nan values
        plx_clrg = plx[~np.isnan(plx)]

        plx_flag_clp = checkPlx(plx_clrg)

        if plx_flag_clp:
            print("Processing parallaxes")

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

            # Weighted average.
            # Source: https://physics.stackexchange.com/a/329412/8514
            plx_w = mp_clp / np.square(e_plx_clp)
            plx_w = plx_w if plx_w > 0. else None
            plx_wa = np.average(plx_clp, weights=plx_w)

            if plx_bayes_flag:
                plx_samples, plx_Bys, plx_bayes_flag_clp, plx_tau_autocorr,\
                    mean_afs, plx_ess = plxBayes(
                        plx_offset, plx_chains, plx_runs, flag_plx_mp, plx_clp,
                        e_plx_clp, mp_clp)

        else:
            print("  WARNING: no valid Plx data found")

    clp.update({
        'plx_flag_clp': plx_flag_clp, 'plx_clrg': plx_clrg,
        'mmag_clp': mmag_clp, 'mp_clp': mp_clp, 'plx_clp': plx_clp,
        'e_plx_clp': e_plx_clp, 'plx_Bys': plx_Bys, 'plx_wa': plx_wa,
        'plx_bayes_flag_clp': plx_bayes_flag_clp, 'plx_samples': plx_samples,
        'plx_tau_autocorr': plx_tau_autocorr, 'mean_afs': mean_afs,
        'plx_ess': plx_ess})
    return clp


def checkPlx(plx_clrg):
    """
    Check that a range of parallaxes is possible.
    """
    if plx_clrg.any() and np.min(plx_clrg) < np.max(plx_clrg):
        return True
    else:
        return False


def plxBayes(
    plx_offset, plx_chains, plx_runs, flag_plx_mp, plx_clp, e_plx_clp,
        mp_clp):
    """
    """
    plx_bayes_flag_clp = True

    # Add offset to parallax data.
    plx_clp += plx_offset

    # If MPs are not to be used, use all 1.
    if flag_plx_mp is False:
        mp_clp = np.ones(mp_clp.shape)

    # Sampler parameters.
    ndim, nwalkers, nruns = 1, plx_chains, plx_runs
    print("  Bayesian Plx model ({} runs)".format(nruns))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        # Define the 'r_i' values used to evaluate the integral.
        int_max = 20.
        N = int(int_max / 0.01)
        x = np.linspace(.1, int_max, N).reshape(-1, 1)
        B1 = ((plx_clp - (1. / x)) / e_plx_clp)**2
        B2 = (np.exp(-.5 * B1) / e_plx_clp)

        # Use DE to estimate the ML
        def DEdist(model):
            return -lnlike(model, x, B2, mp_clp)
        bounds = [[0., 20.]]
        result = DE(DEdist, bounds, popsize=20, maxiter=100)
        # print(result)

    # Prior parameters.
    mu_p = result.x
    # Define the 'r_i' values used to evaluate the integral.
    int_max = mu_p + 5.
    N = int(int_max / 0.01)
    x = np.linspace(.1, int_max, N).reshape(-1, 1)
    B1 = ((plx_clp - (1. / x)) / e_plx_clp)**2
    B2 = (np.exp(-.5 * B1) / e_plx_clp)

    # emcee sampler
    sampler = ensemble.EnsembleSampler(
        nwalkers, ndim, lnprob, args=(x, B2, mp_clp, mu_p))

    # Random initial guesses.
    # pos0 = [np.random.uniform(0., 1., ndim) for i in range(nwalkers)]
    # Ball of initial guesses around 'mu_p'
    pos0 = [mu_p + .5 * np.random.normal() for i in range(nwalkers)]

    tau_index, autocorr_vals, afs = 0, np.empty(nruns), np.empty(nruns)

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            old_tau, N_conv = np.inf, 1000
            for i, _ in enumerate(sampler.sample(pos0, iterations=nruns)):
                # Only check convergence every X steps
                if i % 10 and i < (nruns - 1):
                    continue

                afs[tau_index] = np.mean(sampler.acceptance_fraction)

                tau = sampler.get_autocorr_time(tol=0)
                autocorr_vals[tau_index] = np.mean(tau)
                tau_index += 1

                # Check convergence
                converged = np.all(tau * N_conv < i)
                converged &= np.all(np.abs(old_tau - tau) / tau < 0.01)
                if converged:
                    print("")
                    break
                old_tau = tau
                # tau_mean = np.nanmean(sampler.get_autocorr_time(tol=0))
                update_progress.updt(nruns, i + 1)

        mean_afs = afs[:tau_index]
        tau_autocorr = autocorr_vals[:tau_index]
        # Remove burn-in (25% of chain)
        nburn = int(i * .25)
        samples = sampler.get_chain(discard=nburn, flat=True)
        # 16th, 84th in Kpc
        p16, p84 = np.percentile(samples, (16, 84))
        plx_Bys = np.array([p16, np.mean(samples), p84])
        tau_mean = np.mean(sampler.get_autocorr_time(tol=0))
        plx_ess = samples.size / tau_mean

        # For plotting
        plx_samples = sampler.get_chain()[:, :, 0]

        print("Bayesian plx estimated: " +
              "{:.3f} (ESS={:.0f}, tau={:.0f})".format(
                  1. / plx_Bys[1], plx_ess, tau_mean))
    except Exception as e:
        print(e)
        print("\n  ERROR: could not process Plx data with emcee")
        plx_samples, plx_Bys, plx_bayes_flag_clp, plx_ess = [], np.array([]),\
            False, np.nan

    return plx_samples, plx_Bys, plx_bayes_flag_clp, tau_autocorr,\
        mean_afs, plx_ess


def lnprob(mu, x, B2, MPs, mu_p):
    lp = lnprior(mu, mu_p)
    if np.isinf(lp):
        return lp
    return lp + lnlike(mu, x, B2, MPs)


def lnprior(mu, mu_p, std_p=1.):
    """
    Log prior.
    """
    if mu < 0.:
        return -np.inf

    # Gaussian > 0
    return -0.5 * ((mu - mu_p) / std_p)**2

    # Exponential prior proposed by Bailer-Jones.
    # return (.5 / (.5 * mu_p)**3) * mu**2 * np.exp(-mu / (.5 * mu_p))

    # # Uniform prior
    # return 0.


def lnlike(mu, x, B2, MPs):
    """
    Model defined in Bailer-Jones (2015), Eq (20), The shape parameter s_c
    is marginalized.
    """

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

    def distFunc(r_i):
        """
        Eq (20) of Bailer-Jones (2015) with everything that can be calculated
        outside, moved outside.
        """
        sc_int = .5 * exp1(.5 * ((r_i - mu) / lim_u)**2)
        sc_int.T[np.isinf(sc_int.T)] = 0.
        return B2 * sc_int

    # Double integral
    int_exp = np.trapz(distFunc(x), x, axis=0)

    # Apply MPs
    MP_intexp = MPs * int_exp
    # Mask 'bad' values
    msk = np.logical_or(
        MP_intexp <= 0., MP_intexp >= np.inf, MP_intexp <= -np.inf)
    MP_intexp[msk] = np.nan

    return np.nansum(np.log(MP_intexp))

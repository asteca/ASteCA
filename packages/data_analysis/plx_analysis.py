
import numpy as np
from scipy.optimize import differential_evolution as DE
from scipy.special import exp1
try:
    # If this import is not done outside main(), then eval() fails in the
    # definition of the moves
    from emcee import moves
except ImportError:
    pass
import warnings
from ..best_fit.bf_common import modeKDE
from .. import update_progress


def main(
    clp, plx_bayes_flag, plx_offset, plx_chains, plx_runs, plx_burn,
        plx_emcee_moves, flag_make_plot, outlr_std=3., **kwargs):
    """
    Bayesian parallax distance using the Bailer-Jones (2015) model with the
    'shape parameter' marginalized.

    Hardcoded choices:

    * outlr_std sigma outliers are rejected
    * Bayesian prior is a Gaussian with a fixed standard deviation
    """

    plx_clrg, mmag_clp, plx_clp, e_plx_clp, plx_samples,\
        plx_tau_autocorr, mean_afs = [[] for _ in range(7)]
    plx_flag_clp, plx_bayes_flag_clp, plx_wa, plx_Bayes_kde, plx_Bys,\
        plx_ess = False, False, np.nan, np.array([]), np.array([]), np.nan

    if ('C2' in flag_make_plot) or plx_bayes_flag:

        # Extract parallax data.
        plx = np.array(list(zip(*list(zip(*clp['cl_reg_fit']))[7]))[0])
        # Array with no nan values
        plx_clrg = plx[~np.isnan(plx)]

        plx_flag_clp = checkPlx(plx_clrg)

        if plx_flag_clp:
            print("Processing parallaxes")

            # Reject outlr_std*\sigma outliers.
            max_plx = np.nanmedian(plx) + outlr_std * np.nanstd(plx)
            min_plx = np.nanmedian(plx) - outlr_std * np.nanstd(plx)

            # Suppress Runtimewarning issued when 'plx' contains 'nan'
            # values.
            with np.warnings.catch_warnings():
                np.warnings.filterwarnings('ignore')
                plx_2s_msk = (plx < max_plx) & (plx > min_plx)

            # Prepare masked data.
            mmag_clp = np.array(
                list(zip(*list(zip(
                    *clp['cl_reg_fit']))[3]))[0])[plx_2s_msk]
            plx_clp = plx[plx_2s_msk]
            e_plx_clp = np.array(
                list(zip(*list(zip(
                    *clp['cl_reg_fit']))[8]))[0])[plx_2s_msk]
            # Take care of possible zero values that can produce issues
            # since errors are in the denominator.
            e_plx_clp[e_plx_clp == 0.] = 10.

            # Weighted average.
            # Source: https://physics.stackexchange.com/a/329412/8514
            plx_w = 1. / np.square(e_plx_clp)
            plx_w = plx_w if plx_w.sum() > 0. else None
            plx_wa = np.average(plx_clp, weights=plx_w)

            if plx_bayes_flag:
                plx_samples, plx_Bayes_kde, plx_Bys, plx_bayes_flag_clp,\
                    plx_tau_autocorr, mean_afs, plx_ess = plxBayes(
                        plx_offset, plx_chains, plx_runs, plx_burn,
                        plx_emcee_moves, plx_clp, e_plx_clp)

        else:
            print("  WARNING: no valid Plx data found")

    clp.update({
        'plx_flag_clp': plx_flag_clp, 'plx_clrg': plx_clrg,
        'mmag_clp': mmag_clp, 'plx_clp': plx_clp,
        'e_plx_clp': e_plx_clp, 'plx_Bys': plx_Bys, 'plx_wa': plx_wa,
        'plx_bayes_flag_clp': plx_bayes_flag_clp, 'plx_samples': plx_samples,
        'plx_Bayes_kde': plx_Bayes_kde, 'plx_tau_autocorr': plx_tau_autocorr,
        'mean_afs': mean_afs, 'plx_ess': plx_ess})
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
    plx_offset, plx_chains, plx_runs, plx_burn, plx_emcee_moves,
        plx_clp, e_plx_clp, N_conv=1000, tau_stable=0.05):
    """

    HARDCODED
    N_conv
    tau_stable
    """

    from emcee import ensemble
    # Move used by emcee
    mv = [(eval("(moves." + _ + ")")) for _ in plx_emcee_moves]

    plx_bayes_flag_clp = True

    # Add offset to parallax data.
    plx_clp += plx_offset

    # Sampler parameters.
    ndim, nwalkers, nruns = 1, plx_chains, plx_runs
    print("  Bayesian Plx model ({} runs)".format(nruns))

    # DE initial mu position
    # mu_p = DE_mu_sol(plx_clp, e_plx_clp)
    mu_p = np.mean(plx_clp)

    # Define the 'r_i' values used to evaluate the integral.
    int_max = mu_p + 5.
    x, B2 = r_iVals(int_max, plx_clp, e_plx_clp)

    # emcee sampler
    sampler = ensemble.EnsembleSampler(
        nwalkers, ndim, lnprob, args=(x, B2, mu_p), moves=mv)

    # Ball of initial guesses around 'mu_p'
    # pos0 = np.clip(
    #     np.array([[mu_p + .05 * np.random.normal()] for i in range(nwalkers)]),
    #     a_min=0., a_max=None)
    # Random initial guesses
    pos0 = np.random.uniform(0., 2. * np.mean(plx_clp), (nwalkers, 1))

    tau_index, autocorr_vals, afs = 0, np.empty(nruns), np.empty(nruns)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # HARDCODED
            old_tau = np.inf
            for i, _ in enumerate(sampler.sample(pos0, iterations=nruns)):
                # Only check convergence every X steps
                if i % 10 and i < (nruns - 1):
                    continue

                afs[tau_index] = np.mean(sampler.acceptance_fraction)

                tau = sampler.get_autocorr_time(tol=0)
                autocorr_vals[tau_index] = np.mean(tau)
                tau_index += 1

                # Check convergence
                converged = tau * (N_conv / nwalkers) < i * plx_burn
                converged &= np.all(np.abs(old_tau - tau) / tau < tau_stable)
                if converged:
                    print("")
                    break
                old_tau = tau
                update_progress.updt(nruns, i + 1)

        mean_afs = afs[:tau_index]
        tau_autocorr = autocorr_vals[:tau_index]

        nburn = int(i * plx_burn)
        samples = sampler.get_chain(discard=nburn, flat=True)

        # Mode and KDE to plot
        # This simulates the 'fundam_params and 'varIdxs' arrays.
        fp, vi = [[-np.inf, np.inf], [-np.inf, np.inf]], [0, 1]
        plx_Bys_mode, plx_Bayes_kde = modeKDE(fp, vi, 1. / samples.T)
        plx_Bys_mode, plx_Bayes_kde = 1. / plx_Bys_mode[0], plx_Bayes_kde[0]

        # 16th, median, 84th, mean, mode in Kpc
        p16, p50, p84 = np.percentile(samples, (16, 50, 84))
        plx_Bys = np.array([p16, p50, p84, np.mean(samples), plx_Bys_mode])

        tau = sampler.get_autocorr_time(tol=0)[0]
        plx_ess = samples.size / tau

        # For plotting, (nsteps, nchains, ndim)
        plx_samples = sampler.get_chain()[:, :, 0]

        print("Bayesian plx estimated: "
              + "{:.3f} (ESS={:.0f}, tau={:.0f})".format(
                  1. / plx_Bys[3], plx_ess, tau))
    except Exception as e:
        print(e)
        print("\n  ERROR: could not process Plx data with emcee")
        plx_samples, plx_Bayes_kde, plx_Bys, plx_bayes_flag_clp, plx_ess,\
            tau_autocorr, mean_afs = [], np.array([]), np.array([]), False,\
            np.nan, np.nan, np.nan

    return plx_samples, plx_Bayes_kde, plx_Bys, plx_bayes_flag_clp,\
        tau_autocorr, mean_afs, plx_ess


def DE_mu_sol(plx_clp, e_plx_clp, int_max=20., psize=20, maxi=100):
    """
    Use the Differential Evolution algorithm to approximate the best solution
    used as the mean of the prior.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        # Define the 'r_i' values used to evaluate the integral.
        x, B2 = r_iVals(int_max, plx_clp, e_plx_clp)

        # Use DE to estimate the ML
        def DEdist(model):
            return -lnlike(model, x, B2)
        bounds = [[0., 20.]]
        result = DE(DEdist, bounds, popsize=psize, maxiter=maxi)

    return result.x[0]


def r_iVals(int_max, plx, e_plx):
    """
    The 'r_i' values used to evaluate the integral.
    """
    N = int(int_max / 0.01)
    x = np.linspace(.1, int_max, N).reshape(-1, 1)
    B1 = ((plx - (1. / x)) / e_plx)**2
    B2 = (np.exp(-.5 * B1) / e_plx)
    return x, B2


def lnprob(mu, x, B2, mu_p):
    lp = lnprior(mu, mu_p)
    if np.isinf(lp):
        return -np.inf
    return lp + lnlike(mu, x, B2)


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


def lnlike(mu, x, B2):
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

    # Mask 'bad' values
    msk = np.logical_or(
        int_exp <= 0., int_exp >= np.inf, int_exp <= -np.inf)
    int_exp[msk] = np.nan

    return np.nansum(np.log(int_exp))

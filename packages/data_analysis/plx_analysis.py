
import numpy as np
import warnings
from ..best_fit.bf_common import modeKDE
from .. import update_progress


def main(
    clp, plx_bayes_flag, plx_offset, plx_chains, plx_runs, plx_burn,
        flag_make_plot, outlr_std=3., **kwargs):
    """
    Bayesian parallax distance using the simple likelihood shown in
    Cantat-Gaudin et al. (2018), 'A Gaia DR2 view of the Open Cluster
    population in the Milky Way'

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
                        plx_clp, e_plx_clp)

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
    plx_offset, plx_chains, plx_runs, plx_burn,
        plx_clp, e_plx_clp, N_conv=1000, tau_stable=0.05):
    """

    HARDCODED
    N_conv
    tau_stable
    """

    from emcee import ensemble

    plx_bayes_flag_clp = True

    # Add offset to parallax data.
    plx_clp += plx_offset

    # Sampler parameters.
    ndim, nwalkers, nruns = 1, plx_chains, plx_runs
    print("  Bayesian Plx model ({} runs)".format(nruns))

    e_plx2 = e_plx_clp**2
    # emcee sampler
    sampler = ensemble.EnsembleSampler(
        nwalkers, ndim, lnprob, args=(plx_clp, e_plx2))

    # Ball of initial guesses around 'mu_p'
    # pos0 = np.clip(
    #   np.array([[mu_p + .05 * np.random.normal()] for i in range(nwalkers)]),
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


def lnprob(d, plx, e_plx2):
    if d <= 0:
        return -np.inf
    return lnlike(d, plx, e_plx2)


def lnlike(d, plx, e_plx2):
    """
    Simple likelihood used in Cantat-Gaudin et al. (2018), 'A Gaia DR2 view of
    the Open Cluster population in the Milky Way'

    The final estimated value is almost always equivalent to the weighted
    average.
    """
    return -np.sum((plx - 1 / d)**2 / e_plx2)

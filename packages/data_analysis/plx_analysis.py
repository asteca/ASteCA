
import numpy as np
from scipy import optimize
from ..best_fit.emcee3rc2 import ensemble


def main(clp):
    """
    """
    plx_flag = False
    mmag_clp, mp_clp, plx_clp, e_plx_clp, pl_plx, plx_bay, ph_plx, plx_wa =\
        [], [], [], [], np.nan, np.nan, np.nan, np.nan

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

        # Use optimum likelihood value as mean of the prior.
        def pstv_lnlike(w_t, w_i, s_i, mp):
            return -lnlike(w_t, w_i, s_i, mp)
        plx_lkl = optimize.minimize_scalar(pstv_lnlike, args=(
            plx_clp, e_plx_clp, 1.))

        # Prior parameters.
        w_p, s_p = plx_lkl.x, .5
        # Sampler parameters.
        ndim, nwalkers, nruns, nburn = 1, 100, 2000, 1000
        sampler = ensemble.EnsembleSampler(
            nwalkers, ndim, lnprob,
            args=(plx_clp, e_plx_clp, mp_clp, w_p, s_p))
        # Random initial guesses.
        pos = [np.random.uniform(0., 1., ndim) for i in range(nwalkers)]
        sampler.run_mcmc(pos, nruns)
        # Remove burn-in
        samples = sampler.chain[:, nburn:, :].reshape((-1, ndim))
        # 16th, median, 84th percentiles
        pl_plx, plx_bay, ph_plx = np.percentile(samples, [16, 50, 84])
        print("Bayesian plx estimated: {:.3f} (ESS={:.0f})".format(
            plx_bay, samples.size / sampler.get_autocorr_time()[0]))

        # Weighted average and its error.
        # Source: https://physics.stackexchange.com/a/329412/8514
        plx_w = mp_clp / np.square(e_plx_clp)
        # e_plx_w = np.sqrt(np.sum(np.square(e_plx * plx_w))) / np.sum(plx_w)
        plx_wa = np.average(plx_clp, weights=plx_w)

    clp.update({
        'plx_flag': plx_flag, 'plx_clrg': plx_clrg, 'mmag_clp': mmag_clp,
        'mp_clp': mp_clp, 'plx_clp': plx_clp, 'e_plx_clp': e_plx_clp,
        'pl_plx': pl_plx, 'plx_bay': plx_bay, 'ph_plx': ph_plx,
        'plx_wa': plx_wa})
    return clp


def lnprob(w_t, w_i, s_i, mp, w_p, s_p):
    lp = lnprior(w_t, w_p, s_p)
    return lp + lnlike(w_t, w_i, s_i, mp)


def lnprior(w_t, w_p, s_p):
    """
    Log prior, Gaussian > 0.
    """
    if w_t < 0.:
        return -np.inf
    return -0.5 * ((w_p - w_t)**2 / s_p**2)


def lnlike(w_t, w_i, s_i, mp):
    """
    Log likelihood, product of Gaussian functions.
    """
    return -0.5 * (np.sum(mp * (w_i - w_t)**2 / s_i**2))

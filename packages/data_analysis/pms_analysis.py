
import numpy as np
from scipy.stats import pearsonr
from scipy.optimize import differential_evolution as DE
from ..best_fit.emcee3rc2 import ensemble
import warnings
from .. import update_progress


def main(clp, coords, **kwargs):
    """
    """
    PM_flag = False
    pmMP, pmRA_DE, e_pmRA_DE, pmDE, e_pmDE, DE_pm, mmag_pm, pmRA_Bys,\
        pmDE_Bys, pmRA_std_Bys, pmDE_std_Bys = [[]] * 11

    # Extract RA PMs data.
    pmRA = np.array(list(zip(*list(zip(*clp['cl_reg_fit']))[7]))[1])
    # Array with no nan values
    pmRA_clrg = pmRA[~np.isnan(pmRA)]

    # Check that PMs were defined within the cluster region.
    if pmRA_clrg.any():
        PM_flag = True
        print("  Bayesian PMs model")

        # Cluster region data.
        pmMP, pmRA, e_pmRA, pmDE, e_pmDE =\
            np.array(list(zip(*clp['cl_reg_fit']))[9]),\
            np.array(list(zip(*list(zip(*clp['cl_reg_fit']))[7]))[1]),\
            np.array(list(zip(*list(zip(*clp['cl_reg_fit']))[8]))[1]),\
            np.array(list(zip(*list(zip(*clp['cl_reg_fit']))[7]))[2]),\
            np.array(list(zip(*list(zip(*clp['cl_reg_fit']))[8]))[2])
        DE_pm = np.array(list(zip(*clp['cl_reg_fit']))[2]) if coords == 'deg'\
            else np.zeros(pmRA.size)
        mmag_pm = np.array(list(zip(*list(zip(*clp['cl_reg_fit']))[3]))[0])

        # Remove nan values from cluster region
        msk = ~np.isnan(pmRA) & ~np.isnan(e_pmRA) & ~np.isnan(pmDE) &\
            ~np.isnan(e_pmDE)
        pmMP, pmRA, e_pmRA, pmDE, e_pmDE, DE_pm, mmag_pm = pmMP[msk],\
            pmRA[msk], e_pmRA[msk], pmDE[msk], e_pmDE[msk], DE_pm[msk],\
            mmag_pm[msk]

        pmRA_DE = pmRA * np.cos(np.deg2rad(DE_pm))
        # Propagate error in RA*cos(delta)
        e_pmRA_DE = np.sqrt(
            (e_pmRA * pmRA * np.sin(np.deg2rad(DE_pm)))**2 +
            (e_pmDE * np.cos(np.deg2rad(DE_pm)))**2)
        # Pearson corr coeff
        rho = pearsonr(pmRA_DE, pmDE)[0]

        # Use DE to estimate the ML
        def DEdist(model):
            return -lnlike(
                model, pmRA_DE, e_pmRA_DE, pmDE, e_pmDE, pmMP, rho)
        bounds = [
            [pmRA_DE.min(), pmRA_DE.max()], [pmDE.min(), pmDE.max()],
            [0., 5.], [0., 5.]]
        result = DE(DEdist, bounds, popsize=10, maxiter=1000)
        # print(result)

        # Prior parameters.
        mu_std_p = result.x
        # Sampler parameters.
        ndim, nwalkers, nruns, nburn = 4, 100, 5000, 2500
        sampler = ensemble.EnsembleSampler(
            nwalkers, ndim, lnprob,
            args=(pmRA_DE, e_pmRA_DE, pmDE, e_pmDE, pmMP, rho, mu_std_p))
        # Random initial guesses.
        pos0 = [np.random.uniform(0., 1., ndim) for i in range(nwalkers)]
        old_tau = np.inf
        for i, _ in enumerate(sampler.sample(pos0, iterations=nruns)):

            # Only check convergence every 100 steps
            if i % 50:
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

            update_progress.updt(nruns, i + 1)

        # Remove burn-in
        samples = sampler.get_chain(discard=nburn, flat=True)
        # Medians
        pmRA_Bys, pmDE_Bys, pmRA_std_Bys, pmDE_std_Bys = np.percentile(
            samples, 50, axis=0).T
        tau_mean = np.mean(sampler.get_autocorr_time(tol=0))
        print(("Bayesian PMs estimated: ({:.3f}, {:.3f})"
               " (ESS={:.0f}, tau={:.0f})").format(
            pmRA_Bys, pmDE_Bys, samples.size / tau_mean, tau_mean))

        # import matplotlib.pyplot as plt
        # import corner
        # corner.corner(samples)
        # plt.show()
        # plt.subplot(411)
        # plt.plot(samples.T[0])
        # plt.subplot(412)
        # plt.plot(samples.T[1])
        # plt.subplot(413)
        # plt.plot(samples.T[2])
        # plt.subplot(414)
        # plt.plot(samples.T[3])
        # plt.show()

    clp.update({
        'PM_flag': PM_flag, 'pmMP': pmMP, 'pmRA_DE': pmRA_DE,
        'e_pmRA_DE': e_pmRA_DE, 'pmDE': pmDE, 'e_pmDE': e_pmDE, 'DE_pm': DE_pm,
        'mmag_pm': mmag_pm, 'pmRA_Bys': (pmRA_Bys, pmRA_std_Bys),
        'pmDE_Bys': (pmDE_Bys, pmDE_std_Bys)})
    return clp


def lnprob(mu_std_xy, x, ex, y, ey, mp, rho_term, mu_std_p):
    lp = lnprior(mu_std_xy, mu_std_p)
    return lp + lnlike(mu_std_xy, x, ex, y, ey, mp, rho_term)


def lnprior(mu_std_xy, mu_std_p, std_p=5.):
    """
    Log prior.
    """
    if mu_std_xy[2] < 0. or mu_std_xy[3] < 0.:
        return -np.inf
    return -0.5 * np.sum(
        ((mu_std_xy[0] - mu_std_p[0]) / std_p)**2 +
        ((mu_std_xy[1] - mu_std_p[1]) / std_p)**2 +
        ((mu_std_xy[2] - mu_std_p[2]) / std_p)**2 +
        ((mu_std_xy[3] - mu_std_p[3]) / std_p)**2)


def lnlike(mu_std_xy, x, ex, y, ey, mp, rho):
    """
    Log likelihood, product of Gaussian functions.

    Source:
    https://en.wikipedia.org/wiki/
    Multivariate_normal_distribution#Bivariate_case
    """
    sigma_x = np.sqrt(ex**2 + mu_std_xy[2]**2)
    sigma_y = np.sqrt(ey**2 + mu_std_xy[3]**2)

    A = -np.sum(np.log(2. * np.pi * sigma_x * sigma_y * np.sqrt(1. - rho**2)))
    diff_x, diff_y = x - mu_std_xy[0], y - mu_std_xy[1]
    B_x = (diff_x / sigma_x)**2
    B_y = (diff_y / sigma_y)**2
    B_cov = -2. * rho * (diff_x * diff_y) / (sigma_x * sigma_y)
    C = -(1 / (2. * (1 - rho**2))) * np.sum(mp * (B_x + B_y + B_cov))

    return A + C

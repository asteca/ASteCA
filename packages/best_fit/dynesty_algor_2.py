
import numpy as np
from ..synth_clust import synth_cluster
from . import likelihood
from .bf_common import varPars
from dynesty import DynamicNestedSampler
from dynesty import utils as dyfunc
from dynesty import plotting as dyplot
import time as tm
import matplotlib.pyplot as plt
# from .histodd import histogramdd
# from scipy.stats import wasserstein_distance, energy_distance, anderson_ksamp


"""
Trying to replicate the method described in the MiMO article
"""


def main(
    completeness, err_lst, max_mag_syn, obs_clust, lkl_method,
    pt_ntemps, pt_adapt, pt_tmax, nsteps_mcee, nwalkers_mcee, mins_max,
    priors_mcee, ext_coefs, binar_flag, mean_bin_mr, N_fc, m_ini_idx,
    theor_tracks, err_norm_rand, binar_probs, ext_unif_rand, fundam_params,
        st_dist_mass, **kwargs):
    """
    """

    varIdxs, ndim, ranges = varPars(fundam_params)
    # Pack synthetic cluster arguments.
    synthcl_args = [
        completeness, err_lst, max_mag_syn, ext_coefs, binar_flag,
        mean_bin_mr, N_fc, m_ini_idx, st_dist_mass, theor_tracks,
        err_norm_rand, binar_probs, ext_unif_rand]

    #
    def loglike(model):
        # Generate synthetic cluster.
        synth_clust = synth_cluster.main(
            fundam_params, varIdxs, model, *synthcl_args, transpose_flag=False)
        res = integLkl(synth_clust, obs_clust)
        # print(model, res)
        # breakpoint()
        return res

    deltas = ranges[:, 1][:-1] - ranges[:, 0][:-1]

    # Define our uniform prior via the prior transform.
    def ptform(u):
        return np.array(u) * deltas + ranges[:, 0][:-1]

    start = tm.time()

    #
    dsampler = DynamicNestedSampler(loglike, ptform, ndim)
    dsampler.run_nested(maxcall=20000)
    print("\nTime", tm.time() - start)
    print("ESS", dsampler.n_effective)
    res = dsampler.results
    # initialize our nested sampler
    # sampler = NestedSampler(loglike, ptform, ndim)
    # sampler.run_nested(maxiter=20000, maxcall=150000)
    # res = sampler.results

    fig, axes = dyplot.traceplot(
        res, show_titles=True, trace_cmap='viridis', connect=True)
    plt.show()

    samples, weights = res.samples, np.exp(res.logwt - res.logz[-1])
    mean, cov = dyfunc.mean_and_cov(samples, weights)
    print("\n", mean)

    breakpoint()

    return


def integLkl(synth, obs, alpha=-2.):
    """
    obs:   (mag, e_mag, color, e_color)
    synth: (mag, color, mass)
    """
    from scipy.integrate import dblquad

    if synth.shape[1] < 10:
        return -np.inf

    # synth_clust = [mag, c1, mag_b, c1_b, m_ini_1, m_ini_1]
    smag, scol = synth[:2]
    mass = synth[-2]
    mag, col, e_mag, e_col = obs

    dx = (smag[1:] - smag[:-1]) * (scol[1:] - scol[:-1])
    IMF_m = mass[:-1]**alpha

    def deltaIMF(m, c, eps=0.05):
        d = (smag - m)**2 + (scol - c)**2
        if min(d) > eps:
            return 0.

        i = np.argmin(d)
        return mass[i]**alpha

    mmin, mmax, cmin, cmax = smag.min(), smag.max(), scol.min(), scol.max()

    def func(c, m, mi, ci, emi, eci):
        return deltaIMF(m, c) * np.exp(-.5 * ((mi - m) / emi)**2)\
            * np.exp(-.5 * ((ci - c) / eci)**2)

    Phi = np.ones(len(mag))
    for i, magi in enumerate(mag):
        Phi0 = dblquad(func, mmin, mmax, cmin, cmax, args=(
            magi, col[i], e_mag[i], e_col[i]))[0]

        y = IMF_m * np.exp(-.5 * ((magi - smag[:-1]) / e_mag[i])**2)\
            * np.exp(-.5 * ((col[i] - scol[:-1]) / e_col[i])**2)
        Phi[i] = (dx * y).sum()

        print(Phi0, Phi[i])
    breakpoint()

    logLkl = np.log(Phi + 1).sum()

    return logLkl

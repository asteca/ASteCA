
import numpy as np
import datetime as dt
from ..synth_clust import synth_cluster
from . import likelihood
from .bf_common import initPop, varPars, rangeCheck, fillParams
from dynesty import NestedSampler
from dynesty import DynamicNestedSampler
from dynesty import utils as dyfunc
from dynesty import plotting as dyplot
import time as tm
import matplotlib.pyplot as plt
from .histodd import histogramdd
from scipy.stats import wasserstein_distance, energy_distance, anderson_ksamp


"""
See also: https://stats.stackexchange.com/q/404775/10416
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
    # For the other distances
    bin_edges, cl_histo_f_z, cl_z_idx = obs_clust
    dist_pars = (lkl_method, bin_edges, cl_z_idx)

    def loglike(model):
        # Generate synthetic cluster.
        synth_clust = synth_cluster.main(
            fundam_params, varIdxs, model, *synthcl_args)
        return statDist(synth_clust, cl_histo_f_z, dist_pars)

    deltas = ranges[:, 1][:-1] - ranges[:, 0][:-1]

    # Define our uniform prior via the prior transform.
    def ptform(u):
        return np.array(u) * deltas + ranges[:, 0][:-1]

    start = tm.time()

    #
    dsampler = DynamicNestedSampler(loglike, ptform, ndim)
    dsampler.run_nested(maxcall=200000)
    print("\nTime", tm.time() - start)
    print("ESS", dsampler.n_effective)
    res = dsampler.results
    # initialize our nested sampler
    # sampler = NestedSampler(loglike, ptform, ndim)
    # sampler.run_nested(maxiter=20000, maxcall=150000)
    # res = sampler.results

    fig, axes = dyplot.traceplot(
        res, show_titles=True, trace_cmap='viridis', connect=True)

                                 # connect_highlight=range(5))
    plt.show()

    samples, weights = res.samples, np.exp(res.logwt - res.logz[-1])
    mean, cov = dyfunc.mean_and_cov(samples, weights)
    print("\n", mean)

    breakpoint()


    isoch_fit_params = {
        'varIdxs': varIdxs, 'ndim': ndim, 'Tmax': str(Tmax),
        'map_sol': map_sol, 'map_lkl': map_lkl, 'map_lkl_final': map_lkl_final,
        'prob_mean': prob_mean, 'bf_elapsed': elapsed, 'maf_allT': maf_allT,
        'tswaps_afs': tswaps_afs, 'betas_pt': betas_pt, 'N_steps': N_steps,
        'cold_chain': cold_chain
    }

    return isoch_fit_params


def statDist(synth, obs, dist_pars):
    """
    """
    lkl_method, bin_edges, cl_z_idx = dist_pars

    # The likelihood must be maximized so we invert (to minimize) and add a
    # constant to make the result positive
    # return 500 - 
    return likelihood.main(lkl_method, synth, [bin_edges, obs, cl_z_idx])

    # Histogram of the synthetic cluster, using the bin edges calculated
    # with the observed cluster.
    syn_histo = histogramdd(synth, bins=bin_edges)

    # Flatten N-dimensional histogram.
    syn_histo_f = syn_histo.ravel()
    # Remove all bins where n_i = 0 (no observed stars).
    syn_histo_f_z = syn_histo_f[cl_z_idx]

    # return wasserstein_distance(syn_histo_f_z, obs)
    # return energy_distance(syn_histo_f_z, obs)
    return anderson_ksamp([syn_histo_f_z, obs])[0]


# def loglkl(
#     model, fundam_params, lkl_method, obs_clust, ranges, varIdxs, priors,
#         synthcl_args):
#     """
#     """
#     rangeFlag = rangeCheck(model, ranges, varIdxs)

#     logpost = -1e9
#     if rangeFlag:
#         # Generate synthetic cluster.
#         synth_clust = synth_cluster.main(
#             fundam_params, varIdxs, model, *synthcl_args)

#         # Call likelihood function for this model.
#         lkl = likelihood.main(lkl_method, synth_clust, obs_clust)
#         log_p = 0.
#         for i, pr in enumerate(priors):
#             # Gaussian prior. If the prior is uniform, simply pass.
#             if pr[0] == 'g':
#                 log_p += np.log(1 / pr[2]) - .5 * np.square(
#                     (model[i] - pr[1]) / pr[2])

#         # The negative likelihood is returned since Dolphin requires a
#         # minimization of the PLR. Here we are maximizing, hence the minus.
#         logpost = log_p + (-lkl)

#     return logpost

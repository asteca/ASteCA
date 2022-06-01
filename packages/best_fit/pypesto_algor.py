
import numpy as np
from ..synth_clust import synth_cluster
from . import likelihood
from .bf_common import initPop, varPars, rangeCheck, fillParams
from .histodd import histogramdd
from scipy.stats import wasserstein_distance, energy_distance, anderson_ksamp
import pypesto
import pypesto.sample as sample
import pypesto.optimize as optimize
import pypesto.visualize as visualize
import matplotlib.pyplot as plt


def main(
    completeness, err_lst, max_mag_syn, obs_clust, lkl_method,
    pt_ntemps, pt_adapt, pt_tmax, nsteps_mcee, nwalkers_mcee, mins_max,
    priors_mcee, ext_coefs, binar_flag, mean_bin_mr, N_fc, m_ini_idx,
    theor_tracks, err_norm_rand, binar_probs, ext_unif_rand, fundam_params,
        st_dist_mass, **kwargs):
    """
    https://github.com/ICB-DCM/pyPESTO
    """

    varIdxs, ndim, ranges = varPars(fundam_params)
    # Pack synthetic cluster arguments.
    synthcl_args = [
        completeness, err_lst, max_mag_syn, ext_coefs, binar_flag,
        mean_bin_mr, N_fc, m_ini_idx, st_dist_mass, theor_tracks,
        err_norm_rand, binar_probs, ext_unif_rand]

    bin_edges, cl_histo_f_z, cl_z_idx = obs_clust
    obs_data = cl_histo_f_z
    dist_pars = (lkl_method, bin_edges, cl_z_idx)

    def distance(model):
        synth_clust = synth_cluster.main(
            fundam_params, varIdxs, model, *synthcl_args)
        return -statDist(synth_clust, obs_data, dist_pars)

    # first type of objective
    objective = pypesto.Objective(fun=distance)

    dim_full = 6
    lb = ranges[:-1][:, 0].reshape([dim_full, 1])
    ub = ranges[:-1][:, 1].reshape([dim_full, 1])

    problem = pypesto.Problem(objective=objective, lb=lb, ub=ub)

    res_opt = optimize.minimize(
        problem, n_starts=100, progress_bar=False, filename=None)
    init_model = res_opt.optimize_result.list[0].x
    # visualize.waterfall(res_opt)
    # plt.show()
    # breakpoint()

    sampler = sample.AdaptiveParallelTemperingSampler(
        internal_sampler=sample.AdaptiveMetropolisSampler(), n_chains=20)
    result = sample.sample(problem, 1.5e4, sampler, x0=init_model)

    sample.geweke_test(result=result)
    sample.effective_sample_size(result=result)

    breakpoint()
    visualize.sampling_fval_traces(result)
    visualize.sampling_parameter_traces(result, full_trace=True)
    alpha = [99, 95, 90]
    visualize.sampling_parameter_cis(result, alpha=alpha)
    visualize.sampling_scatter(result)
    # i_chain=0?
    visualize.sampling_1d_marginals(result)

    aa = result.optimize_result.list[0]

    breakpoint()


def statDist(synth, obs, dist_pars):
    """
    """
    lkl_method, bin_edges, cl_z_idx = dist_pars

    # bin_edges, cl_histo_f_z, cl_z_idx = obs_clust
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

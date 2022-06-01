
import numpy as np
import datetime as dt
from ..synth_clust import synth_cluster
from . import likelihood
from .bf_common import varPars, rangeCheck, fillParams, random_population
import kombine
import time as tm
import matplotlib.pyplot as plt
from .histodd import histogramdd
from scipy.stats import wasserstein_distance, energy_distance, anderson_ksamp


"""
NOT working until this issue is fixed:
https://github.com/bfarr/kombine/issues/29
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

    #
    nwalkers = 100
    sampler = kombine.Sampler(nwalkers, ndim, loglike, args=(
        fundam_params, varIdxs, synthcl_args, cl_histo_f_z, dist_pars))

    p0 = np.random.uniform(-10, 10, size=(nwalkers, ndim))
    breakpoint()

    p0 = random_population(fundam_params, varIdxs, nwalkers)
    p, post, q = sampler.burnin(p0)

    breakpoint()


def loglike(
        model, fundam_params, varIdxs, synthcl_args, cl_histo_f_z, dist_pars):
    # print("here")
    # breakpoint()
    # Generate synthetic cluster.
    synth_clust = synth_cluster.main(
        fundam_params, varIdxs, model, *synthcl_args)
    return statDist(synth_clust, cl_histo_f_z, dist_pars)


def statDist(synth, obs, dist_pars):
    """
    """
    lkl_method, bin_edges, cl_z_idx = dist_pars

    # The likelihood must be maximized so we invert (to minimize) and add a
    # constant to make the result positive
    return likelihood.main(lkl_method, synth, [bin_edges, obs, cl_z_idx])

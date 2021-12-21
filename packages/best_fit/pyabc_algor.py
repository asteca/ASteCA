
import numpy as np
import datetime as dt
from ..synth_clust import synth_cluster
from . import likelihood
from .bf_common import initPop, varPars, rangeCheck, fillParams
import pyabc
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

    import os
    import tempfile

    varIdxs, ndim, ranges = varPars(fundam_params)
    # Pack synthetic cluster arguments.
    synthcl_args = [
        completeness, err_lst, max_mag_syn, ext_coefs, binar_flag,
        mean_bin_mr, N_fc, m_ini_idx, st_dist_mass, theor_tracks,
        err_norm_rand, binar_probs, ext_unif_rand]

    prior = pyabc.Distribution(
        z=pyabc.RV("uniform", ranges[0][0], ranges[0][1] - ranges[0][0]),
        a=pyabc.RV("uniform", ranges[1][0], ranges[1][1] - ranges[1][0]),
        e=pyabc.RV("uniform", ranges[2][0], ranges[2][1] - ranges[2][0]),
        d=pyabc.RV("uniform", ranges[3][0], ranges[3][1] - ranges[3][0]),
        m=pyabc.RV("uniform", ranges[4][0], ranges[4][1] - ranges[4][0]),
        b=pyabc.RV("uniform", ranges[5][0], ranges[5][1] - ranges[5][0])
    )

    #
    # For the other distances
    bin_edges, cl_histo_f_z, cl_z_idx = obs_clust
    obs_data = {'y': cl_histo_f_z}
    # lkl_method = ''
    dist_pars = (lkl_method, bin_edges, cl_z_idx)

    def model(parameter):
        model_p = (
            parameter['z'], parameter['a'], parameter['e'],
            parameter['d'], parameter['m'], parameter['b'])
        # Generate synthetic cluster.
        synth_clust = synth_cluster.main(
            fundam_params, varIdxs, model_p, *synthcl_args)
        return {"y": synth_clust}

    def distance(synth, obs):
        return statDist(synth['y'], obs['y'], dist_pars)

    abc = pyabc.ABCSMC(
        model, prior, distance, population_size=25,
        sampler=pyabc.sampler.MulticoreEvalParallelSampler(n_procs=1))
    db_path = ("sqlite:///" + os.path.join(tempfile.gettempdir(), "pyABC.db"))
    print(db_path)
    abc.new(db_path, obs_data)

    history = abc.run(
        max_nr_populations=50,  # min_eps_diff=0.001, 
        max_walltime=dt.timedelta(hours=0, minutes=15))

    import matplotlib.pyplot as plt
    max_r = history.max_t

    def makePlot(par):
        fig, ax = plt.subplots()
        for t in range(max_r - 3, max_r + 1):
            df, w = history.get_distribution(m=0, t=t)
            pyabc.visualization.plot_kde_1d(
                df, w,
                x=par, ax=ax,
                label="PDF t={}".format(t))
        ax.legend()
        plt.show()

    makePlot("z")
    breakpoint()

    pyabc.visualization.plot_effective_sample_sizes(history)
    pyabc.visualization.plot_acceptance_rates_trajectory(history)
    pyabc.visualization.plot_eps_walltime(history)
    pyabc.visualization.plot_epsilons(history)
    # m (int, optional (default = 0)) – Id of the model to plot for.
    # t (int, optional (default = None, i.e. the last time)) – Time point to plot for.
    pyabc.visualization.plot_histogram_1d(history, m=0, t=25, x='z')

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
    return 500 - likelihood.main(lkl_method, synth, [bin_edges, obs, cl_z_idx])

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

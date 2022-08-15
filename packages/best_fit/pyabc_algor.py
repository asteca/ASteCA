
import numpy as np
import datetime as dt
from . import likelihood
from .bf_common import getSynthClust
import pyabc
from .histodd import histogramdd
from scipy.stats import wasserstein_distance, energy_distance, anderson_ksamp


"""
Works fine but the uncertainties appear to be larger than those reported
by ptemcee.
I'm not sure how to extract the samples properly to estimate the summary
statistics (mean, median, stddev, etc), I think I need to apply the weights
but I'm not sure how.

See: https://stats.stackexchange.com/q/555473/10416
See also: https://stats.stackexchange.com/q/404775/10416
"""


def main(lkl_method, obs_clust, varIdxs, ndim, ranges, syntClustArgs):
    """
    """

    import os
    import tempfile

    # parameters = np.array(['z', 'a', 'B', 'Av', 'dr', 'Rv', 'dm'])
    rvar = ranges[varIdxs]

    prior = pyabc.Distribution(
        z=pyabc.RV("uniform", rvar[0][0], rvar[0][1] - rvar[0][0]),
        a=pyabc.RV("uniform", rvar[1][0], rvar[1][1] - rvar[1][0]),
        B=pyabc.RV("uniform", rvar[2][0], rvar[2][1] - rvar[2][0]),
        Av=pyabc.RV("uniform", rvar[3][0], rvar[3][1] - rvar[3][0]),
        # dr=pyabc.RV("uniform", rvar[4][0], rvar[4][1] - rvar[4][0]),
        # Rv=pyabc.RV("uniform", rvar[5][0], rvar[5][1] - rvar[5][0]),
        dm=pyabc.RV("uniform", rvar[4][0], rvar[4][1] - rvar[4][0])
    )

    #
    # For the other distances
    bin_edges, cl_histo_f_z, cl_z_idx = obs_clust
    obs_data = {'y': cl_histo_f_z}
    # lkl_method = ''
    dist_pars = (lkl_method, bin_edges, cl_z_idx)

    def model(parameter):
        model_p = (
            parameter['z'], parameter['a'], parameter['B'],
            parameter['Av'], parameter['dm'])
        # Generate synthetic cluster.
        synth_clust = getSynthClust(model_p, True, syntClustArgs)[0]
        return {"y": synth_clust}

    def distance(synth, obs):
        return statDist(synth['y'], obs['y'], dist_pars)

    abc = pyabc.ABCSMC(
        model, prior, distance, population_size=1000,
        sampler=pyabc.sampler.MulticoreEvalParallelSampler(n_procs=2))
    db_path = ("sqlite:///" + os.path.join(tempfile.gettempdir(), "pyABC.db"))
    abc.new(db_path, obs_data)

    history = abc.run(
        max_nr_populations=50,  # min_eps_diff=0.001, 
        max_walltime=dt.timedelta(hours=0, minutes=10))

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

    pyabc.visualization.plot_histogram_matrix(history)
    plt.show()
    breakpoint()

    df, w = history.get_distribution()
    makePlot("z")
    pyabc.visualization.plot_effective_sample_sizes(history)
    pyabc.visualization.plot_acceptance_rates_trajectory(history)
    pyabc.visualization.plot_eps_walltime(history)
    pyabc.visualization.plot_epsilons(history)
    pyabc.visualization.plot_histogram_1d(history, x='z')


def statDist(synth, obs, dist_pars):
    """
    """
    lkl_method, bin_edges, cl_z_idx = dist_pars

    # The likelihood must be maximized so we invert (to minimize) and add a
    # constant to make the result positive
    return 1e6 - likelihood.main(lkl_method, synth, [bin_edges, obs, cl_z_idx])

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

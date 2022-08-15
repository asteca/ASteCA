
import numpy as np
import datetime as dt
from ..synth_clust import synth_cluster
from . import likelihood
from .bf_common import initPop, varPars, rangeCheck, fillParams
from .bf_common import getSynthClust
import pocomc as pc
from .histodd import histogramdd
from scipy.stats import wasserstein_distance, energy_distance, anderson_ksamp


"""
Could not make it work, not sure if the double likelihood is correct
"""


def main(lkl_method, obs_clust, varIdxs, ndim, ranges, synthcl_args, priors):
    """
    """

    # Number of particles to use
    n_particles = 1000

    def log_likelihood(x):
        lkl1 = loglkl(x[0], lkl_method, obs_clust, synthcl_args)
        lkl2 = loglkl(x[1], lkl_method, obs_clust, synthcl_args)
        return np.array([lkl1, lkl2])

    def log_pr(model):
        return logp(model, ranges, varIdxs)

    # Initialise sampler
    sampler = pc.Sampler(n_particles,
                         ndim,
                         log_likelihood=log_likelihood,
                         log_prior=log_pr,
                         bounds=ranges[varIdxs])

    # Initialise particles' positions using samples from the prior (this is
    # very important, other initialisation will not work).
    prior_samples = np.random.uniform(
        low=ranges[varIdxs].T[0], high=ranges[varIdxs].T[1],
        size=(n_particles, ndim))

    # Start sampling
    sampler.run(prior_samples)
    # # We can add more samples at the end
    # sampler.add_samples(1000)
    # Get results
    results = sampler.results

    import matplotlib.pyplot as plt

    breakpoint()


def loglkl(model, lkl_method, obs_clust, synthcl_args):
    """
    """
    # Generate synthetic cluster.
    synth_clust = getSynthClust(model, True, synthcl_args)[0]

    if not synth_clust.any():
        return -np.inf

    # Call likelihood function for this model.
    lkl = likelihood.main(lkl_method, synth_clust, obs_clust)

    return lkl


def logp(model, ranges, varIdxs):
    rangeFlag = rangeCheck(model, ranges, varIdxs)
    if rangeFlag:
        return 0.
    else:
        return -np.inf

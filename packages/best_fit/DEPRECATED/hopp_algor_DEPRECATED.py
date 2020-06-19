
import numpy as np
from .ptemcee_algor import varPars, loglkl, initPop
from .hopp.hoppMCMC import hoppMCMC, chainMCMC


def main(
    lkl_method, e_max, err_lst, completeness, max_mag_syn,
    fundam_params, obs_clust, theor_tracks, R_V,
    ext_coefs, st_dist_mass, N_fc, cmpl_rnd, err_rnd, init_mode_ptm,
    popsize_ptm, maxiter_ptm, ntemps, nwalkers_ptm, nsteps_ptm, nburn_ptm,
        pt_adapt, tmax_ptm, priors_ptm, hmax_ptm):
    """
    It "works" (at least it does not crash), but it seems to get stuck
    at a single point in parameter space. Don't have the time right now to
    inspect what's not working.
    """

    varIdxs, ndim, ranges = varPars(fundam_params)
    # Pack synthetic cluster arguments.
    synthcl_args = [
        theor_tracks, e_max, err_lst, completeness, max_mag_syn, st_dist_mass,
        R_V, ext_coefs, N_fc, cmpl_rnd, err_rnd]

    ntemps, nwalkers_ptm = 1, 1
    # Initial population.
    p0 = initPop(
        ranges, varIdxs, lkl_method, obs_clust, fundam_params, synthcl_args,
        ntemps, nwalkers_ptm, init_mode_ptm, popsize_ptm, maxiter_ptm)[0][0]

    def fitness(model):
        return loglkl(
            model, fundam_params, synthcl_args, lkl_method, obs_clust, ranges,
            varIdxs, priors_ptm)

    # run for 50 hopp-steps
    num_hopp = 10
    # each hopp-step comprises 100 adaptation steps
    num_adapt = 20
    # run with 3 parallel chains
    num_chain = 3
    # each chain is 10 iterations long
    chain_length = 10

    results = hoppMCMC(
        fitness, param=p0, varmat=np.diag([1e-14] * 4), gibbs=True,
        rangeT=[1, 10], model_comp=10, num_hopp=num_hopp, num_adapt=num_adapt,
        num_chain=num_chain, chain_length=chain_length)

    # This will select the final hopp-step and run an MCMC chain to sample
    # from around the identified posterior mode
    n = len(results.parmats) - 1
    parmat = np.array(results.parmats[n])
    param = parmat[parmat[:, 0] == min(parmat[:, 0]), 1:][0]
    varmat = np.cov(parmat[:, 1:], rowvar=False)

    # Parameters to infer by the MCMC chain
    inferpar = list(range(ndim))

    # The following will run the chainMCMC algorithm
    mcmc = chainMCMC(
        fitness, param=param, varmat=varmat, inferpar=inferpar, gibbs=True,
        chain_id=0, pulsevar=1, anneal=1, accthr=0.25, varmat_change=0,
        pulse_change=10, pulse_change_ratio=2, print_iter=100)

    # This will iterate the chainMCMC for 1000 steps recording every 25th
    # iteration while discarding the initial 500 steps
    samples = []
    for m in range(1000):
        mcmc.iterate()
        if m > 100:  # and m%25==0:
            samples.append(mcmc.getParam())
    samples = np.array(samples)

    import pdb; pdb.set_trace()  # breakpoint 7374b507 //

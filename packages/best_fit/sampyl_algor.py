
import sampyl as smp
from sampyl import np
# import numpy as np
from .emcee_algor import varPars, log_posterior, convergenceVals


def main(
    lkl_method, e_max, err_lst, completeness, max_mag_syn,
        fundam_params, obs_clust, theor_tracks, R_V, ext_coefs, st_dist_mass,
        N_fc, cmpl_rnd, err_rnd, nsteps, nburn):
    """
    """
    varIdxs, ndim, ranges = varPars(fundam_params)

    synthcl_args = [
        theor_tracks, e_max, err_lst, completeness, max_mag_syn, st_dist_mass,
        R_V, ext_coefs, N_fc, cmpl_rnd, err_rnd]

    def logp(met, age, ext, dist, mass, binar):
        # met, age = 0.0305, 9.8
        model_done = {}
        model = [met, age, ext, dist, mass, binar]
        lp = log_posterior(
            model, varIdxs, ranges, fundam_params, synthcl_args, lkl_method,
            obs_clust, model_done)
        # The log_posterior returns a value (log_prior+log_posterior) that
        # needs to be maximized.
        return lp

    # bounds = {
    #     'met': ranges[0], 'age': ranges[1], 'ext': ranges[2],
    #     'dist': ranges[3], 'mass': ranges[4], 'binar': ranges[5]}
    # start = smp.find_MAP(logp, {
    #     'met': .015, 'age': 9., 'ext': 0., 'dist': 13., 'mass': 5000.,
    #     'binar': .3}, verbose=True, bounds=bounds)

    start = {
        'met': .03, 'age': 9.7, 'ext': 0.3, 'dist': 12.7, 'mass': 5000.,
        'binar': .3}
    # start = {'ext': 0.3, 'dist': 12.7, 'mass': 5000., 'binar': .3}

    # sampler = smp.NUTS(logp, start)
    sampler = smp.Slice(logp, start)
    chain = sampler.sample(int(nsteps), burn=nburn)

    chains_nruns = np.array([
        chain.met, chain.age, chain.ext, chain.dist, chain.mass,
        chain.binar]).T
    print(chains_nruns.shape)
    mcmc_trace = chains_nruns.T

    import seaborn ; import matplotlib.pyplot as plt
    import pdb; pdb.set_trace()  # breakpoint 01ee2136 //


    # Convergence parameters.
    acorr_t, max_at_10c, geweke_z, emcee_acorf, pymc3_ess, minESS,\
        mESS, mESS_epsilon = convergenceVals(ndim, chains_nruns, mcmc_trace)

    isoch_fit_params = {
        'nsteps': nsteps, 'best_sol': best_sol, 'map_sol': map_sol,
        'pars_chains_bi': pars_chains_bi, 'pars_chains': chains_nruns.T,
        'm_accpt_fr': m_accpt_fr, 'varIdxs': varIdxs,
        'emcee_trace': emcee_trace, 'autocorr_time': acorr_t,
        'max_at_10c': max_at_10c, 'minESS': minESS, 'mESS': mESS,
        'mESS_epsilon': mESS_epsilon, 'emcee_acorf': emcee_acorf,
        'geweke_z': geweke_z, 'pymc3_ess': pymc3_ess
    }




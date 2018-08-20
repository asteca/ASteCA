
import time as t
import sampyl as smp
# from sampyl import np
# from ..core_imp import np
# import numpy as np
from .emcee_algor import varPars, log_posterior, convergenceVals, closeSol


def main(
    lkl_method, e_max, err_lst, completeness, max_mag_syn,
        fundam_params, obs_clust, theor_tracks, R_V, ext_coefs, st_dist_mass,
        N_fc, cmpl_rnd, err_rnd, nwalkers, nsteps, nburn, N_burn,
        emcee_a, priors):
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
            model, priors, varIdxs, ranges, fundam_params, synthcl_args,
            lkl_method, obs_clust, model_done)
        # The log_posterior returns a value (log_prior+log_posterior) that
        # needs to be maximized.
        return lp

    import autograd.numpy as np  # Thinly-wrapped numpy
    from autograd import grad    # The only autograd function you may ever need
    #
    grad_logp = grad(logp)       # Obtain its gradient function
    import pdb; pdb.set_trace()  # breakpoint 544a7e21 //

    grad_logp(.03, 9.7, 0.3, 12.7, 5000., .3)

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

    s = t.time()
    sampler = smp.NUTS(logp, start)
    # sampler = smp.Slice(logp, start)
    chain = sampler.sample(int(nsteps), n_chains=1, burn=nburn)

    # TODO DOES NOT ACCEPT MORE THAN 1 CHAINS
    # (nsteps - nburn, nwalkers, ndim)
    chains_nruns = np.array([
        chain.met, chain.age, chain.ext, chain.dist, chain.mass,
        chain.binar]).T
    chains_nruns = chains_nruns[:, np.newaxis, :]
    print(chains_nruns.shape)
    # (ndim, (nsteps - nburn) * nwalkers)
    mcmc_trace = chains_nruns.reshape(-1, ndim).T
    elapsed = t.time() - s

    N_conv, N_steps_conv, tol_conv, runs = 0., 0., 0., nsteps - nburn
    map_sol, map_lkl, map_lkl_final = [np.nan] * 6, [[0.], [0.]], np.nan
    maf_steps = [[0.], [0.]]
    tau_index, tau_autocorr = 0, np.empty(nsteps)

    # Convergence parameters.
    acorr_t, max_at_5c, min_at_5c, geweke_z, emcee_acorf, pymc3_ess, minESS,\
        mESS, mESS_epsilon, mcmc_halves = convergenceVals(
            ndim, varIdxs, N_conv, chains_nruns, mcmc_trace)
    import pdb; pdb.set_trace()  # breakpoint dd76f2fd //
    

    # Pass the mean as the best model fit found.
    best_sol = closeSol(fundam_params, varIdxs, np.mean(mcmc_trace, axis=1))

    import pdb; pdb.set_trace()  # breakpoint 2d909a3b //

    isoch_fit_params = {
        'varIdxs': varIdxs, 'nsteps': runs, 'best_sol': best_sol,
        'map_sol': map_sol, 'map_lkl': map_lkl, 'map_lkl_final': map_lkl_final,
        'mcmc_elapsed': elapsed, 'mcmc_trace': mcmc_trace,
        'pars_chains_bi': pars_chains_bi, 'pars_chains': chains_nruns.T,
        'maf_steps': maf_steps, 'autocorr_time': acorr_t,
        'max_at_5c': max_at_5c, 'min_at_5c': min_at_5c,
        'minESS': minESS, 'mESS': mESS, 'mESS_epsilon': mESS_epsilon,
        'emcee_acorf': emcee_acorf, 'geweke_z': geweke_z,
        'pymc3_ess': pymc3_ess, 'mcmc_halves': mcmc_halves,
        'N_steps_conv': N_steps_conv, 'N_conv': N_conv, 'tol_conv': tol_conv,
        'tau_index': tau_index, 'tau_autocorr': tau_autocorr
    }

    return isoch_fit_params

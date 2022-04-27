
import numpy as np
import warnings
import time as t
from .bf_common import getSynthClust, initPop, rangeCheck, fillParams
from .ptemcee import sampler
from . import likelihood
from .. import update_progress


def main(
    lkl_method, pt_ntemps, pt_adapt, pt_tmax, nsteps_mcee, nwalkers_mcee,
    mins_max, priors_mcee, fundam_params, obs_clust, varIdxs, ndim, ranges,
        synthcl_args):
    """
    """

    if pt_tmax in ('n', 'none', 'None'):
        Tmax = None
    elif pt_tmax == 'inf':
        Tmax = np.inf
    else:
        Tmax = float(pt_tmax)
    if pt_ntemps in ('n', 'none', 'None'):
        pt_ntemps = None
    else:
        pt_ntemps = int(float(pt_ntemps))

    # Start timing.
    max_secs = mins_max * 60.
    available_secs = max(30, max_secs)

    # Define Parallel tempered sampler
    ptsampler = sampler.Sampler(
        nwalkers_mcee, ndim, loglkl, logp,
        loglargs=[lkl_method, obs_clust, ranges, varIdxs,
                  priors_mcee, synthcl_args], Tmax=Tmax, ntemps=pt_ntemps)

    ntemps = ptsampler.ntemps
    # Initial population.
    pos0 = initPop(
        ranges, varIdxs, lkl_method, obs_clust, fundam_params, synthcl_args,
        ntemps, nwalkers_mcee, 'random', None, None)

    # Track how the acceptance fractions, and temperature swaps acceptance
    # fractions.
    afs, tswaps = [], []
    # Store for Lkl values for plotting.
    prob_mean, map_lkl, map_sol_old = [], [], [[], -np.inf]

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        N_steps_store, runs = 50, 0

        elapsed, start = 0., t.time()
        for i, (pos, lnprob, lnlike) in enumerate(ptsampler.sample(
                pos0, iterations=nsteps_mcee, adapt=pt_adapt)):

            elapsed += t.time() - start
            start = t.time()
            remaining_time = max(0, min(
                available_secs, (nsteps_mcee * elapsed) / (i + 1)) - elapsed)
            m, s = divmod(remaining_time, 60)
            h, m = divmod(m, 60)
            if m > 0:
                tt, hm, ms = "hm", h, m
            else:
                tt, hm, ms = "ms", m, s
            txt = " [{:.0f} models/sec | {:.0f}{}{:.0f}{}]".format(
                (ntemps * nwalkers_mcee * (i + 1)) / elapsed, hm, tt[0], ms,
                tt[1])
            update_progress.updt(nsteps_mcee, i + 1, txt)

            # Only check convergence every 'N_steps_store' steps
            if (i + 1) % N_steps_store:
                continue
            runs += 1

            # Temperature swap acceptance fractions.
            tswaps.append(ptsampler.tswap_acceptance_fraction)
            # Mean acceptance fractions for all temperatures.
            afs.append(np.mean(ptsampler.acceptance_fraction, axis=1))

            # Store MAP solution in this iteration.
            prob_mean.append(np.mean(lnprob[0]))
            idx_best = np.argmax(lnprob[0])
            # Update if a new optimal solution was found.
            if lnprob[0][idx_best] > map_sol_old[1]:
                map_sol_old = [
                    fillParams(fundam_params, varIdxs, pos[0][idx_best]),
                    lnprob[0][idx_best]]
            map_lkl.append(map_sol_old[1])

            # Stop when available time is consumed.
            if elapsed >= available_secs:
                print("  Time consumed")
                break

    # Total number of steps
    N_steps = N_steps_store * np.arange(1, runs + 1)

    # Mean acceptance fractions for all replicas.
    maf_allT = np.asarray(afs).T
    # Temperature swaps acceptance fractions.
    tswaps_afs = np.asarray(tswaps).T
    # Betas history
    betas_pt = ptsampler.beta_history[:, (N_steps_store - 1)::N_steps_store]

    # Final MAP fit.
    map_sol, map_lkl_final = map_sol_old

    # # This number should be between approximately 0.25 and 0.5 if everything
    # # went as planned.
    # m_accpt_fr = np.mean(ptsampler.acceptance_fraction[0])
    # if m_accpt_fr > .5 or m_accpt_fr < .25:
    #     print("  WARNING: mean acceptance fraction is outside of the\n"
    #           "  recommended range.")

    # ptsampler.chain.shape: (ntemps, nchains, nsteps, ndim)
    # cold_chain.shape: (i, nchains, ndim)
    cold_chain = ptsampler.chain[0, :, :i, :].transpose(1, 0, 2)

    isoch_fit_params = {
        'Tmax': str(Tmax), 'map_sol': map_sol, 'map_lkl': map_lkl,
        'map_lkl_final': map_lkl_final,
        'prob_mean': prob_mean, 'bf_elapsed': elapsed, 'maf_allT': maf_allT,
        'tswaps_afs': tswaps_afs, 'betas_pt': betas_pt, 'N_steps': N_steps,
        'cold_chain': cold_chain
    }

    return isoch_fit_params


def loglkl(
        model, lkl_method, obs_clust, ranges, varIdxs, priors, synthcl_args):
    """
    """
    rangeFlag = rangeCheck(model, ranges, varIdxs)

    logpost = -1e9
    if rangeFlag:
        # Generate synthetic cluster.
        synth_clust = getSynthClust(model, *synthcl_args, True)[0]

        # Call likelihood function for this model.
        lkl = likelihood.main(lkl_method, synth_clust, obs_clust)
        log_p = 0.
        for i, pr in enumerate(priors):
            # Gaussian prior. If the prior is uniform, simply pass.
            if pr[0] == 'g':
                log_p += np.log(1 / pr[2]) - .5 * np.square(
                    (model[i] - pr[1]) / pr[2])

        logpost = log_p + lkl

    return logpost


def logp(_):
    """
    Just here as a place holder for 'ptemcee'.

    The prior is moved inside the log-likelihood to save computation time
    because 'ptemcee' calls both functions and then adds the result. But if
    the model is outside of the allowed range, there is no point in checking
    twice (ie: calling rangeCheck() twice). As I can not "send" information
    from here to the likelihood (as can be done in 'emcee'), it's better to
    just put everything inside the loglkl() function.
    """
    return 0.

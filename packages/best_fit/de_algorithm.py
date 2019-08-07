
import time as t
import numpy as np
from scipy.optimize import differential_evolution as DE
import warnings
from .bf_common import synthClust, varPars, fillParams
from . import likelihood
from .. import update_progress


def main(
    e_max, err_lst, completeness, max_mag_syn, st_dist_mass, ext_coefs, N_fc,
    cmpl_rnd, err_rnd, obs_clust, popsize, maxiter, available_secs,
    flag_print_perc, fundam_params, lkl_method, theor_tracks, R_V,
        **kwargs):
    """
    """

    # Start timing.
    start_t = t.time()

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        varIdxs, ndim, ranges = varPars(fundam_params)

        # Pack synthetic cluster arguments.
        synthcl_args = [
            theor_tracks, e_max, err_lst, completeness, max_mag_syn,
            st_dist_mass, R_V, ext_coefs, N_fc, cmpl_rnd, err_rnd]

        def get_callback(start_t):
            def callbackF(current_params, convergence):
                time_elapsed = t.time() - start_t
                if time_elapsed > available_secs:
                    return True
            return callbackF

        def DEdist(model, info):
            # Print progressbar
            if flag_print_perc and\
                    info['Nfeval'] < popsize * (maxiter + 1) * len(model):
                update_progress.updt(
                    popsize * (maxiter + 1) * len(model), info['Nfeval'] + 1)
                info['Nfeval'] += 1

            synth_clust = synthClust(
                fundam_params, varIdxs, model, synthcl_args)
            if synth_clust:
                lkl = likelihood.main(lkl_method, synth_clust, obs_clust)
                return lkl
            return np.inf

        result = DE(
            DEdist, ranges[varIdxs], popsize=popsize, maxiter=maxiter,
            disp=False, args=({'Nfeval': 0},), callback=get_callback(start_t))

    # If this is a bootstrap run, return the best model found only.
    if not flag_print_perc:
        return fillParams(fundam_params, varIdxs, result.x)

    elapsed = t.time() - start_t
    if elapsed > available_secs:
        print("\n  Time consumed (steps={})".format(result.nit + 1))

    map_sol_filled = fillParams(fundam_params, varIdxs, result.x)
    map_sol = []
    for i, par in enumerate(fundam_params):
        map_sol.append(min(par, key=lambda x: abs(x - map_sol_filled[i])))

    # Store the ML solution as 'map_sol' for consistency.
    isoch_fit_params = {
        'OF_elapsed': elapsed, 'map_sol': map_sol, 'OF_steps': result.nit,
        'OF_models': popsize * (maxiter + 1) * len(result.x)}

    return isoch_fit_params


import time as t
import numpy as np
import random
from . import genetic_algorithm, de_algorithm
from . import obs_clust_prepare
from .. import update_progress
from .bf_common import varPars, modeKDE, fillParams, r2Dist


def main(
    pd, clp, cl_max_mag, max_mag_syn, obs_clust, ext_coefs, st_dist_mass,
        N_fc, cmpl_rnd, err_rnd):
    '''
    Non-parametric bootstrap process, runs the selected algorithm a number of
    times each time generating a new observed cluster representation through
    resampling with replacement.
    '''

    start_t = t.time()
    max_secs = pd['hmax'] * 60. * 60.

    if pd['hperc_btstrp'] == 0.:
        # If bootstrap is set not to run, use all the time.
        available_secs = max_secs
    else:
        available_secs = max_secs * (1. - pd['hperc_btstrp'])

    # # TODO not yet fully implemented
    # pd['best_fit_algor'] = 'boot+DE'

    if pd['best_fit_algor'] == 'boot+GA':
        # Use random initial population
        init_pop, flag_print_perc = [], True
        # Initial number of steps for the numerical optimizer applied over
        # the observed data.
        N_pop_init, N_gen_init = pd['N_pop'], pd['N_gen']
        argsOF = [
            pd, clp, max_mag_syn, obs_clust, ext_coefs, st_dist_mass, N_fc,
            cmpl_rnd, err_rnd, available_secs, init_pop, N_pop_init,
            N_gen_init, flag_print_perc]

    elif pd['best_fit_algor'] == 'boot+DE':
        popsize, maxiter = pd['N_pop'], pd['N_gen']

        argsOF = [
            pd, clp, max_mag_syn, st_dist_mass, ext_coefs, N_fc, cmpl_rnd,
            err_rnd, obs_clust, popsize, maxiter, available_secs, True]

    # First run of the numerical optimizer function with the observed data.
    isoch_fit_params = optimizerFunc(pd['best_fit_algor'], argsOF)

    # Check time available for bootstrap after the previous run.
    available_secs = max_secs - isoch_fit_params['OF_elapsed']

    # Holds the parameter values obtained by the bootstrap process.
    params_boot, N_steps_availb, btstrp_start_t = [], 0, t.time()
    if pd['hperc_btstrp'] > 0.:

        print('Begin bootstrap process ({} | {}).'.format(
            pd['N_pop_btstrp'], pd['N_step_btstrp']))

        btstrp_init, btstrp_t, btstrp_runs = t.time(), 0., 0
        while btstrp_t < available_secs:

            # Resample cluster with replacement.
            cl_max_mag_ran = resample_replacement(cl_max_mag)
            # Prepare observed cluster according to the likelihood
            # method selected.
            obs_cl = obs_clust_prepare.main(
                cl_max_mag_ran, pd['lkl_method'], pd['lkl_binning'],
                pd['lkl_weight'])

            if pd['best_fit_algor'] == 'boot+GA':
                argsOF = [
                    pd, clp, max_mag_syn, obs_cl, ext_coefs, st_dist_mass,
                    N_fc, cmpl_rnd, err_rnd, np.inf, isoch_fit_params[
                        'OF_final_generation'][:pd['N_pop_btstrp']],
                    pd['N_pop_btstrp'], pd['N_step_btstrp'], False]

            elif pd['best_fit_algor'] == 'boot+DE':
                argsOF = [
                    pd, clp, max_mag_syn, st_dist_mass, ext_coefs, N_fc,
                    cmpl_rnd, err_rnd, obs_cl, pd['N_pop_btstrp'],
                    pd['N_step_btstrp'], np.nan, False]

            params_boot.append(optimizerFunc(pd['best_fit_algor'], argsOF))

            btstrp_runs += 1
            btstrp_t = t.time() - btstrp_init
            update_progress.updt(available_secs, btstrp_t)
            if btstrp_t > available_secs:
                break

        # As array with (params, btstrp_runs) shape
        isoch_fit_params['params_boot'] = np.array(params_boot).T
        isoch_fit_params['varIdxs'] = varPars(pd['fundam_params'])[0]
        # Remove fixed parameters from array.
        isoch_fit_params['params_boot'] = isoch_fit_params['params_boot'][
            isoch_fit_params['varIdxs']]

        # Used by cornerPlot()
        isoch_fit_params['param_r2'] = r2Dist(
            pd['fundam_params'], isoch_fit_params['varIdxs'],
            isoch_fit_params['params_boot'])

        # Point estimates (mean, median, mode).
        isoch_fit_params = getSols(
            pd, isoch_fit_params['params_boot'], isoch_fit_params,
            isoch_fit_params['varIdxs'])

        print("Total bootstrap runs: {}".format(btstrp_runs))
        isoch_fit_params['N_bootstrap'] = btstrp_runs

    else:
        print('Skip bootstrap process.')
        isoch_fit_params['N_bootstrap'] = 0
        isoch_fit_params['params_boot'] = np.array(params_boot).T
        isoch_fit_params['mean_sol'] = [np.nan] * 6
        isoch_fit_params['median_sol'] = [np.nan] * 6
        isoch_fit_params['mode_sol'] = [np.nan] * 6
        isoch_fit_params['param_r2'] = [np.nan] * 6

    isoch_fit_params['btstrp_t'] = t.time() - btstrp_start_t
    isoch_fit_params['bf_elapsed'] = t.time() - start_t
    isoch_fit_params['N_total'] = float(
        isoch_fit_params['OF_steps'] + int(N_steps_availb))

    return isoch_fit_params


def optimizerFunc(best_fit_algor, args):
    """
    TODO fully implement the DE and other methods
    """
    if best_fit_algor == 'boot+GA':
        pd, clp, max_mag_syn, obs_clust, ext_coefs, st_dist_mass, N_fc,\
            cmpl_rnd, err_rnd, available_secs, init_pop, N_pop, N_gen,\
            flag_print_perc = args

        if flag_print_perc:
            # Print advances.
            isoch_fit_params = genetic_algorithm.main(
                available_secs, init_pop, flag_print_perc, N_pop, N_gen,
                max_mag_syn, obs_clust, ext_coefs, st_dist_mass, N_fc,
                cmpl_rnd, err_rnd, clp['em_float'], clp['err_lst'],
                clp['completeness'], **pd)
        else:
            isoch_fit_params = genetic_algorithm.main(
                available_secs, init_pop, flag_print_perc, N_pop, N_gen,
                max_mag_syn,
                obs_clust, ext_coefs, st_dist_mass, N_fc, cmpl_rnd, err_rnd,
                clp['em_float'], clp['err_lst'], clp['completeness'],
                **pd)

    elif best_fit_algor == 'boot+DE':
        pd, clp, max_mag_syn, st_dist_mass, ext_coefs, N_fc, cmpl_rnd,\
            err_rnd, obs_clust, popsize, maxiter, available_secs,\
            flag_print_perc = args

        if flag_print_perc:
            # Print advances.
            isoch_fit_params = de_algorithm.main(
                clp['em_float'], clp['err_lst'], clp['completeness'],
                max_mag_syn, st_dist_mass, ext_coefs, N_fc, cmpl_rnd, err_rnd,
                obs_clust, popsize, maxiter, available_secs,
                flag_print_perc, **pd)
        else:
            isoch_fit_params = de_algorithm.main(
                clp['em_float'], clp['err_lst'], clp['completeness'],
                max_mag_syn, st_dist_mass, ext_coefs, N_fc, cmpl_rnd, err_rnd,
                obs_clust, popsize, maxiter, available_secs,
                flag_print_perc, **pd)

    return isoch_fit_params


def resample_replacement(cl_max_mag):
    '''
    Resamples the observed cluster with replacement.
    '''
    cl_max_mag_ran = [random.choice(cl_max_mag) for _ in cl_max_mag]

    return cl_max_mag_ran


def getSols(pd, params_boot, isoch_fit_params, varIdxs):
    """
    Save mean, median, and mode solutions for the bootstrap.
    """
    # 'pardist_kde' is used for cornerPlot()
    mode_sol, isoch_fit_params['pardist_kde'] = modeKDE(
        pd['fundam_params'], varIdxs, params_boot)
    # 'mode_sol' comes back with the free parameters only.
    mode_sol = fillParams(pd['fundam_params'], varIdxs, mode_sol)

    # Push values to the closest values in the grid.
    mean_boot_sol, median_boot_sol, mode_boot_sol = [], [], []
    j = 0
    for i, par in enumerate(pd['fundam_params']):
        if i in isoch_fit_params['varIdxs']:
            mean_boot_sol.append(min(
                par, key=lambda x: abs(x - np.mean(params_boot[i - j]))))
            median_boot_sol.append(min(
                par, key=lambda x: abs(x - np.median(params_boot[i - j]))))
            mode_boot_sol.append(min(
                par, key=lambda x: abs(x - mode_sol[i - j])))
        else:
            mean_boot_sol.append(par[0])
            median_boot_sol.append(par[0])
            mode_boot_sol.append(par[0])
            j += 1

    isoch_fit_params['mean_sol'] = mean_boot_sol
    isoch_fit_params['median_sol'] = median_boot_sol
    isoch_fit_params['mode_sol'] = mode_boot_sol

    return isoch_fit_params

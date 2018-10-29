
import numpy as np
import random
from . import genetic_algorithm
from . import obs_clust_prepare
from .. import update_progress


def resample_replacement(cl_max_mag):
    '''
    Resamples the observed cluster with replacement. Used by the bootstrap
    process.
    '''
    cl_max_mag_ran = [random.choice(cl_max_mag) for _ in cl_max_mag]

    return cl_max_mag_ran


def main(lkl_method, e_max, err_lst, completeness, fundam_params,
         cl_max_mag, max_mag_syn, theor_tracks, R_V, ext_coefs, st_dist_mass,
         N_fc, cmpl_rnd, err_rnd, N_pop, N_gen, fit_diff, cross_prob,
         cross_sel, mut_prob, N_el, N_ei, N_es, lkl_binning, lkl_weight,
         N_b, flag_print_perc, isoch_fit_params):
    '''
    Non-parametric bootstrap process, runs the selected algorithm a number of
    times each time generating a new observed cluster representation through
    resampling with replacement.
    '''
    if N_b >= 2:
        print('Begin bootstrap process ({}).'.format(N_b))

        # Holds the parameter values obtained by the bootstrap process.
        params_boot = []

        # Begin bootstrap block (run a minimum of two times).
        for i in range(N_b):

            # Resample cluster with replacement.
            cl_max_mag_ran = resample_replacement(cl_max_mag)
            # Obtain prepared observed cluster according to the likelihood
            # method selected.
            obs_cl = obs_clust_prepare.main(
                cl_max_mag_ran, lkl_method, lkl_binning, lkl_weight)

            params_boot.append(genetic_algorithm.main(
                lkl_method, e_max, err_lst, completeness, max_mag_syn,
                fundam_params, obs_cl, theor_tracks, R_V, ext_coefs,
                st_dist_mass, N_fc, cmpl_rnd, err_rnd, N_pop, N_gen, fit_diff,
                cross_prob, cross_sel, mut_prob, N_el, N_ei, N_es,
                flag_print_perc)['map_sol'])

            update_progress.updt(N_b, i + 1)

        mean_boot_sol = closeSol(fundam_params, np.mean(params_boot, 0))
        # Calculate errors for each fundamental parameter.
        isoch_fit_errors = np.std(params_boot, 0)
        for i, p_er in enumerate(isoch_fit_errors):
            # If any parameter has a single valued range, assign 'nan'
            if len(fundam_params[i]) > 1:
                # Find largest delta in this parameter used values.
                largest_delta = np.diff(fundam_params[i]).max()
                # Errors can not be smaller than the largest step in each
                # parameter
                isoch_fit_errors[i] = max(largest_delta, p_er)
            else:
                isoch_fit_errors[i] = np.nan
    else:
        print('Skip bootstrap process.')
        mean_boot_sol = isoch_fit_params['map_sol']
        # No error assignment.
        isoch_fit_errors = [np.nan] * len(isoch_fit_params['map_sol'])

    return isoch_fit_errors, mean_boot_sol


def closeSol(fundam_params, model):
    """
    Find the closest value in the parameters list for the discrete parameters
    metallicity, age, and mass.
    """
    model_proper = []
    for i, par in enumerate(fundam_params):
        # If it is the parameter metallicity, age or mass.
        if i in [0, 1, 4]:
            # Select the closest value in the array of allowed values.
            model_proper.append(min(
                par, key=lambda x: abs(x - model[i])))
        else:
            model_proper.append(model[i])

        # # Select the closest value in the array of allowed values.
        # model_proper.append(min(par, key=lambda x: abs(x - model[i])))

    return model_proper

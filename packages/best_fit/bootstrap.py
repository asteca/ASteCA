
import numpy as np
import random
import genetic_algorithm
import obs_clust_prepare
from .. import update_progress


def resample_replacement(obs_clust):
    '''
    Resamples the observed cluster with replacement. Used by the bootstrap
    process.
    '''
    obs_cl = [random.choice(obs_clust) for _ in obs_clust]

    return obs_cl


def main(lkl_method, e_max, bin_mr, err_lst, completeness, fundam_params,
         cl_max_mag, max_mag_syn, theor_tracks, R_V, ext_coefs, st_dist_mass,
         N_fc, N_pop, N_gen, fit_diff, cross_prob, cross_sel, mut_prob, N_el,
         N_ei, N_es, bin_method, best_fit_algor, N_b):
    '''
    Bootstrap process, runs the selected algorithm a number of times each
    time generating a new observed cluster representation through resampling
    with replacement.
    '''
    print('Begin bootstrap process ({}).'.format(N_b))

    # Holds the parameter values obtained by the bootstrap process.
    params_boot = []

    # Begin bootstrap block (run a minimum of two times).
    for i in range(N_b):

        # Resample cluster with replacement.
        obs_cl_r = resample_replacement(cl_max_mag)
        # Obtain prepared observed cluster according to the likelihood method
        # selected.
        obs_cl = obs_clust_prepare.main(obs_cl_r, lkl_method, bin_method)

        # Algorithm selected.
        if best_fit_algor == 'genet':
            # Let the GA algorithm know this call comes from the bootstrap
            # process so it will not print percentages to screen.
            flag_print_perc = False
            params_boot.append(genetic_algorithm.main(
                lkl_method, e_max, bin_mr, err_lst, completeness, max_mag_syn,
                fundam_params, obs_cl, theor_tracks, R_V, ext_coefs,
                st_dist_mass, N_fc, N_pop, N_gen, fit_diff, cross_prob,
                cross_sel, mut_prob, N_el, N_ei, N_es, flag_print_perc)[0])

        update_progress.updt(N_b, i + 1)

    # Calculate errors for each fundamental parameter.
    isoch_fit_errors = np.std(params_boot, 0)
    for i, p_er in enumerate(isoch_fit_errors):
        # If any parameter has a single valued range, assign an error of -1.
        if len(fundam_params[i]) > 1:
            # Find largest delta in this parameter used values.
            largest_delta = np.diff(fundam_params[i]).max()
            # Errors can not be smaller than the largest step in each parameter
            isoch_fit_errors[i] = max(largest_delta, p_er)
        else:
            isoch_fit_errors[i] = np.nan

    return isoch_fit_errors

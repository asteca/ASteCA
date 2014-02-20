# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 16:10:05 2014

@author: gabriel
"""

import numpy as np
import random

from genetic_algorithm import gen_algor as g_a
from brute_force_algorithm import brute_force as b_f
from synth_cluster import synth_clust as s_c
from get_IMF_PDF import IMF_PDF as i_p
from move_isochrone import move_isoch


def boostrap_resample(stars_list):
    '''
    Resamples the observed cluster with replacement. Used by the bootstrap
    process.
    '''
    obs_clust = [random.choice(stars_list) for _ in stars_list]
    return obs_clust
    


def best_fit(err_lst, memb_prob_avrg_sort, completeness, ip_list, bf_params,
             sc_params, ga_params, ps_params):
    '''
    Perform a best fitting process to find the cluster's parameters:
    E(B-V), distance (distance modulus), age and metallicity.
    '''
    
    bf_flag, best_fit_algor, N_b = bf_params
    
    # Check if algorithm should run.
    if bf_flag:
    
        # Obtain the selected IMF's PDF. We run it once because the array only
        # depends on the IMF selected.
        imf_pdf = i_p(sc_params[1])
        # Replace the name of the IMF by its PDF.
        sc_params[1] = imf_pdf
    
        # Call algorithm to calculate the likelihoods for the set of
        # isochrones and return the best fitting parameters.
        if best_fit_algor == 'brute':
        # Brute force algorithm.
            isoch_fit_params = b_f(err_lst, memb_prob_avrg_sort, completeness,
                                   ip_list, sc_params, ga_params)
            # Assign errors as the steps in each parameter.
            isoch_fit_errors = [ps_params[i+3][2] for i in range(4)]

            
        elif best_fit_algor == 'genet':
            # Genetic algorithm.
            # The first pass is done with no resampling to calculate the final
            # values. After that we resample to get the uncertainty in each
            # parameter.  
            isoch_fit_params = g_a(err_lst, memb_prob_avrg_sort, completeness,
                                   ip_list, sc_params, ga_params)    
        
            # List that holds the parameters values obtained by the bootstrap
            # process.
            params_boot = []
            # Begin bootstrap block (run a minimum of two times).
            for i in range(max(N_b, 2)):
                # Resample cluster with replacement.
                obs_clust = boostrap_resample(memb_prob_avrg_sort)
                # Genetic algorithm.
                params_boot.append(g_a(err_lst, obs_clust, completeness, ip_list,
                                       sc_params, ga_params)[0])
    
            # Calculate errors for each parameter.
            isoch_fit_errors = np.std(params_boot, 0)
            # Errors can not be smaller than the steps in each parameter.
            for i in range(4):
                isoch_fit_errors[i] = max(ps_params[i+3][2], isoch_fit_errors[i])
        
        
        # For plotting purposes.
        # Get list of stored isochrones and their parameters.
        isoch_list, isoch_ma = ip_list[0], ip_list[1]
        # Read best fit values for all parameters.
        m, a, e, d = isoch_fit_params[0]
        # Find indexes for metallixity and age.
        m_indx, a_indx = next(((i,j) for i,x in enumerate(isoch_ma) for j,y\
        in enumerate(x) if y == [m, a]), None)
        # Generate shifted best fit isochrone.
        shift_isoch = move_isoch(sc_params[0],
                                 isoch_list[m_indx][a_indx][:2], e, d)
        # Generate best fit synthetic cluster.
        synth_clst = s_c(err_lst, completeness, sc_params,
                         isoch_list[m_indx][a_indx], [-1., -1., e, d])

            
    else:
        # Pass empty lists to make_plots.
        print 'Skipping synthetic cluster fitting process.'
        isoch_fit_params, isoch_fit_errors, shift_isoch, synth_clst = \
        [[-1.,-1.,-1.,-1.]], [-1.,-1.,-1.,-1.], [], []

    bf_return = [isoch_fit_params, isoch_fit_errors, shift_isoch, synth_clst]

    return bf_return
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 16:10:05 2014

@author: gabriel
"""

import numpy as np
import random

from genetic_algorithm import gen_algor as g_a
from get_IMF_PDF import IMF_PDF as i_p
from move_isochrone import move_isoch


def boostrap_resample(in_list):
    '''
    Resamples the observed cluster to use in the bootstrap process.
    '''
    obs_clust = [random.choice(in_list) for _ in in_list]
    return obs_clust
    


def bf(err_lst, memb_prob_avrg_sort, completeness, ip_list, sc_params):
    '''
    Main function.
    
    Perform a best fitting process to find the cluster's parameters:
    E(B-V), distance (distance modulus), age and metallicity.
    '''
    
    # Obtain the selected IMF's PDF. We run it once because the array only
    # depends on the IMF selected.
    imf_pdf = i_p(sc_params[1])
    # Replace the name of the IMF by its PDF.
    sc_params[1] = imf_pdf
    
    
    # Genetic algorithm parameters.
    # n_pop, n_gen, fdif, p_cross, cr_sel, p_mut, n_el, n_ei, n_es
    params_ga = [4, 20, 1., 0.85, '2P', 0.01, 1, 25, 6]
    
    # Number of times to run the bootstrap block.
    N_B = 2
   
    # List that holds the parameters values obtained by the bootstrap
    # process.
    params_boot = []
    # Begin bootstrap block.
    for i in range(N_B):
        
        # The first pass is done with no resampling to calculate the final
        # values. After that we resample to get the uncertainty in each
        # parameter.
        if i == 0:
            obs_clust = memb_prob_avrg_sort
        else:
            obs_clust = boostrap_resample(memb_prob_avrg_sort)

        # Call algorithm to calculate the likelihoods for the set of
        # isochrones and return the best fitting parameters.
        
        if i==0:
            # Brute force algorithm.
#            isoch_fit_params = brute_force(sys_select, isoch_params, iso_path,
#                                           line_start, indexes, obs_clust,\
#                                           mass_dist)
            # Genetic algorithm.
            ga_return = g_a(err_lst, obs_clust, completeness, ip_list, sc_params,
                                   params_ga)
        else:
            # Brute force.
#            params_boot.append(brute_force(sys_select, isoch_params, iso_path,
#                                           line_start, indexes, obs_clust,\
#                                           mass_dist))
            # Genetic algorithm algorithm.
            params_boot.append(g_a(err_lst, obs_clust, completeness, ip_list,
                                   sc_params, params_ga)[0])

        
    # Calculate errors for each parameter.
    isoch_fit_errors = np.std(params_boot, 0)
    
    # For plotting purposes: generate shifted isochrone.
    isoch_list, isoch_ma = ip_list[0], ip_list[1]
    m, a, e, d = ga_return[0]
    m_indx, a_indx = next(((i,j) for i,x in enumerate(isoch_ma) for j,y in enumerate(x) if y == [m, a]), None)
    shift_isoch = move_isoch(sc_params[0], [isoch_list[m_indx][a_indx][0], isoch_list[m_indx][a_indx][1]], e, d)

    return shift_isoch, ga_return, isoch_fit_errors
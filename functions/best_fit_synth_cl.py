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


def boostrap_resample(stars_list):
    '''
    Resamples the observed cluster with replacement. Used by the bootstrap
    process.
    '''
    obs_clust = [random.choice(stars_list) for _ in stars_list]
    return obs_clust
    


def best_fit(err_lst, memb_prob_avrg_sort, completeness, ip_list, N_b, sc_params,
            ga_params):
    '''
    Perform a best fitting process to find the cluster's parameters:
    E(B-V), distance (distance modulus), age and metallicity.
    '''
    
    # Obtain the selected IMF's PDF. We run it once because the array only
    # depends on the IMF selected.
    imf_pdf = i_p(sc_params[1])
    # Replace the name of the IMF by its PDF.
    sc_params[1] = imf_pdf

    # List that holds the parameters values obtained by the bootstrap
    # process.
    params_boot = []
    # Begin bootstrap block.
    for i in range(max(N_b, 2)):

        # Call algorithm to calculate the likelihoods for the set of
        # isochrones and return the best fitting parameters.
        # The first pass is done with no resampling to calculate the final
        # values. After that we resample to get the uncertainty in each
        # parameter.        
        if i==0:
            # Brute force algorithm.
#            isoch_fit_params = brute_force(sys_select, isoch_params, iso_path,
#                                           line_start, indexes, obs_clust,\
#                                           mass_dist)
            # Genetic algorithm.
            ga_return = g_a(err_lst, memb_prob_avrg_sort, completeness, ip_list,
                            sc_params, ga_params)
        else:
            # Brute force.
#            params_boot.append(brute_force(sys_select, isoch_params, iso_path,
#                                           line_start, indexes, obs_clust,\
#                                           mass_dist))
            # Resample cluster with replacement.
            obs_clust = boostrap_resample(memb_prob_avrg_sort)
            # Genetic algorithm.
            params_boot.append(g_a(err_lst, obs_clust, completeness, ip_list,
                                   sc_params, ga_params)[0])

        
    # Calculate errors for each parameter.
    isoch_fit_errors = np.std(params_boot, 0)
    
    # For plotting purposes: generate shifted best fit isochrone.
    isoch_list, isoch_ma = ip_list[0], ip_list[1]
    m, a, e, d = ga_return[0]
    m_indx, a_indx = next(((i,j) for i,x in enumerate(isoch_ma) for j,y in \
    enumerate(x) if y == [m, a]), None)
    shift_isoch = move_isoch(sc_params[0], [isoch_list[m_indx][a_indx][0],
                                       isoch_list[m_indx][a_indx][1]], e, d)

    return shift_isoch, ga_return, isoch_fit_errors
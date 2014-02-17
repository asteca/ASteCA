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
    


def best_fit(err_lst, memb_prob_avrg_sort, completeness, ip_list, bf_params,
             sc_params, ga_params):
    '''
    Perform a best fitting process to find the cluster's parameters:
    E(B-V), distance (distance modulus), age and metallicity.
    '''
    
    bf_flag, N_b = bf_params
    
    # Check if algorithm should run.
    if bf_flag:
    
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
                isoch_fit_params = g_a(err_lst, memb_prob_avrg_sort,
                                       completeness, ip_list, sc_params,
                                       ga_params)
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
        m, a, e, d = isoch_fit_params[0]
        m_indx, a_indx = next(((i,j) for i,x in enumerate(isoch_ma) for j,y in \
        enumerate(x) if y == [m, a]), None)
        shift_isoch = move_isoch(sc_params[0], [isoch_list[m_indx][a_indx][0],
                                           isoch_list[m_indx][a_indx][1]], e, d)
    else:
        # Pass empty lists to make_plots.
        print 'Skipping synthetic cluster fitting process.'
        shift_isoch, isoch_fit_params, isoch_fit_errors = [], \
        [-1.,-1.,-1.,-1.], [-1.,-1.,-1.,-1.]

    return shift_isoch, isoch_fit_params, isoch_fit_errors
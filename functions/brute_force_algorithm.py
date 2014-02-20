# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 16:38:10 2014

@author: gabriel
"""

import numpy as np
from isoch_likelihood import isoch_likelihood as i_l


def brute_force(err_lst, obs_clust, completeness, ip_list, sc_params,
                ga_params):
    '''
    Brute force algorithm that computes the likelihoods for *all* the defined
    isochrones.
    '''
    
    isoch_list, isoch_ma, isoch_ed, ranges_steps = ip_list
    
    # Create lists with all possible 
    e_lst = isoch_ed[0]
    d_lst = isoch_ed[1]
    
    # Initiate list that will hold the likelihood values telling us how well
    # each isochrone (syhtnetic cluster) fits the observed data.
    score = []
    
    tot_sols, i = len(isoch_list)*len(isoch_list[0])*len(e_lst)*len(d_lst), 0
    flag_25, flag_50, flag_75 = False, False, False
    # Iterate through all metallicities.
    for m, m_isochs in enumerate(isoch_list):
        
        # Iterate through all ages.
        for a, isochs in enumerate(m_isochs):
            
            # Iterate through all extinction values.
            for e in e_lst:
                
                # Iterate through all distance modulus.
                for d in d_lst:

                    # Pass metallicity and age values for plotting purposes.
                    params = [isoch_ma[m][a][0], isoch_ma[m][a][1], e, d]
        
                    # Call likelihood function with m,a,e,d values.
                    likel_val = i_l(err_lst, obs_clust, completeness, sc_params, 
                                    isoch_list[m][a], params)
                                    
                    # Store the likelihood for each synthetic cluster.
                    score.append([likel_val, isoch_ma[m][a][0],
                                  isoch_ma[m][a][1], e, d])
                    i += 1

                    if i+1 >= tot_sols/4 and flag_25 == False:
                        print '  25% done'
                        flag_25 = True
                    elif i+1 >= tot_sols/2 and flag_50 == False:
                        print '  50% done'
                        flag_50 = True
                    elif i+1 >= (tot_sols/2 + tot_sols/4) and flag_75 == False:
                        print '  75% done'
                        flag_75 = True
                    elif i+1 == tot_sols:
                        print '  100% done'    
    
    # Find index of function with smallest likelihood value.
    # This index thus points to the isochrone that best fits the observed
    # group of stars.
    best_fit_indx = np.argmin(zip(*score)[0])
    
    isoch_fit_params = [score[best_fit_indx][1:], score[best_fit_indx][0]]
        
    return isoch_fit_params
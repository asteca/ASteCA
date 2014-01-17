# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 16:38:10 2014

@author: gabriel
"""

def brute_force(sys_select, isoch_list, isoch_params, obs_clust, mass_dist):
    '''
    Brute force algorithm that computes the likelihoods for *all* the defined
    isochrones.
    '''
    # Initiate list that will hold the values (scores) which defines how well
    # each isochrone fits the observed data.
    score = []
    # Iterate through all the tracks defined and stored.
    tik = time.time()
    for indx, isoch in enumerate(isoch_params):

        # Get parameters value from this isochrone.
        m, a, e, d = isoch
        
        # Get stored isochrone corresponding to these values of 'm' and 'a'.
        isoch_ma = isoch_list[indx]
        
        # Call function that returns the score for a given track.
        isoch_score = isoch_likelihood(e, d, isoch_ma, sys_select, obs_clust, mass_dist)
        print isoch, isoch_score
        # Store the scores for each function/track into list.
        score.append(isoch_score)    
        
    # Find index of function with smallest value of the likelihoods.
    # This index thus points to the isochrone that best fits the observed
    # group of stars.
    best_func= np.argmin(score)
#    min_score = score[best_func]

    z_met, age_gyr, e_bv, dis_mod = [i for i in isoch_params[best_func]]
    dist_kpc = round(10**(0.2*(dis_mod+5.))/1000., 2)
    
    print z_met, age_gyr, e_bv, dis_mod
    print 'time: ', time.time()-tik
    raw_input()
    
    isoch_fit_params = [z_met, age_gyr, e_bv, dis_mod, dist_kpc]         
        
    return isoch_fit_params
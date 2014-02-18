# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 17:21:58 2014

@author: gabriel
"""

from synth_cluster import synth_clust as s_c
import numpy as np


def likelihood(synth_clust, obs_clust):
    '''
    Takes a synthetic cluster, compares it to the observed cluster and
    returns the weighted (log) likelihood value.
    This function follows the recipe given in Monteiro, Dias & Caetano (2010).
    '''
    
    if not synth_clust.any():
        isoch_score = 10000.
    else:
        # Store observed and synthetic clusters as arrays.
        obs_arr = np.array(obs_clust)
        syn_arr = np.array(zip(*synth_clust))
        clust_stars_probs = []
        
        for star in obs_arr:
            # Get probability for this cluster star.
            
            # Avoid numeric errors if one of the errors is 0.            
            e_col, e_mag = max(star[6], 1e-10), max(star[4], 1e-10)
            
            A = 1./(e_col*e_mag)
            B = -0.5*((star[5]-syn_arr[:,0])/e_col)**2
            C = -0.5*((star[3]-syn_arr[:,1])/e_mag)**2
            synth_stars = A*np.exp(B+C)

            # The final prob for this cluster star is the sum over all synthetic
            # stars. Use 1e-10 to avoid nan and inf values in the calculations
            # that follow.
            sum_synth_stars = max(synth_stars.sum(), 1e-10)
            clust_stars_probs.append(sum_synth_stars)
        
        # Store weights data (membership probabilities) into array.
        weights = np.array([zip(*obs_clust)[7]], dtype=float)   
        # Weight probabilities for each cluster star.
        weighted_probs = clust_stars_probs*weights/len(synth_clust[0])
        
        # Final score: sum log likelihoods for each star in cluster.
        isoch_score = min(-sum(np.log(np.asarray(weighted_probs[0]))), 10000.)
    
    return isoch_score

    
    
def isoch_likelihood(err_lst, obs_clust, completeness, sc_params, isochrone,
                     params):
    '''
    Main function.
    
    Call with an isochrone of given values for metallicity and age and supply
    the extinction and distance modulus values to move that isochrone. Use
    that isochrone to generate a synthetic cluster with those parameters and
    finally compare it wiht the observed cluster.
    
    e, d = extinction, distance modulus.
    '''
        
    # Generate synthetic cluster using this "moved" isochrone and a mass
    # distribution.
    synth_clust = s_c(err_lst, completeness, sc_params, isochrone, params)
    
    # Call function to obtain the likelihood by comparing the synthetic cluster
    # with the observed cluster.
    isoch_lik = likelihood(synth_clust, obs_clust)
    
    return isoch_lik
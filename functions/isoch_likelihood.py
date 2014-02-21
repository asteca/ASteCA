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
    
        # Small value used to replace zeros.
        epsilon = 1e-10
    
        # Split obs_arr.
        P = np.split(obs_arr, 8, axis=1)
        # Square errors in color and magnitude.
        P[6] = np.square(np.maximum(P[6], epsilon)) # color
        P[4] = np.square(np.maximum(P[4], epsilon)) # magnitude
        P = np.hstack(P)
    
        # Split syn_arr.
        Q = np.split(syn_arr, 5, axis=1)
        # Square synthetic photometric errors.
        Q[1] = np.square(Q[1])
        Q[3] = np.square(Q[3])
    
        for star in P:
            # Squares sum of errors.
            e_col_2 = star[6] + Q[1]
            e_mag_2 = star[4] + Q[3]
            # star[5] & Q[0] = colors
            # star[3] & Q[2] = magnitude
            B = np.square(star[5] - Q[0]) / e_col_2
            C = np.square(star[3] - Q[2]) / e_mag_2
            synth_stars = np.exp(-0.5 * (B + C)) / np.sqrt(e_col_2 * e_mag_2)
            # The final prob for this cluster star is the sum over all synthetic
            # stars. Use 1e-10 to avoid nan and inf values in the calculations
            # that follow.
            clust_stars_probs.append(max(synth_stars.sum(), epsilon))
        
        
        # Store weights data (membership probabilities) into array.
        weights = np.array([zip(*obs_clust)[7]], dtype=float)   
        # Weight probabilities for each cluster star.
        weighted_probs = clust_stars_probs*weights/len(synth_clust[0])
        
        # Final score: sum log likelihoods for each star in cluster.
        isoch_score = min(-sum(np.log(np.asarray(weighted_probs[0]))), 10000.)
    
    return isoch_score

    
    
def isoch_likelihood(err_lst, obs_clust, completeness, sc_params, isochrone,
                     params, sys_sel):
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
    synth_clust = s_c(err_lst, completeness, sc_params, isochrone, params,
                      sys_sel)
    
    # Call function to obtain the likelihood by comparing the synthetic cluster
    # with the observed cluster.
    isoch_lik = likelihood(synth_clust, obs_clust)
    
    return isoch_lik
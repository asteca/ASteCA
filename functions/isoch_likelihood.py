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
            
            # e_col, e_mag = star[6], star[4]
            A = 1./(star[6]*star[4])
            B = -0.5*((star[5]-syn_arr[:,0])/star[6])**2
            C = -0.5*((star[3]-syn_arr[:,1])/star[4])**2
            synth_stars = A*np.exp(B+C)
        
            # The final prob for this cluster star is the sum over all synthetic
            # stars. Use 1e-10 to avoid nan and inf values in the calculations that
            # follow.
            sum_synth_stars = max(synth_stars.sum(), 1e-10) 
            clust_stars_probs.append(sum_synth_stars)
        
        # Store weights data (membership probabilities) into array.
        weights = np.array([zip(*obs_clust)[7]], dtype=float)   
        # Weight probabilities for each cluster star.
        weighted_probs = clust_stars_probs*weights/len(synth_clust[0])
        
        # Final score: sum log likelihoods for each star in cluster.
        isoch_score = -sum(np.log(np.asarray(weighted_probs[0])))
    
    return isoch_score

    
    
def isoch_likelihood(sys_select, isochrone, e, d, obs_clust, mass_dist,
                     completeness, f_bin, q_bin, popt_mag, popt_col1):
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
    synth_clust = s_c(sys_select, isochrone, e, d, mass_dist, completeness,
                      f_bin, q_bin, popt_mag, popt_col1)
    
    # Call function to obtain the likelihood by comparing the synthetic cluster
    # with the observed cluster.
    isoch_lik = likelihood(synth_clust, obs_clust)
    
    return isoch_lik
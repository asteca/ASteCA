# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 17:21:58 2014

@author: gabriel
"""

from move_isochrone import move_isoch
import numpy as np
import time

def likelihood(synth_clust, obs_clust):
    '''
    Takes a synthetic cluster, compares it to the observed cluster and
    returns the weighted (log) likelihood value.
    This function follows the recipe given in Monteiro, Dias & Caetano (2010).
    '''
    
    # Store cluster data in new arrays: color & magnitude, their errors and 
    # each star's membership probabilities (weights)
    # Store color and magnitude into array.
#    data0 = zip(*obs_clust)
#    col_mag = np.array([data0[5], data0[3]], dtype=float)
#    # Store color and magnitude errors into array.
#    err_col_mag = np.array([data0[6], data0[4]], dtype=float)
#    print 'synth_clust', synth_clust, '\n'

    obs_arr = np.array(obs_clust)
    syn_arr = np.array(zip(*synth_clust))
    
    clust_stars_probs = []
    for star in obs_arr:
        # Get probability for this cluster star.
        
        # sigma_c, sigma_m = star[6], star[4]
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
    weighted_probs = clust_stars_probs*weights
    
    # Final score: sum log likelihoods for each star in cluster.
    isoch_score = -sum(np.log(np.asarray(weighted_probs[0])))
    
    return isoch_score



def synthetic_clust(isochrone, mass_dist):
    '''
    Takes an isochrone and returns a synthetic cluster created according to
    a certain mass distribution.
    '''
    
    # Interpolate masses from the mass distribution into the isochrone to
    # obtain its magnitude and color values. Reject masses that fall outside of
    # the isochrone mass range.
    
    # Randomly select a fraction of these stars to be binaries.
    
    # Calculate the secondary masses of these binary stars.
    
    # Add masses and update array.
    
    # Randomly move stars according to given error distributions.
    
    # Remove stars according to a completness limit function.

    synth_clust = isochrone
    
    return synth_clust
    
    
    
def isoch_likelihood(sys_select, isochrone, e, d, obs_clust, mass_dist):
    '''
    Main function.
    
    Call with an isochrone of given values for metallicity and age and supply
    the extinction and distance modulus values to move that isochrone. Use
    that isochrone to generate a synthetic cluster with those parameters and
    finally compare it wiht the observed cluster.
    
    e, d = extinction, distance modulus.
    '''
    
#    tik = time.time()
    # Store isochrone moved by the values 'e' and 'd'.
    isoch_final = move_isoch(sys_select, isochrone, e, d)
    
    # Generate synthetic cluster using this "moved" isochrone and a mass
    # distribution.
    synth_clust = synthetic_clust(isoch_final, mass_dist)
    
    # Call function to obtain the likelihood by comparing the synthetic cluster
    # with the observed cluster.
    
    isoch_lik = likelihood(synth_clust, obs_clust)
#    print '  likel', time.time()-tik
    
    return isoch_lik
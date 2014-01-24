# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 17:21:58 2014

@author: gabriel
"""

from move_isochrone import move_isoch
import numpy as np
#import time

def likelihood(synth_clust, obs_clust):
    '''
    Takes a synthetic cluster, compares it to the observed cluster and
    returns the weighted (log) likelihood value.
    This function follows the recipe given in Monteiro, Dias & Caetano (2010).
    '''
    
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
    N = 100
    col, mag, mass = np.linspace(0, 1, len(isochrone[0])), np.linspace(0, 1, N),\
    np.linspace(0, 1, N)
    # One-dimensional linear interpolation.
    col_i, mag_i, mass_i = (np.interp(mag, col, isochrone[i]) for i in range(3))
#    mag_i = np.interp(mag, col, isochrone[1])
#    mass_i = np.interp(mag, col, isochrone[2])
    # Store track interpolated values.
    isoch_inter = np.asarray([col_i, mag_i, mass_i])
    
    print isoch_inter[0], '\n'
    print isoch_inter[1], '\n'
    print isoch_inter[2], '\n'
#    raw_input()


#    import bisect
#    # Reject masses out of range.
#    masses = bisect.bisect(mass_i, mass_dist)

    print mass_dist, '\n'
    
    isoch_m_d = [[], [], []]
    min_m, max_m = min(isochrone[2]), max(isochrone[2])
    for m in mass_dist:
        if min_m <= m <= max_m:
            indx, m_i = min(enumerate(isoch_inter[2]), key=lambda x:abs(x[1]-m))
            isoch_m_d[0].append(isoch_inter[0][indx])
            isoch_m_d[1].append(isoch_inter[1][indx])
            isoch_m_d[2].append(m_i)
        
    print isoch_m_d[0], '\n'
    print isoch_m_d[1], '\n'
    print isoch_m_d[2], '\n'
    raw_input()
    
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
    # Store isochrone moved by the values 'e' and 'd'.
    isoch_moved = move_isoch(sys_select, isochrone, e, d)
    isoch_final = isoch_moved + [isochrone[2]]
        
    # Generate synthetic cluster using this "moved" isochrone and a mass
    # distribution.
    synth_clust = synthetic_clust(isoch_final, mass_dist)
    
    # Call function to obtain the likelihood by comparing the synthetic cluster
    # with the observed cluster.
    isoch_lik = likelihood(synth_clust, obs_clust)
    
    return isoch_lik
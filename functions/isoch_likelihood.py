# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 17:21:58 2014

@author: gabriel
"""

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
   
    clust_stars_probs = []
    for indx,star in enumerate(obs_clust):
        # The first term does not depend on the synth stars.
        sigma_c, sigma_m = star[6], star[4]
        A = 1/(sigma_c*sigma_m)
        
        # Get probability for this cluster star.
        synth_stars = []
        for synth_st in zip(*synth_clust):
            # synth_st[0] = color ; synth_st[1] = mag
            B = np.exp(-0.5*((star[5]-synth_st[0])/sigma_c)**2)
            C = np.exp(-0.5*((star[3]-synth_st[1])/sigma_m)**2)
            synth_stars.append(A*B*C)
            
        # The final prob for this cluster star is the sum over all synthetic
        # stars.
        sum_synth_stars = sum(synth_stars) if sum(synth_stars)>0. else 1e-06
        clust_stars_probs.append(sum_synth_stars)
        
    # Store weights data (membership probabilities) into array.
    weights = np.array([zip(*obs_clust)[7]], dtype=float)   
    # Weight probabilities for each cluster star.
    weighted_probs = clust_stars_probs*weights
    
    # Get weighted likelihood.
#    L_x = reduce(lambda x, y: x*y, weighted_probs)
    
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
    
    
    
def isoch_likelihood(m, a, e, d, sys_select, obs_clust, mass_dist):
    '''
    Main function.
    
    Call with given values for metallicity, age, extinction and distance modulus
    to generate a synthetic cluster with those parameters and compare it wiht
    the observed cluster.
    
    m, a, e, d = metallicity, age, extinction, distance modulus.
    '''
    
    # Store isochrone of metallicity value 'm' and age 'a' moved by the
    # values 'e' and 'd'.
    isoch_final = move_track(isoch_list, sys_select, e, d)
#    isoch_final = read_isoch(m, a, e, d, sys_select, iso_path, line_start,
#                             indexes)
    
    # Generate synthetic cluster using this "moved" isochrone and a mass
    # distribution.
    synth_clust = synthetic_clust(isoch_final, mass_dist)
    
    # Call function to obtain the likelihood by comparing the synthetic cluster
    # with the observed cluster.
    isoch_lik = likelihood(synth_clust, obs_clust)
    
    return isoch_lik
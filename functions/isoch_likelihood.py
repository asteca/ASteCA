# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 17:21:58 2014

@author: gabriel
"""

from move_isochrone import move_isoch
import numpy as np
import itertools
#import time

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
        weighted_probs = clust_stars_probs*weights
        
        # Final score: sum log likelihoods for each star in cluster.
        isoch_score = -sum(np.log(np.asarray(weighted_probs[0])))
    
    return isoch_score



def synthetic_clust(isochrone, mass_dist, completeness):
    '''
    Takes an isochrone and returns a synthetic cluster created according to
    a certain mass distribution.
    '''
    
    # Interpolate extra color, magnitude and masses into the isochrone.
    N = 1500
    col, mag, mass = np.linspace(0, 1, len(isochrone[0])), np.linspace(0, 1, N),\
    np.linspace(0, 1, N)
    # One-dimensional linear interpolation.
    col_i, mag_i, mass_i = (np.interp(mag, col, isochrone[i]) for i in range(3))
    # Store isochrone's interpolated values.
    isoch_inter = np.asarray([col_i, mag_i, mass_i])

    # Remove stars from isochrone with magnitude values larger that the maximum
    # value found in the observation (entire field, not just the cluster region).
    #
    # Sort isochrone first according to magnitude values (min to max).
    isoch_sort = zip(*sorted(zip(*isoch_inter), key=lambda x: x[1]))
    # Now remove values beyond max_mag (= completeness[0]).
    # Get index of closest mag value to max_mag.
    max_indx = min(range(len(isoch_sort[1])), key=lambda i: abs(isoch_sort[1][i]-completeness[0]))
    # Remove elements.
    isoch_cut = [isoch_sort[i][0:max_indx] for i in range(3)]
    
    
    # Find the mass in the isochrone closest to each mass in dist_mass while
    # rejecting those masses that fall outside of the isochrone's mass range. 
    isoch_m_d = [[], [], []]
    min_m, max_m = min(isoch_cut[2]), max(isoch_cut[2])
    for m in mass_dist:
        if min_m <= m <= max_m:
            indx, m_i = min(enumerate(isoch_cut[2]), key=lambda x:abs(x[1]-m))
            isoch_m_d[0].append(isoch_cut[0][indx])
            isoch_m_d[1].append(isoch_cut[1][indx])
            isoch_m_d[2].append(m_i)

#    print isoch_m_d[1], '\n'
    

    ### Completeness limit removal of stars. ###
    if not isoch_m_d[0]:
        # If the isochrone is empty after removing stars outside of the observed
        # ranges, then pass empty array.
        synth_clust = np.asarray([])
    else:
        # If stars exist in the isochrone beyond the completeness magnitude
        # level, then apply the removal of stars. Otherwise, skip it.
        if max(isoch_m_d[1]) > completeness[1][completeness[2]]:
    
    #        print 'peak mag', completeness[1][completeness[2]], '\n'
    #        print 'mag', len(isoch_m_d[1]), isoch_m_d[1], '\n'
            
            # Get histogram. completeness[2] = bin_edges of the observed region
            # histogram.
            synth_mag_hist, bin_edges = np.histogram(isoch_m_d[1], completeness[1])
            pi = completeness[3]
            n1, p1 = synth_mag_hist[completeness[2]], pi[0]
            di = np.around((synth_mag_hist[completeness[2]:]-(n1/p1)*np.asarray(pi)), 0)
    #        print 'di:', di, '\n'
                
                
            # Store indexes of *all* elements in isoch_m_d whose magnitude value
            # falls between the ranges given.
    #        print 'rango', completeness[1][completeness[2]:], '\n'
            rang_indx = [[] for _ in range(len(completeness[1][completeness[2]:])-1)]
            for indx,elem in enumerate(isoch_m_d[1]):
                for i in range(len(completeness[1][completeness[2]:])-1):
                    if completeness[1][completeness[2]+i] < elem <= completeness[1][completeness[2]+(i+1)]:
                        rang_indx[i].append(indx)
                    
    #        print 'rang_indx', rang_indx, '\n'
            
            # Pick a number (given by the list 'di') of random elements in each
            # range. Those are the indexes of the elements that should be removed
            # from the three sub-lists.
            rem_indx = []
            for indx,num in enumerate(di):
                if rang_indx[indx] and len(rang_indx[indx]) >= num:
                    rem_indx.append(np.random.choice(rang_indx[indx], num, replace=False))
#                    rem_indx.append(rang_indx[indx][:int(num)])
                else:
                    rem_indx.append(rang_indx[indx])
            
    #        print 'rem_indx', rem_indx, '\n'
            # Remove items from list.
            # itertools.chain() flattens the list of indexes and
            # sorted() with reverse=True inverts them so we don't change the indexes of
            # the elements in the lists after removing them.
            d = sorted(list(itertools.chain(*rem_indx)), reverse=True)
    #        print 'd', d, '\n'
            # Remove those selected indexes from the three sub-lists.
            synth_clust = np.delete(np.asarray(isoch_m_d), d, axis=1)
    #        print 'syn', len(synth_clust[1]), synth_clust[1]
    #        raw_input()
        else:
            synth_clust = np.asarray(isoch_m_d)
            
    #    print synth_clust[0], '\n'
    #    print synth_clust[1], '\n'
    #    print synth_clust[2], '\n'
    #    raw_input()
    
        # Randomly select a fraction of these stars to be binaries.
        
        # Calculate the secondary masses of these binary stars.
        
        # Add masses and update array.
        
        # Randomly move stars according to given error distributions.
    
    #    synth_clust = isochrone
    
    return synth_clust
    
    
    
def isoch_likelihood(sys_select, isochrone, e, d, obs_clust, mass_dist, completeness):
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
    synth_clust = synthetic_clust(isoch_final, mass_dist, completeness)
    
    # Call function to obtain the likelihood by comparing the synthetic cluster
    # with the observed cluster.
    isoch_lik = likelihood(synth_clust, obs_clust)
    
    return isoch_lik
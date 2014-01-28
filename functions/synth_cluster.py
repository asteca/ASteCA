# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 15:22:10 2014

@author: gabriel
"""

import numpy as np
import itertools
import time


def find_closest(A, target):
    '''
    Helping function for mass interpolating into the isochrone.
    Find closest target element for elements in A.
    '''
    #A must be sorted
    idx = A.searchsorted(target)
    idx = np.clip(idx, 1, len(A)-1)
    left = A[idx-1]
    right = A[idx]
    # target - left < right - target is True (or 1) when target is closer to
    # left and False (or 0) when target is closer to right
    idx -= target - left < right - target
    return idx


def synth_clust(isochrone, mass_dist, completeness):
    '''
    Takes an isochrone and returns a synthetic cluster created according to
    a certain mass distribution.
    '''
    
    # Interpolate extra color, magnitude and masses into the isochrone.
    N = 1000
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


#    isoch_m_d0 = [[], [], []]
#    min_m, max_m = min(isoch_cut[2]), max(isoch_cut[2])
#    print 'min, max mass', min(isoch_cut[2]), max(isoch_cut[2]), '\n'
#    for m in mass_dist:
#        if min_m <= m <= max_m:
#            indx, m_i = min(enumerate(isoch_cut[2]), key=lambda x:abs(x[1]-m))
#            isoch_m_d0[0].append(isoch_cut[0][indx])
#            isoch_m_d0[1].append(isoch_cut[1][indx])
#            isoch_m_d0[2].append(m_i)
#    print len(isoch_m_d0[2]), isoch_m_d0[2], '\n'
          
    # Find the mass in the isochrone closest to each mass in dist_mass while
    # rejecting those masses that fall outside of the isochrone's mass range.
    # Convert to arrays.            
    data = np.array(isoch_cut)
    target = np.array(mass_dist)    
    # Returns the indices that would sort the array.
    order = data[2, :].argsort()
    key = data[2, order]
    target = target[(target >= key[0]) & (target <= key[-1])]
    # Call function to return closest elements (indexes)
    closest = find_closest(key, target)
    # Store values in array.
    isoch_m_d = data[:, order[closest]]
            
    
    # Completeness limit removal of stars. Remove a number of stars according
    # to the percentages of star loss find in get_completeness for the
    # real observation.
    # Check for an empty array.
    if not isoch_m_d.any():
        # If the isochrone is empty after removing stars outside of the observed
        # ranges, then pass an empty array.
        synth_clust = np.asarray([])
    else:
        # If stars exist in the isochrone beyond the completeness magnitude
        # level, then apply the removal of stars. Otherwise, skip it.
        if max(isoch_m_d[1]) > completeness[1][completeness[2]]:
           
            # Get histogram. completeness[2] = bin_edges of the observed region
            # histogram.
            synth_mag_hist, bin_edges = np.histogram(isoch_m_d[1], completeness[1])
            pi = completeness[3]
            n1, p1 = synth_mag_hist[completeness[2]], pi[0]
            di = np.around((synth_mag_hist[completeness[2]:]-(n1/p1)*np.asarray(pi)), 0)

            # Store indexes of *all* elements in isoch_m_d whose magnitude value
            # falls between the ranges given.
            rang_indx = [[] for _ in range(len(completeness[1][completeness[2]:])-1)]
            for indx,elem in enumerate(isoch_m_d[1]):
                for i in range(len(completeness[1][completeness[2]:])-1):
                    if completeness[1][completeness[2]+i] < elem <= completeness[1][completeness[2]+(i+1)]:
                        rang_indx[i].append(indx)
#            print 'r_indx', rang_indx, '\n'
            
            # Pick a number (given by the list 'di') of random elements in each
            # range. Those are the indexes of the elements that should be removed
            # from the three sub-lists.
            rem_indx = []
            for indx,num in enumerate(di):
                if rang_indx[indx] and len(rang_indx[indx]) >= num:
                    rem_indx.append(np.random.choice(rang_indx[indx], num, 
                                                     replace=False))
#                    rem_indx.append(rang_indx[indx][:int(num)])
                else:
                    rem_indx.append(rang_indx[indx])
#            print 'rem_indx', rem_indx, '\n'
#            raw_input()
            
            # Remove items from list.
            # itertools.chain() flattens the list of indexes and sorted() with
            # reverse=True inverts them so we don't change the indexes of
            # the elements in the lists after removing them.
            d = sorted(list(itertools.chain(*rem_indx)), reverse=True)
            # Remove those selected indexes from the three sub-lists.
            synth_clust = np.delete(np.asarray(isoch_m_d), d, axis=1)
        else:
            synth_clust = np.asarray(isoch_m_d)
            
#        print synth_clust[0], '\n'
#        print synth_clust[1], '\n'
#        print synth_clust[2], '\n'
#        raw_input()
    
        # Randomly select a fraction of the stars left to be binaries.
        
        # Calculate the secondary masses of these binary stars.
        
        # Add masses and update array.
        
        # Randomly move stars according to given error distributions.
    
    #    synth_clust = isochrone
    
    return synth_clust
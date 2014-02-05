# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 15:22:10 2014

@author: gabriel
"""

from move_isochrone import move_isoch

import numpy as np
import random
import itertools

#import matplotlib.pyplot as plt


def exp_func(x, a, b, c):
    '''
    Exponential function.
    '''
    return a * np.exp(b * x) + c
    
    
def gauss_error(col, e_col, mag, e_mag):
    '''
    Randomly move mag and color through a Gaussian function.
    '''
    col_gauss = col + np.random.normal(0, 1, len(col))*e_col
    mag_gauss = mag + np.random.normal(0, 1, len(col))*e_mag
    
    return col_gauss, mag_gauss
    

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


def mass_interp(isochrone, mass_dist):
    '''
    For each mass in the mass distribution, find the mass in the isochrone 
    closest to it while rejecting those masses that fall outside of the
    isochrone's mass range.
    '''
    # Convert to arrays.            
    data = np.array(isochrone)
    target = np.array(mass_dist)    
    # Returns the indices that would sort the array.
    order = data[2, :].argsort()
    key = data[2, order]
    target = target[(target >= key[0]) & (target <= key[-1])]
    # Call function to return closest elements (indexes)
    closest = find_closest(key, target)
    # Store values in array and return.
    isoch_interp = data[:, order[closest]]
    
    return isoch_interp
    


def synth_clust(sys_select, isochrone, m_p, a_p, e, d, mass_dist, completeness, f_bin,
                q_bin, popt_mag, popt_col1):
    '''
    Main function.
    
    Takes an isochrone and returns a synthetic cluster created according to
    a certain mass distribution.
    '''
    
#    f, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8)) = plt.subplots(2, 4)
#    for ax in [ax2, ax3, ax4, ax5, ax6, ax7]:
#        ax.set_xlabel('$(C-T_1)$', fontsize=15)
#        ax.set_ylabel('$T_1$', fontsize=15)
    
    
    # Interpolate extra color, magnitude and masses into the isochrone.
    N = 1000
    col, mag, mass = np.linspace(0, 1, len(isochrone[0])), np.linspace(0, 1, N),\
    np.linspace(0, 1, N)
    # One-dimensional linear interpolation.
    col_i, mag_i, mass_i = (np.interp(mag, col, isochrone[i]) for i in range(3))
    # Store isochrone's interpolated values.
    isoch_inter = np.asarray([col_i, mag_i, mass_i])
    
#    ax1.set_title('Isochrone (interp)')
#    ax1.invert_yaxis()
#    ax1.set_xlabel('$(C-T_1)_o$', fontsize=15)
#    ax1.set_ylabel('$M_{T_1}$', fontsize=15)
#    text1 = '$z = %0.4f$' '\n' % m_p
#    text2 = '$log(age) = %0.2f$' % a_p
#    text=text1+text2
#    plt.text(0.1, 0.1, text, transform=ax1.transAxes,
#             bbox=dict(facecolor='white', alpha=0.5), fontsize=15)
#    plt.text(0.1, 0.9, 'a)', transform=ax1.transAxes,
#             bbox=dict(facecolor='white', alpha=0.5), fontsize=16)
#    ax1.scatter(isoch_inter[0], isoch_inter[1], s=10, c='aqua')


    # Move synth cluster with the values 'e' and 'd'.
    isoch_moved = move_isoch(sys_select, [isoch_inter[0], isoch_inter[1]], e, d) + [isoch_inter[2]]
    
#    ax2.set_title('Shifted isochrone')
#    ax2.invert_yaxis()
#    text1 = '$E_{(B-V)} = %0.2f$' '\n' % e
#    text2 = '$(m-M)_o = %0.2f$' % d
#    text=text1+text2
#    plt.text(0.1, 0.1, text, transform=ax2.transAxes,
#             bbox=dict(facecolor='white', alpha=0.5), fontsize=15)
#    plt.text(0.1, 0.9, 'b)', transform=ax2.transAxes,
#             bbox=dict(facecolor='white', alpha=0.5), fontsize=16)
#    ax2.scatter(isoch_moved[0], isoch_moved[1], s=15, c='azure')           
         
         
    # Remove stars from isochrone with magnitude values larger that the maximum
    # value found in the observation (entire field, not just the cluster region).
    #
    # Sort isochrone first according to magnitude values (min to max).
    isoch_sort = zip(*sorted(zip(*isoch_moved), key=lambda x: x[1]))
    # Now remove values beyond max_mag (= completeness[0]).
    # Get index of closest mag value to max_mag.
    max_indx = min(range(len(isoch_sort[1])), key=lambda i: abs(isoch_sort[1][i]-completeness[0]))
    # Remove elements.
    isoch_cut = np.array([isoch_sort[i][0:max_indx] for i in range(3)])

#    ax3.set_title('Max magnitude cut')
#    ax3.invert_yaxis()
#    plt.text(0.1, 0.9, 'c)', transform=ax3.transAxes,
#             bbox=dict(facecolor='white', alpha=0.5), fontsize=16)
#    ax3.scatter(isoch_cut[0], isoch_cut[1], s=15, c='bisque')     
         

    # Interpolate masses in mass_dist into the isochrone rejecting those
    # masses that fall outside of the isochrone's mass range.
    isoch_m_d = mass_interp(isoch_cut, mass_dist)

#    ax4.set_title('IMF masses')
#    ax4.invert_yaxis()
#    plt.text(0.1, 0.9, 'd)', transform=ax4.transAxes,
#             bbox=dict(facecolor='white', alpha=0.5), fontsize=16)
#    ax4.scatter(isoch_m_d[0], isoch_m_d[1], s=15, c='teal')     
    

    # Assignment of binarity.
    # Randomly select a fraction of stars to be binaries.
    # Indexes of the randomly selected stars in isoch_m_d.
    bin_indxs = random.sample(range(len(isoch_m_d[0])), int(f_bin*len(isoch_m_d[0])))
    
    # Calculate the secondary masses of these binary stars between q_bin*m1
    # and m1, where m1 is the primary mass.
    # Primary masses.
    m1 = np.asarray(isoch_m_d[2][bin_indxs])
    # Secondary masses.
    mass_bin0 = np.random.uniform(q_bin*m1, m1)
    # If any secondary mass falls outside of the lower isochrone's mass range,
    # change its value to the min value.
    if not type(mass_bin0) is float:
        # This prevents a rare error where apparently mass_bin0 is a float.
        mass_bin = [i if i >= min(isoch_m_d[2]) else min(isoch_m_d[2]) for i in mass_bin0]

        # Find color and magnitude values for each secondary star. This will
        # slightly change the values of the masses since they will be assigned
        # to the closest value found in the interpolated isochrone.
        bin_isoch = mass_interp(isoch_cut, mass_bin)
        
        # Obtain color, magnitude and masses for each binary system.
        # Transform color to the first filter's magnitude before obtaining the
        # new binary magnitude.
        col_mag_bin = -2.5*np.log10(10**(-0.4*(isoch_m_d[0][bin_indxs]+isoch_m_d[1][bin_indxs]))+10**(-0.4*(bin_isoch[0]+bin_isoch[1])))
        mag_bin = -2.5*np.log10(10**(-0.4*isoch_m_d[1][bin_indxs])+10**(-0.4*bin_isoch[1]))
        # Transform back first filter's magnitude into color.
        col_bin = col_mag_bin - mag_bin
        
        # Add masses to obtain the system's mass.
        mass_bin = isoch_m_d[2][bin_indxs] + bin_isoch[2]
        
        # Update array with new values of color, magnitude and masses.
        for indx,i in enumerate(bin_indxs):
            isoch_m_d[0][i] = col_bin[indx]
            isoch_m_d[1][i] = mag_bin[indx]
            isoch_m_d[2][i] = mass_bin[indx]


#    ax5.set_title('Binarity')
#    ax5.invert_yaxis()
#    plt.text(0.1, 0.9, 'e)', transform=ax5.transAxes,
#             bbox=dict(facecolor='white', alpha=0.5), fontsize=16)
#    ax5.scatter(isoch_m_d[0], isoch_m_d[1], s=15, c='steelblue')        
    
    
    # Completeness limit removal of stars. Remove a number of stars according
    # to the percentages of star loss find in get_completeness for the
    # real observation.
    # Check for an empty array.
    if not isoch_m_d.any():
        # If the isochrone is empty after removing stars outside of the observed
        # ranges, then pass an empty array.
        clust_compl = np.asarray([])
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
            
            # Pick a number (given by the list 'di') of random elements in each
            # range. Those are the indexes of the elements that should be removed
            # from the three sub-lists.
            rem_indx = []
            for indx,num in enumerate(di):
                if rang_indx[indx] and len(rang_indx[indx]) >= num:
                    rem_indx.append(np.random.choice(rang_indx[indx], num, 
                                                     replace=False))
                else:
                    rem_indx.append(rang_indx[indx])
            
            # Remove items from list.
            # itertools.chain() flattens the list of indexes and sorted() with
            # reverse=True inverts them so we don't change the indexes of
            # the elements in the lists after removing them.
            d = sorted(list(itertools.chain(*rem_indx)), reverse=True)
            # Remove those selected indexes from the three sub-lists.
            clust_compl = np.delete(np.asarray(isoch_m_d), d, axis=1)
        else:
            clust_compl = np.asarray(isoch_m_d)

#        ax6.set_title('Completeness')
#        ax6.invert_yaxis()
#        plt.text(0.1, 0.9, 'f)', transform=ax6.transAxes,
#                 bbox=dict(facecolor='white', alpha=0.5), fontsize=16)
#        ax6.scatter(clust_compl[0], clust_compl[1], s=15, c='blue')

        
        # Randomly move stars according to given error distributions.
        # Get errors according to errors distribution.
        sigma_mag = np.array(exp_func(clust_compl[1], *popt_mag))
        sigma_col = np.array(exp_func(clust_compl[1], *popt_col1))
        col_gauss, mag_gauss = gauss_error(clust_compl[0], 2*sigma_col, clust_compl[1], 2*sigma_mag)
        clust_error = [col_gauss, mag_gauss]
        
#        ax7.set_title('Add errors')
#        ax7.invert_yaxis()
#        plt.text(0.1, 0.9, 'g)', transform=ax7.transAxes,
#                 bbox=dict(facecolor='white', alpha=0.5), fontsize=16)
#        ax7.scatter(clust_error[0], clust_error[1], s=15, c='black')

#        plt.show()
#        raw_input()
       
        # Append masses.
        synth_clust = clust_error + [clust_compl[2]]
    
    return np.array(synth_clust)
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 15:22:10 2014

@author: gabriel
"""

from move_isochrone import move_isoch
from get_mass_dist import mass_dist as m_d
#from synth_plot import synth_clust_plot as s_c_p

import numpy as np
import random
import itertools


def exp_func(x, a, b, c):
    '''
    Exponential function.
    '''
    return a * np.exp(b * x) + c


def gauss_error(col, e_col, mag, e_mag):
    '''
    Randomly move mag and color through a Gaussian function.
    '''
    col_gauss = col + np.random.normal(0, 1, len(col)) * e_col
    mag_gauss = mag + np.random.normal(0, 1, len(col)) * e_mag

    return col_gauss, mag_gauss


def find_closest(A, target):
    '''
    Helping function for mass interpolating into the isochrone.
    Find closest target element for elements in A.
    '''
    #A must be sorted
    idx = A.searchsorted(target)
    idx = np.clip(idx, 1, len(A) - 1)
    left = A[idx - 1]
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


def synth_clust(err_lst, completeness, sc_params, isochrone, params, sys_sel):
    '''
    Main function.

    Takes an isochrone and returns a synthetic cluster created according to
    a certain mass distribution.
    '''

    mass_params, f_bin, q_bin = sc_params[0:2], sc_params[2], sc_params[3]

    # Store mass distribution used to produce a synthetic cluster based on
    # a given theoretic isochrone.
    mass_dist = m_d(mass_params)

    # Interpolate extra color, magnitude and masses into the isochrone.
    N = 1000
    col, mag, mass = np.linspace(0, 1, len(isochrone[0])), \
    np.linspace(0, 1, N), np.linspace(0, 1, N)
    # One-dimensional linear interpolation.
    col_i, mag_i, mass_i = (np.interp(mag, col, isochrone[i]) for i in range(3))
    # Store isochrone's interpolated values.
    isoch_inter = np.asarray([col_i, mag_i, mass_i])

    # Move synth cluster with the values 'e' and 'd'.
    e, d = params[2], params[3]
    isoch_moved = move_isoch(sys_sel, [isoch_inter[0], isoch_inter[1]], e, d) +\
    [isoch_inter[2]]

    # Remove stars from isochrone with magnitude values larger that the maximum
    # value found in the observation (entire field, not just the cluster
    # region).
    #
    # Sort isochrone first according to magnitude values (min to max).
    isoch_sort = zip(*sorted(zip(*isoch_moved), key=lambda x: x[1]))
    # Now remove values beyond max_mag (= completeness[0]).
    # Get index of closest mag value to max_mag.
    max_indx = min(range(len(isoch_sort[1])), key=lambda i:
    abs(isoch_sort[1][i] - completeness[0]))
    # Remove elements.
    isoch_cut = np.array([isoch_sort[i][0:max_indx] for i in range(3)])

    # Check for an empty array.
    if isoch_cut.any():
        # Interpolate masses in mass_dist into the isochrone rejecting those
        # masses that fall outside of the isochrone's mass range.
        isoch_m_d = mass_interp(isoch_cut, mass_dist)
    else:
        # If the isochrone is empty after the magnitude cut, pass empty array.
        isoch_m_d = np.asarray([])

    # For plotting purposes: store a copy of this list before adding binaries.
#    from copy import deepcopy
#    isoch_m_d0 = deepcopy(isoch_m_d)

    # Check for an empty array.
    if isoch_m_d.any():

        # Assignment of binarity.
        # Randomly select a fraction of stars to be binaries.
        # Indexes of the randomly selected stars in isoch_m_d.
        bin_indxs = random.sample(range(len(isoch_m_d[0])),
                                  int(f_bin * len(isoch_m_d[0])))

        # Calculate the secondary masses of these binary stars between q_bin*m1
        # and m1, where m1 is the primary mass.
        # Primary masses.
        m1 = np.asarray(isoch_m_d[2][bin_indxs])
        # Secondary masses.
        mass_bin0 = np.random.uniform(q_bin * m1, m1)
        # This prevents a rare error where apparently mass_bin0 is a float.
        if not type(mass_bin0) is float:

            # If any secondary mass falls outside of the lower isochrone's mass
            # range, change its value to the min value.
            mass_bin = np.maximum(min(isoch_m_d[2]), mass_bin0)

            # Find color and magnitude values for each secondary star. This will
            # slightly change the values of the masses since they will be
            # assigned to the closest value found in the interpolated isochrone.
            bin_isoch = mass_interp(isoch_cut, mass_bin)

            # Obtain color, magnitude and masses for each binary system.
            # Transform color to the first filter's magnitude before obtaining
            # the new binary magnitude.
            col_mag_bin = -2.5 * np.log10(10 ** (-0.4 *
            (isoch_m_d[0][bin_indxs] + isoch_m_d[1][bin_indxs])) +
            10 ** (-0.4 * (bin_isoch[0] + bin_isoch[1])))
            mag_bin = -2.5 * np.log10(10 ** (-0.4 * isoch_m_d[1][bin_indxs]) +
            10 ** (-0.4 * bin_isoch[1]))
            # Transform back first filter's magnitude into color.
            col_bin = col_mag_bin - mag_bin

            # Add masses to obtain the system's mass.
            mass_bin = isoch_m_d[2][bin_indxs] + bin_isoch[2]

            # Update array with new values of color, magnitude and masses.
            for indx, i in enumerate(bin_indxs):
                isoch_m_d[0][i] = col_bin[indx]
                isoch_m_d[1][i] = mag_bin[indx]
                isoch_m_d[2][i] = mass_bin[indx]

        # Check for an empty array.
        if isoch_m_d.any():

            # Completeness limit removal of stars. Remove a number of stars
            # according to the percentages of star loss find in
            # get_completeness for the real observation.

            # If stars exist in the isochrone beyond the completeness magnitude
            # level, then apply the removal of stars. Otherwise, skip it.
            if max(isoch_m_d[1]) > completeness[1][completeness[2]]:

                # Get histogram. completeness[2] = bin_edges of the observed
                # region histogram.
                synth_mag_hist, bin_edges = np.histogram(isoch_m_d[1],
                                                         completeness[1])
                pi = completeness[3]
                n1, p1 = synth_mag_hist[completeness[2]], pi[0]
                di = np.around((synth_mag_hist[completeness[2]:] -
                (n1 / p1) * np.asarray(pi)), 0)

                # Store indexes of *all* elements in isoch_m_d whose magnitude
                # value falls between the ranges given.
                c_indx = np.searchsorted(completeness[1][completeness[2]:],
                                         isoch_m_d[1], side='left')
                N = len(completeness[1][completeness[2]:])
                mask = (c_indx > 0) & (c_indx < N)
                elements = c_indx[mask]
                indices = np.arange(c_indx.size)[mask]
                sorting_idx = np.argsort(elements, kind='mergesort')
                ind_sorted = indices[sorting_idx]
                x = np.searchsorted(elements, range(N), side='right',
                    sorter=sorting_idx)
                # Indexes.
                rang_indx = [ind_sorted[x[i]:x[i + 1]] for i in range(N - 1)]

                # Pick a number (given by the list 'di') of random elements in
                # each range. Those are the indexes of the elements that
                # should be removed from the three sub-lists.
                rem_indx = []
                for indx, num in enumerate(di):
                    if rang_indx[indx].any() and len(rang_indx[indx]) >= num:
                        rem_indx.append(np.random.choice(rang_indx[indx], num,
                                                         replace=False))
                    else:
                        rem_indx.append(rang_indx[indx])

                # Remove items from list.
                # itertools.chain() flattens the list of indexes and sorted()
                # with reverse=True inverts them so we don't change the
                # indexes of the elements in the lists after removing them.
                d_i = sorted(list(itertools.chain(*rem_indx)), reverse=True)
                # Remove those selected indexes from the three sub-lists.
                clust_compl = np.delete(np.asarray(isoch_m_d), d_i, axis=1)
            else:
                clust_compl = np.asarray(isoch_m_d)

            # Randomly move stars according to given error distributions.
            # Get errors according to errors distribution.
            popt_mag, popt_col1, e_max = err_lst
            sigma_mag = np.array(exp_func(clust_compl[1], *popt_mag))
            sigma_col = np.array(exp_func(clust_compl[1], *popt_col1))
            # Replace all error values greater than e_max with e_max.
            sigma_mag[sigma_mag > e_max] = e_max
            sigma_col[sigma_col > e_max] = e_max
            # Do the same but for a minimum error value. This ensures a
            # reasonable error spread even for low magnitudes.
            sigma_mag[sigma_mag < 0.05] = 0.05
            sigma_col[sigma_col < 0.05] = 0.05
            # Call function to shift stars around these errors.
            col_gauss, mag_gauss = gauss_error(clust_compl[0], sigma_col,
                                               clust_compl[1], sigma_mag)
            clust_error = [col_gauss, sigma_col, mag_gauss, sigma_mag]

            # Append masses.
            synth_clust = np.array(clust_error + [clust_compl[2]])

        else:
            # If the isochrone is empty after removing stars outside of the
            # observed ranges, then pass an empty array.
            synth_clust = np.asarray([])
    else:
        # If the isochrone is empty after the masses interpolation, pass an
        # empty array.
        synth_clust = np.asarray([])

         #Plot diagrams.
#        s_c_p(mass_dist, isoch_inter, params, isoch_moved, isoch_cut,
#              isoch_m_d0, isoch_m_d, clust_compl, clust_error)

    return synth_clust
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 17:21:58 2014

@author: gabriel
"""

from synth_cluster import synth_clust as s_c
import numpy as np


def tolstoy(Q, obs_clust):
    '''
    Takes a synthetic cluster, compares it to the observed cluster and
    returns the weighted (log) likelihood value.
    This function follows the recipe given in Tolstoy & Saha (1996),
    Hernandez & Valls-Gabaud (2008) andMonteiro, Dias & Caetano (2010).
    '''

    if not Q.any():
        # If synthetic cluster is empty, assign high likelihood value.
        likelihood = 10000.
    else:

        # Unpack observed cluster with squared errors and membership
        # probabilities separated into list.
        P, mem_probs = obs_clust[1:]

        # Store synthetic clusters as array.
        #syn_arr = np.array(zip(*(zip(*obs_arr)[:-2])))  # Observed cluster.
        syn_arr = np.array(zip(*Q))
        cl_stars_probs = []

        # Split syn_arr.
        Q = np.split(syn_arr, 5, axis=1)
        # Square synthetic photometric errors.
        Q[1] = np.square(Q[1])
        Q[3] = np.square(Q[3])

        # Small value used to replace zeros.
        epsilon = 1e-10
        for star in P:
            # Squares sum of errors.
            e_col_2 = np.maximum(star[5] + Q[1], epsilon)
            e_mag_2 = np.maximum(star[3] + Q[3], epsilon)
            # star[4] & Q[0] = colors
            # star[2] & Q[2] = magnitudes
            B = np.square(star[4] - Q[0]) / e_col_2
            C = np.square(star[2] - Q[2]) / e_mag_2
            star_prob = np.exp(-0.5 * (B + C)) / np.sqrt(e_col_2 * e_mag_2)
            # The final prob for this cluster star is the sum over all synthetic
            # stars. Use 1e-10 to avoid nan and inf values in the calculations
            # that follow.
            cl_stars_probs.append(max(star_prob.sum(), epsilon))

        # Weight probabilities for each cluster star.
        clust_prob = cl_stars_probs * mem_probs / len(syn_arr)

        # Final score: sum log likelihoods for each star in cluster.
        likelihood = -sum(np.log(np.asarray(clust_prob[0])))

        #n, p = len(P), len(syn_arr)

        # BIC
        #likelihood = 2 * likelihood + p * np.log(n)

        # AIC_c
        #if (n - p) != 1:
            #likelihood = 2 * likelihood + 2 * p + \
            #(2 * p) * (p + 1) / (n - p - 1)
        #else:
            #likelihood = 2 * likelihood + 2 * p

        #print n, p, likelihood
        #import matplotlib.pyplot as plt
        #fig = plt.figure()
        #ax1 = fig.add_subplot(1, 2, 1)
        #ax2 = fig.add_subplot(1, 2, 2)
        #ax1.scatter(zip(*P)[4], zip(*P)[2], c='r')
        #ax2.scatter(Q[0], Q[2], c='b')
        #text = 'N = {}'.format(len(zip(*P)[4]))
        #ax1.text(0.6, 0.9, text, transform=ax1.transAxes)
        #text1 = 'L = {:.2f}\n'.format(likelihood)
        #text2 = 'N = {}'.format(len(syn_arr))
        #text = text1 + text2
        #ax2.text(0.5, 0.9, text, transform=ax2.transAxes)
        #ax1.invert_yaxis()
        #ax2.invert_yaxis()
        #fig.subplots_adjust(hspace=1)
        #plt.show()

    return likelihood


def dolphin(Q, P):
    '''
    Takes a synthetic cluster, compares it to the observed cluster and
    returns the weighted (log) likelihood value.
    This function follows the recipe given in Dolphin (2002).
    '''

    if not Q.any():
        # If synthetic cluster is empty, assign high likelihood value.
        poiss_lkl = 10000.
    else:

        # Observed cluster's histogram.
        cl_histo = P[1]
        # Bin edges for each dimension.
        b_rx, b_ry = P[2]

        # Magnitude and color for the synthetic cluster.
        syn_mags_cols = np.array(zip(*[Q[0], Q[2]]))
        # Histogram of the synthetic cluster, using the bin edges calculated
        # with the observed cluster.
        syn_histo = np.histogramdd(syn_mags_cols, bins=[b_rx, b_ry])[0]

        # Small value used to replace zeros.
        epsilon = 1e-10
        # Obtain inverse logarithmic 'Poisson likelihood ratio'.
        poiss_lkl = len(Q[0])
        for el1 in zip(*(cl_histo, syn_histo)):
            for el2 in zip(*(el1[0], el1[1])):
                c = -1. * el2[0] * np.log(max(el2[1], epsilon))
                poiss_lkl += c

        #print len(Q[0]), poiss_lkl
        #import matplotlib.pyplot as plt
        #fig = plt.figure()
        #ax1 = fig.add_subplot(1, 3, 1)
        #ax2 = fig.add_subplot(1, 3, 2)
        #ax3 = fig.add_subplot(1, 3, 3)
        #ax1.imshow(d_1.transpose(), origin='lower', aspect='auto')
        #ax2.imshow(d_2.transpose(), origin='lower', aspect='auto')
        #ax3.scatter(P[4], P[2], c='r')
        #ax3.scatter(Q[0], Q[2], c='b')
        #text1 = 'chi = {:.2f}\n'.format(poiss_lkl)
        #text2 = 'N = {}'.format(len(Q[0]))
        #text = text1 + text2
        #ax3.text(0.05, 0.9, text, transform=ax3.transAxes)
        #fig.subplots_adjust(hspace=1)
        #plt.show()

    return poiss_lkl


def isoch_likelihood(err_lst, obs_clust, completeness, st_d_bin_mr, isochrone,
                     params, sys_sel):
    '''
    Call with an isochrone of given values for metallicity and age and supply
    the extinction and distance modulus values to move that isochrone. Use
    that isochrone to generate a synthetic cluster with those parameters and
    finally compare it wiht the observed cluster.
    '''

    # Generate synthetic cluster using this "moved" isochrone and a mass
    # distribution.
    synth_clust = s_c(err_lst, completeness, st_d_bin_mr, isochrone, params,
                      sys_sel)

    # Call function to obtain the likelihood by comparing the synthetic cluster
    # with the observed cluster.
    if obs_clust[0] == 'tolstoy':
        likelihood = tolstoy(synth_clust, obs_clust)
    else:
        likelihood = dolphin(synth_clust, obs_clust)

    return likelihood
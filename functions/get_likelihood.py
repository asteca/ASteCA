# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 17:21:58 2014

@author: gabriel
"""

from synth_cluster import synth_clust as s_c
import numpy as np


def lk_func(synth_clust, obs_clust):
    '''
    Takes a synthetic cluster, compares it to the observed cluster and
    returns the weighted (log) likelihood value.
    This function follows the recipe given in Monteiro, Dias & Caetano (2010).
    '''

    if not synth_clust.any():
        likelihood = 10000.
    else:

        # Unpack observed cluster with squared errors and membership
        # probabilities separated into list.
        P, mem_probs = obs_clust

        # Store synthetic clusters as array.
        #syn_arr = np.array(zip(*(zip(*obs_arr)[:-2])))  # Observed cluster.
        syn_arr = np.array(zip(*synth_clust))
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

    return likelihood


def mighell(Q, P):
    '''
    '''

    if not Q.any():
        chi = 10000.
    else:

        # Unpack observed cluster with squared errors and membership
        # probabilities separated into list.
        #P, mem_probs = obs_clust

        # Store synthetic clusters as array.
        #syn_arr = np.array(zip(*(zip(*obs_arr)[:-2])))  # Observed cluster.
        #Q = np.array(zip(*synth_clust))

        # Split syn_arr.
        #Q = np.split(syn_arr, 5, axis=1)

        #P = np.split(P, 7, axis=1)
        #print np.shape(P), np.shape(Q)
        #print len(P[0]), len(Q[0])

        d1 = np.array(zip(*[P[4], P[2]]))
        d2 = np.array(zip(*[Q[0], Q[2]]))

        #print np.shape(zip(*[P[4], P[2]])), np.shape(zip(*[Q[0], Q[2]]))
        #print np.shape(d1), np.shape(d2)

        # Number of bins.
        b = np.sqrt(len(P[0])) * 4

        # Range for the histograms.
        x_min, x_max = min(P[4]), max(P[4])
        y_min, y_max = min(P[2]), max(P[2])
        rang = [np.linspace(x_min, x_max, b), np.linspace(y_min, y_max, b)]

        d_1 = np.histogramdd(d1, bins=rang)[0]
        d_2 = np.histogramdd(d2, bins=rang)[0]

        chi = 0.
        for el1 in zip(*(d_1, d_2)):
            for el2 in zip(*(el1[0], el1[1])):
                c = np.square(el2[0] + min(el2[0], 1) - el2[1]) / (el2[0] + 1)
                chi += c

        print len(Q[0]), chi
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax1 = fig.add_subplot(1, 3, 1)
        ax2 = fig.add_subplot(1, 3, 2)
        ax3 = fig.add_subplot(1, 3, 3)
        ax1.imshow(d_1.transpose(), origin='lower', aspect='auto')
        ax2.imshow(d_2.transpose(), origin='lower', aspect='auto')
        ax3.scatter(P[4], P[2], c='r')
        ax3.scatter(Q[0], Q[2], c='b')
        text1 = 'chi = {:.2f}\n'.format(chi)
        text2 = 'N = {}'.format(len(Q[0]))
        text = text1 + text2
        ax3.text(0.05, 0.9, text, transform=ax3.transAxes)
        fig.subplots_adjust(hspace=1)
        plt.show()

    return chi


def dolphin(Q, P):
    '''
    '''

    if not Q.any():
        poiss_lkl = 10000.
    else:

        d1 = np.array(zip(*[P[4], P[2]]))
        d2 = np.array(zip(*[Q[0], Q[2]]))

        # Number of bins.
        b = np.sqrt(len(P[0])) * 0.7
        #b = np.sqrt(len(Q[0])) * 2

        # Range for the histograms.
        x_min, x_max = min(P[4]), max(P[4])
        y_min, y_max = min(P[2]), max(P[2])
        rang = [np.linspace(x_min, x_max, b), np.linspace(y_min, y_max, b)]

        d_1 = np.histogramdd(d1, bins=rang)[0]
        d_2 = np.histogramdd(d2, bins=rang)[0]

        # Small value used to replace zeros.
        epsilon = 1e-10
        poiss_lkl = 0.
        for el1 in zip(*(d_1, d_2)):
            for el2 in zip(*(el1[0], el1[1])):
                c = el2[1] - el2[0] + el2[0] * np.log(max(el2[0], epsilon) /
                    max(el2[1], epsilon))
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
    Main function.

    Call with an isochrone of given values for metallicity and age and supply
    the extinction and distance modulus values to move that isochrone. Use
    that isochrone to generate a synthetic cluster with those parameters and
    finally compare it wiht the observed cluster.

    e, d = extinction, distance modulus.
    '''

    # Generate synthetic cluster using this "moved" isochrone and a mass
    # distribution.
    synth_clust = s_c(err_lst, completeness, st_d_bin_mr, isochrone, params,
                      sys_sel)

    # Call function to obtain the likelihood by comparing the synthetic cluster
    # with the observed cluster.
    #likelihood = lk_func(synth_clust, obs_clust)
    #likelihood = mighell(synth_clust, obs_clust)
    likelihood = dolphin(synth_clust, obs_clust)

    return likelihood
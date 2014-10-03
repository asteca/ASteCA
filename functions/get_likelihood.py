# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 17:21:58 2014

@author: gabriel
"""

from synth_cluster import synth_clust as s_c
import numpy as np


def lk_func(synth_clust, obs_arr):
    '''
    Takes a synthetic cluster, compares it to the observed cluster and
    returns the weighted (log) likelihood value.
    This function follows the recipe given in Monteiro, Dias & Caetano (2010).
    '''

    if not synth_clust.any():
        likelihood = 10000.
    else:
        # Store observed and synthetic clusters as arrays.
        #obs_arr = np.array(obs_clust)
        syn_arr = np.array(zip(*synth_clust))
        cl_stars_probs = []

        # Small value used to replace zeros.
        epsilon = 1e-10

        # Split obs_arr.
        P = np.split(obs_arr, 7, axis=1)
        # Square errors in color and magnitude.
        P[5] = np.square(P[5])  # color
        P[3] = np.square(P[3])  # magnitude
        P = np.hstack(P)

        # Split syn_arr.
        Q = np.split(syn_arr, 5, axis=1)
        # Square synthetic photometric errors.
        Q[1] = np.square(Q[1])
        Q[3] = np.square(Q[3])

        for star in P:
            # Squares sum of errors.
            e_col_2 = np.maximum(star[5] + Q[1], epsilon)
            e_mag_2 = np.maximum(star[3] + Q[3], epsilon)
            # star[4] & Q[0] = colors
            # star[2] & Q[2] = magnitude
            B = np.square(star[4] - Q[0]) / e_col_2
            C = np.square(star[2] - Q[2]) / e_mag_2
            star_prob = np.exp(-0.5 * (B + C)) / np.sqrt(e_col_2 * e_mag_2)
            # The final prob for this cluster star is the sum over all synthetic
            # stars. Use 1e-10 to avoid nan and inf values in the calculations
            # that follow.
            cl_stars_probs.append(max(star_prob.sum(), epsilon))

        # Store weights data (membership probabilities) into array.
        weights = np.array([zip(*obs_arr)[6]], dtype=float)
        # Weight probabilities for each cluster star.
        clust_prob = cl_stars_probs * weights / len(synth_clust[0])

        # Final score: sum log likelihoods for each star in cluster.
        likelihood = min(-sum(np.log(np.asarray(clust_prob[0]))), 10000.)

    return likelihood


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
    likelihood = lk_func(synth_clust, obs_clust)

    return likelihood
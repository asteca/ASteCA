# -*- coding: utf-8 -*-
"""
Created on Tue Nov 25 16:00:00 2014

@author: gabriel
"""
import numpy as np
from astroML.plotting import hist


def prepare(memb_prob_avrg_sort, lkl_method, bin_method):
    '''
    Prepare observed cluster array here to save time when the algorithm to
    find the best synthetic cluster fit is used.
    '''

    if lkl_method == 'tolstoy':

        # Remove IDs and re-pack converted to array.
        obs_cl = np.array(zip(*zip(*memb_prob_avrg_sort)[1:]), dtype=float)

        # Square errors ans separate membership probabilities. Done here so
        # as to not repeat the same calculations each time a new synthetic
        # cluster is checked.
        P = np.split(obs_cl, 7, axis=1)

        # Square errors in color and magnitude. Store membership probabilities.
        P[3], P[5], mem_probs = np.square(P[3]), np.square(P[5]), \
        np.asarray(P[6])

        # Pass observed cluster data.
        obs_clust = [lkl_method, np.hstack(P), mem_probs]

    else:

        # Remove ID's and zip.
        P = np.array(zip(*memb_prob_avrg_sort)[1:])

        # Obtain bin edges for each dimension.
        if bin_method in ['sturges', 'sqrt']:
            if bin_method == 'sturges':
                b_num = 1 + np.log2(len(P[0]))
            else:
                b_num = np.sqrt(len(P[0]))

            b_rx = np.histogram(P[4], bins=b_num)[1]
            b_ry = np.histogram(P[2], bins=b_num)[1]
        else:
            b_rx = hist(P[4], bins=bin_method)[1]
            b_ry = hist(P[2], bins=bin_method)[1]

        # Zip magnitudes and colors into array.
        cl_mags_cols = np.array(zip(*[P[4], P[2]]))

        # Obtain histogram for observed cluster.
        cl_histo = np.histogramdd(cl_mags_cols, bins=[b_rx, b_ry],
            weights=np.asarray(P[6]))[0]

        # Pass observed cluster data.
        obs_clust = [lkl_method, cl_histo, [b_rx, b_ry]]

    return obs_clust

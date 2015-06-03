# -*- coding: utf-8 -*-
"""
Created on Tue Nov 25 16:00:00 2014

@author: gabriel
"""
import numpy as np
from astroML.plotting import hist
from .._in import get_in_params as g


def sum_log_list(max_val=100000):
    '''
    '''

    sum_log_lst = [0]
    for val in np.arange(1, max_val):
        sum_log_lst.append(sum_log_lst[int(val - 1)] + np.log(val))

    return sum_log_lst


def prepare(memb_prob_avrg_sort):
    '''
    Prepare observed cluster array here to save time when the algorithm to
    find the best synthetic cluster fit is used.
    '''

    lkl_method, bin_method = g.bf_params[2:4]

    if lkl_method == 'tolstoy':

        # Remove IDs and re-pack converted to array.
        obs_cl = np.array(zip(*zip(*memb_prob_avrg_sort)[1:]), dtype=float)

        # Square errors and separate membership probabilities. Done here so
        # as to not repeat the same calculations each time a new synthetic
        # cluster is checked.
        P = np.split(obs_cl, 7, axis=1)

        # Square errors in color and magnitude. Store membership probabilities
        # separately.
        P[3], P[5], mem_probs = np.square(P[3]), np.square(P[5]), \
            np.asarray(P[6])

        # Pass observed cluster data.
        obs_clust = [np.hstack(P), mem_probs]

    elif lkl_method == 'saha':

        SLL = sum_log_list()

        # Remove ID's (to make entire array of floats) and zip.
        P = np.array(zip(*memb_prob_avrg_sort)[1:], dtype='float')
        mag_col_cl = [P[4], P[2]]

        # Obtain bin edges for each dimension, defining a grid.
        bin_edges = []
        if bin_method in ['sturges', 'sqrt']:
            if bin_method == 'sturges':
                b_num = 1 + np.log2(len(mag_col_cl[0]))
            else:
                b_num = np.sqrt(len(mag_col_cl[0]))

            for mag_col in mag_col_cl:
                bin_edges.append(np.histogram(mag_col, bins=b_num)[1])
        else:
            for mag_col in mag_col_cl:
                bin_edges.append(hist(mag_col, bins=bin_method)[1])

        # Zip magnitudes and colors into array.
        cl_mags_cols = np.array(zip(*mag_col_cl))

        # Obtain *weighted* histogram for observed cluster.
        cl_histo = np.histogramdd(
            cl_mags_cols, bins=bin_edges, weights=np.asarray(P[6]))[0]

        # Pass observed cluster data.
        # obs_clust = [cl_histo, bin_edges, SLL]

        # Pass this list instead if plotting in get_likelihood.
        obs_clust = [cl_histo, bin_edges, SLL, mag_col_cl, 'first']

        from get_likelihood import saha
        W_O_O = saha(np.asarray([list(P[4]), [], list(P[2])]), obs_clust)
        print 'W(O, O):', W_O_O

        # Pass this list instead if plotting in get_likelihood.
        obs_clust = [cl_histo, bin_edges, SLL, mag_col_cl, W_O_O]

    else:

        # Remove ID's (to make entire array of floats) and zip.
        P = np.array(zip(*memb_prob_avrg_sort)[1:], dtype='float')
        mag_col_cl = [P[4], P[2]]

        # Obtain bin edges for each dimension, defining a grid.
        bin_edges = []
        if bin_method in ['sturges', 'sqrt']:
            if bin_method == 'sturges':
                b_num = 1 + np.log2(len(mag_col_cl[0]))
            else:
                b_num = np.sqrt(len(mag_col_cl[0]))

            for mag_col in mag_col_cl:
                bin_edges.append(np.histogram(mag_col, bins=b_num)[1])
        else:
            for mag_col in mag_col_cl:
                bin_edges.append(hist(mag_col, bins=bin_method)[1])

        # Zip magnitudes and colors into array.
        cl_mags_cols = np.array(zip(*mag_col_cl))

        # Obtain *weighted* histogram for observed cluster.
        cl_histo = np.histogramdd(
            cl_mags_cols, bins=bin_edges, weights=np.asarray(P[6]))[0]

        # Pass observed cluster data.
        obs_clust = [cl_histo, bin_edges]

        # Pass this list instead if plotting in get_likelihood.
        # obs_clust = [cl_histo, bin_edges, mag_col_cl]

    return obs_clust

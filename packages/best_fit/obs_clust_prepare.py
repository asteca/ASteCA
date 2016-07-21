
import numpy as np
from ..inp import input_params as g
from ..decont_algors.local_cell_clean import bin_edges_f


def main(memb_prob_avrg_sort):
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

    else:

        # Remove ID's (to make entire array of floats) and zip.
        P = np.array(zip(*memb_prob_avrg_sort)[1:], dtype='float')
        mag_col_cl = [P[4], P[2]]

        # Obtain bin edges for each dimension, defining a grid.
        bin_edges = bin_edges_f(bin_method, mag_col_cl)

        # Zip magnitudes and colors into array.
        cl_mags_cols = np.array(zip(*mag_col_cl))

        # Obtain *weighted* histogram for observed cluster.
        cl_histo = np.histogramdd(
            cl_mags_cols, bins=bin_edges, weights=np.asarray(P[6]))[0]

        # Pass observed cluster data.
        obs_clust = [cl_histo, bin_edges]

        # # Pass this list instead if plotting in get_likelihood.
        # obs_clust = [cl_histo, bin_edges, mag_col_cl]

    return obs_clust

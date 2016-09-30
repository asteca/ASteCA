
import numpy as np
from ..decont_algors.local_cell_clean import bin_edges_f


def main(memb_prob_avrg_sort, lkl_method, bin_method):
    '''
    Prepare observed cluster array here to save time before the algorithm to
    find the best synthetic cluster fit is used.
    '''
    if lkl_method == 'tolstoy':
        # Extract photometric data, and membership probabilities.
        # Square errors here to not repeat the same calculations each time a
        # new synthetic cluster is checked.
        mags = zip(*zip(*memb_prob_avrg_sort)[3])
        e_mags = np.square(zip(*zip(*memb_prob_avrg_sort)[4]))
        cols = zip(*zip(*memb_prob_avrg_sort)[5])
        e_cols = np.square(zip(*zip(*memb_prob_avrg_sort)[6]))
        mem_probs = zip(*memb_prob_avrg_sort)[7]

        # Pass observed cluster data.
        obs_clust = [mags, e_mags, cols, e_cols, mem_probs]

    else:
        # Remove ID's (to make entire array of floats) and zip.
        # Use first magnitude and color.
        # TODO
        mag = zip(*zip(*memb_prob_avrg_sort)[1:][2])[0]
        col = zip(*zip(*memb_prob_avrg_sort)[1:][4])[0]
        prob = zip(*memb_prob_avrg_sort)[1:][6]
        mag_col_cl = [col, mag]

        # Obtain bin edges for each dimension, defining a grid.
        bin_edges = bin_edges_f(bin_method, mag_col_cl)

        # Zip magnitudes and colors into array.
        cl_mags_cols = np.array(zip(*mag_col_cl))

        # Obtain *weighted* histogram for observed cluster.
        cl_histo = np.histogramdd(
            cl_mags_cols, bins=bin_edges, weights=np.asarray(prob))[0]

        # Pass observed cluster data.
        obs_clust = [cl_histo, bin_edges]

        # # Pass this list instead if plotting in get_likelihood.
        # obs_clust = [cl_histo, bin_edges, mag_col_cl]

    return obs_clust

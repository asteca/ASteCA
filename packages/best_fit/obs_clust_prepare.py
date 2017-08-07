
import numpy as np
from ..decont_algors.local_cell_clean import bin_edges_f


def main(cl_max_mag, lkl_method, bin_method):
    '''
    Prepare observed cluster array here to save time before the algorithm to
    find the best synthetic cluster fit is used.
    '''

    # Extract photometric data, and membership probabilities. Remove ID's to
    # make entire array of floats.
    mags_cols_cl = [[], []]
    for mag in zip(*zip(*cl_max_mag)[1:][2]):
        mags_cols_cl[0].append(mag)
    for col in zip(*zip(*cl_max_mag)[1:][4]):
        mags_cols_cl[1].append(col)

    # Square errors here to not repeat the same calculations each time a
    # new synthetic cluster is matched.
    e_mags_cols = []
    for e_m in zip(*zip(*cl_max_mag)[1:][3]):
        e_mags_cols.append(np.square(e_m))
    for e_c in zip(*zip(*cl_max_mag)[1:][5]):
        e_mags_cols.append(np.square(e_c))

    # Store membership probabilities here.
    memb_probs = np.array(zip(*cl_max_mag)[1:][6])

    if lkl_method == 'tolstoy':
        # Store and pass to use in likelihood function. The 'obs_st' list is
        # made up of:
        # obs_st = [star_1, star_2, ...]
        # star_i = [phot_1, phot_2, phot_3, ...]
        # phot_j = [phot_val, error]
        # Where 'phot_j' is a photometric dimension (magnitude or color), and
        # 'phot_val', 'error' the associated value and error for 'star_i'.
        obs_st = []
        mags_cols = mags_cols_cl[0] + mags_cols_cl[1]
        for st_phot, st_e_phot in zip(zip(*mags_cols), zip(*e_mags_cols)):
            obs_st.append(zip(*[st_phot, st_e_phot]))
        obs_clust = [obs_st, memb_probs]

    else:
        # Obtain bin edges for each dimension, defining a grid.
        bin_edges = bin_edges_f(bin_method, mags_cols_cl)

        # Put all magnitudes and colors into a single list.
        obs_mags_cols = mags_cols_cl[0] + mags_cols_cl[1]
        # Obtain histogram for observed cluster.
        cl_histo = np.histogramdd(obs_mags_cols, bins=bin_edges)[0]

        # Flatten N-dimensional histogram.
        cl_histo_f = np.array(cl_histo).ravel()

        # Index of bins where n_i = 0 (no observed stars). Used by the
        # 'Dolphin' and 'Mighell' likelihoods.
        cl_z_idx = [cl_histo_f != 0]

        # Remove all bins where n_i = 0 (no observed stars). Used by the
        # 'Dolphin' likelihood.
        cl_histo_f_z = cl_histo_f[cl_z_idx]

        # sum(n_i * ln(n_i)) - N
        dolphin_cst = np.sum(cl_histo_f_z * np.log(cl_histo_f_z)) -\
            len(obs_mags_cols[0])

        obs_clust = [bin_edges, cl_histo, cl_histo_f, cl_z_idx, cl_histo_f_z,
                     dolphin_cst]

    return obs_clust

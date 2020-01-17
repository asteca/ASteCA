
import numpy as np
from scipy.stats import binned_statistic_dd, gaussian_kde
from ..decont_algors.local_cell_clean import bin_edges_f


def dataProcess(cl_max_mag):
    """
    Extract photometric data, and membership probabilities. Remove ID's to
    make entire array of floats.
    """
    mags_cols_cl = [[], []]
    for mag in list(zip(*list(zip(*cl_max_mag))[1:][2])):
        mags_cols_cl[0].append(mag)
    for col in list(zip(*list(zip(*cl_max_mag))[1:][4])):
        mags_cols_cl[1].append(col)

    # Store membership probabilities here.
    memb_probs = np.array(list(zip(*cl_max_mag))[1:][8])

    return mags_cols_cl, memb_probs


def main(cl_max_mag, lkl_method, lkl_binning, lkl_manual_bins):
    '''
    Prepare observed cluster array here to save time before the algorithm to
    find the best synthetic cluster fit is used.
    '''

    mags_cols_cl, memb_probs = dataProcess(cl_max_mag)

    if lkl_method in ['dolphin', 'mighell', 'tremmel']:
        # Obtain bin edges for each dimension, defining a grid.
        bin_edges = bin_edges_f(
            lkl_binning, mags_cols_cl, lkl_manual_bins=lkl_manual_bins)

        # Put all magnitudes and colors into a single list.
        obs_mags_cols = mags_cols_cl[0] + mags_cols_cl[1]
        # Obtain histogram for observed cluster.

        # Tried this (15/01/20) for RUP85 to see if I could get the fit
        # to not give that much weight to the low mass section.
        # idx = np.argsort(obs_mags_cols[0])
        # ww = np.zeros(len(obs_mags_cols[0]))
        # for i, w in enumerate(np.linspace(1., 0., len(obs_mags_cols[0]))):
        #     ww[idx[i]] = w
        # memb_probs = np.array(ww)

        cl_histo = np.histogramdd(obs_mags_cols, bins=bin_edges)[0]

        # DEPRECATED 10/01/20
        # Weight each bin by the 'lkl_weight' statistic (mean by default) of
        # the MPs of each star within that bin.
        # bin_w = np.nan_to_num(binned_statistic_dd(
        #     np.array(obs_mags_cols).T, memb_probs, statistic='mean',
        #     bins=bin_edges)[0])

        # Flatten N-dimensional histograms.
        cl_histo_f = np.array(cl_histo).ravel()
        # DEPRECATED 10/01/20
        # bin_weight_f = np.array(bin_w).ravel()

        mighell_ni_cnst = -2. * (cl_histo_f + np.minimum(cl_histo_f, 1.))

        # Index of bins where stars were observed. Used by the
        # 'Dolphin' and 'Mighell' likelihoods.
        cl_z_idx = (cl_histo_f != 0)

        # Remove all bins where n_i=0 (no observed stars). Used by the
        # 'Dolphin' likelihood.
        cl_histo_f_z = cl_histo_f[cl_z_idx]
        # DEPRECATED 10/01/20
        # bin_weight_f_z = bin_weight_f[cl_z_idx]

        # Dolphin n_i dependent constant.
        # N = cl_histo.sum()
        # n_i constant: 2 * [sum(n_i * ln(n_i)) - N] =
        dolphin_cst = 2. * (
            np.sum(cl_histo_f_z * np.log(cl_histo_f_z)) - cl_histo.sum())

        # Value to use when filling bins where n_i!=0 & m_i=0. This is an
        # empirical value that seems to work reasonably well. Low mass clusters
        # will have their masses underestimated, and I don't know how to
        # prevent this.
        fill_factor = min(.9, cl_histo_f.sum() / 1e4)

        if lkl_method == 'tremmel':
            obs_clust = [bin_edges, cl_histo_f_z, cl_z_idx]
        elif lkl_method == 'dolphin':
            obs_clust = [bin_edges, fill_factor, cl_histo_f_z, dolphin_cst,
                         cl_z_idx]
        elif lkl_method == 'mighell':
            obs_clust = [bin_edges, cl_histo_f, mighell_ni_cnst]

    elif lkl_method == 'tolstoy':

        # Square errors here to not repeat the same calculations each time a
        # new synthetic cluster is matched.
        e_mags_cols = []
        for e_m in list(zip(*list(zip(*cl_max_mag))[1:][3])):
            e_mags_cols.append(np.square(e_m))
        for e_c in list(zip(*list(zip(*cl_max_mag))[1:][5])):
            e_mags_cols.append(np.square(e_c))

        # Store and pass to use in likelihood function. The 'obs_st' list is
        # made up of:
        # obs_st = [star_1, star_2, ...]
        # star_i = [phot_1, phot_2, phot_3, ...]
        # phot_j = [phot_val, error]
        # Where 'phot_j' is a photometric dimension (magnitude or color), and
        # 'phot_val', 'error' the associated value and error for 'star_i'.
        obs_st = []
        mags_cols = mags_cols_cl[0] + mags_cols_cl[1]
        for st_phot, st_e_phot in list(
                zip(list(zip(*mags_cols)), list(zip(*e_mags_cols)))):
            obs_st.append(list(zip(*[st_phot, st_e_phot])))

        obs_clust = [[np.array(obs_st), len(obs_st), np.log(memb_probs)]]

        # New
        # shape: (Nstars, Ndims)
        obs_photom = np.array(mags_cols_cl[0] + mags_cols_cl[1]).T

        # Sum of squared photometric errors, for all dimensions. Clip at a
        # minimum of 0.005 to avoid numeric issues below.
        # sigma_sum = np.clip(
        #     obs_st[:, None, :, 1] + syn_st[None, :, :, 1], 0.005, None)

        sigma = np.clip(np.array(e_mags_cols).T, .005, None)

        # Product of summed squared sigmas.
        sigma_prod = np.prod(sigma, axis=-1)

        N = obs_photom.shape[0]
        obs_clust.append([
            obs_photom, sigma, sigma_prod, N, np.log(memb_probs)])

    elif lkl_method == 'isochfit':

        bin_edges = bin_edges_f(
            lkl_binning, mags_cols_cl, lkl_manual_bins=lkl_manual_bins)

        # Put all magnitudes and colors into a single list.
        obs_mags_cols = mags_cols_cl[0] + mags_cols_cl[1]
        # Obtain histogram for observed cluster.
        cl_histo = np.histogramdd(obs_mags_cols, bins=bin_edges)[0]
        # Flatten
        cl_histo = cl_histo.ravel()
        # Downsample
        # cl_histo[cl_histo > 5] = 5

        # idxs = np.array(np.where(cl_histo > 0)).T
        # coords = []
        # for pt in idxs:
        #     coord = []
        #     for i, j in enumerate(pt):
        #         coord.append(bin_edges[i][j])
        #     coords.append(coord)
        # coords = np.array(coords)

        obs_clust = [cl_histo, bin_edges]

    # DEPRECATED 17/01/20
    # elif lkl_method == 'duong':
    #     # Define variables to communicate with package 'R'.
    #     import rpy2.robjects as robjects
    #     from rpy2.robjects.packages import importr
    #     ks = importr('ks')
    #     kde_test = ks.kde_test
    #     hpi_kfe = ks.Hpi_kfe

    #     # CMD for cluster region.
    #     mags_cols = mags_cols_cl[0] + mags_cols_cl[1]
    #     matrix_cl = np.ravel(np.column_stack((mags_cols)))
    #     # matrix_cl = []
    #     # for st in obs_st:
    #     #     matrix_cl.append(st[0][0])
    #     #     matrix_cl.append(st[1][0])
    #     rows_cl = int(len(matrix_cl) / 2)

    #     # Create matrices for these CMDs.
    #     m_cl = robjects.r.matrix(robjects.FloatVector(matrix_cl),
    #                              nrow=rows_cl, byrow=True)

    #     # Bandwidth matrices.
    #     hpic = hpi_kfe(x=m_cl, binned=True)

    #     obs_clust = [kde_test, hpi_kfe, m_cl, hpic]

    elif lkl_method in ['dolphin_kde', 'kdeKL']:

        vals = np.array(mags_cols_cl[0] + mags_cols_cl[1])
        # vals.shape = (# of dims, # of data)
        kernel = gaussian_kde(vals)

        # TODO 'kdeKL' is still terribly slow. Also: is this a likelihood? Is
        # this a log-likelihood? Can I just plug in whatever I want here?
        # Being a KDE and thus normalized, does it fit masses? <--!!?

        # How should the number that defines the grid density be selected?
        Nb = 10
        pts = np.mgrid[[slice(A.min(), A.max(), complex(Nb)) for A in vals]]
        kde_pts = np.vstack([_.ravel() for _ in pts])
        obs_kde = kernel(kde_pts)

        obs_clust = [obs_kde, kde_pts]

    return obs_clust

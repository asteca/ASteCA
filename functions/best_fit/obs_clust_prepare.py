# -*- coding: utf-8 -*-
"""
Created on Tue Nov 25 16:00:00 2014

@author: gabriel
"""
import numpy as np
from astroML.plotting import hist
from .._in import get_in_params as g


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

    else:

        # Remove ID's (to make entire array of floats) and zip.
        P = np.array(zip(*memb_prob_avrg_sort)[1:], dtype='float')
        mag_col_cl = [P[4], P[2]]

        # Obtain bin edges for each dimension.
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

        # # Obtain *weighted* histogram for observed cluster.
        # cl_histo = np.histogramdd(
        #     cl_mags_cols, bins=bin_edges,
        #     weights=np.asarray(P[6]))[0]

        # print cl_histo, '\n'

        # Create extra bins at the edges to ensure that *all* stars in the
        # synthetic clusters generated are considered.

        print bin_edges, '\n'
        # Check that the edges are in fact below the [-50., 50] range for
        # each magnitude/color (which they definitely should).
        for i, phot_dimens_edge in enumerate(bin_edges):

            if phot_dimens_edge[0] > -30.:
                # Generate extra bin for left edge of this photometric
                # dimension.
                bin_edges[i] = np.insert(
                    bin_edges[i], 0, np.arange(-30., phot_dimens_edge[0], 0.5))
            if phot_dimens_edge[-1] < 10.:
                # Generate extra bin for right edge of this photometric
                # dimension.
                bin_edges[i] = np.append(
                    bin_edges[i], np.arange(phot_dimens_edge[-1], 10.,  0.5))

        print bin_edges
        # Zip magnitudes and colors into array.
        # cl_mags_cols = np.array(zip(*mag_col_cl))

        # Obtain *weighted* histogram for observed cluster.
        cl_histo = np.histogramdd(
            cl_mags_cols, range=bin_edges,
            weights=np.asarray(P[6]))[0]

        # print cl_histo2
        # import matplotlib.pyplot as plt
        # from matplotlib.ticker import MultipleLocator
        # b_rx, b_ry = bin_edges
        # fig = plt.figure()
        # ax1 = fig.add_subplot(1, 3, 1)
        # x_extend = [min(mag_col_cl[0]) - 1., max(mag_col_cl[0]) + 1.]
        # y_extend = [max(mag_col_cl[1]) + 5., min(mag_col_cl[1]) - 5.]
        # ax1.set_xlim(x_extend)
        # ax1.set_ylim(y_extend)
        # ax2 = fig.add_subplot(1, 3, 2)
        # ax2.minorticks_on()
        # ax2.yaxis.set_major_locator(MultipleLocator(1.0))
        # ax3 = fig.add_subplot(1, 3, 3)
        # ax3.minorticks_on()
        # ax3.yaxis.set_major_locator(MultipleLocator(1.0))
        # ax1.scatter(mag_col_cl[0], mag_col_cl[1], c='r', label='Obs clust')
        # for x_ed in b_rx:
        #     # vertical lines
        #     ax1.axvline(x_ed, linestyle=':', color='k', zorder=1)
        # for y_ed in b_ry:
        #     # horizontal lines
        #     ax1.axhline(y_ed, linestyle=':', color='k', zorder=1)
        # ax1.legend(fancybox=True, loc='lower left', scatterpoints=1)
        # ax2.imshow(cl_histo.transpose(), origin='lower', aspect='auto',
        #            interpolation="nearest", cmap=plt.cm.jet)
        # ax3.imshow(cl_histo2.transpose(), origin='lower', aspect='auto',
        #            interpolation="nearest", cmap=plt.cm.jet)
        # # ax1.invert_yaxis()
        # ax2.invert_yaxis()
        # ax3.invert_yaxis()
        # ax2.text(0.05, 0.95, 'Obs clust 2D histo', transform=ax2.transAxes,
        #          fontsize=15)
        # # Set limits for imshow.
        # # # ax2.set_xlim(*x_extend)
        # # # ax2.set_ylim(y_extend[0], y_extend[1])
        # # ax2.set_ylim(14., -4)
        # # # ax3.set_xlim(x_extend)
        # # # ax3.set_ylim(y_extend[0], y_extend[1])
        # # ax3.set_ylim(14., -4)
        # fig.subplots_adjust(hspace=1)
        # plt.show()

        # Pass observed cluster data.
        # obs_clust = [cl_histo, bin_edges]

        # Pass this list instead if plotting in get_likelihood.
        obs_clust = [cl_histo, bin_edges, mag_col_cl]

    return obs_clust

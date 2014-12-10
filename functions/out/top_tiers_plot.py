# -*- coding: utf-8 -*-
"""
Created on Wed Dic 10 12:00:00 2014

@author: gabriel
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from .._in import get_in_params as g


def plot_top_tiers():
    '''
    Produce output image for top tier models.
    '''

    # Plot all outputs
    # figsize(x1, y1), GridSpec(y2, x2) --> To have square plots: x1/x2 =
    # y1/y2 = 2.5
    fig = plt.figure(figsize=(20, 40))  # create the top-level container
    gs = gridspec.GridSpec(16, 8)  # create a GridSpec object

    # Best fitting process plots for GA.
    if bf_flag and best_fit_algor == 'genet':

        # Age vs metallicity GA diagram.
        ax20 = plt.subplot(gs[10:12, 0:2])
        # Axis limits.
        plt.xlim(m_min, m_max)
        plt.ylim(a_min, a_max)
        plt.xlabel('$z$', fontsize=16)
        plt.ylabel('$log(age)$', fontsize=16)
        plt.minorticks_on()
        # Plot best fit point.
        plt.scatter(m, a, marker='o', c='r', s=30)
        # Check if errors in both dimensions are defined.
        if all([i > 0. for i in [e_m, e_a]]):
            # Plot ellipse error.
            plt.gca()
            ellipse = Ellipse(xy=(m, a), width=2 * e_m, height=2 * e_a,
                                    edgecolor='r', fc='None', lw=1.)
            ax20.add_patch(ellipse)
        elif e_m < 0.:
            plt.errorbar(m, a, yerr=e_a, color='r')
        elif e_a < 0.:
            plt.errorbar(m, a, xerr=e_m, color='r')
        # Plot density map.
        hist, xedges, yedges = np.histogram2d(zip(*model_done[0])[0],
            zip(*model_done[0])[1], bins=100)
        # H_g is the 2D histogram with a gaussian filter applied
        h_g = gaussian_filter(hist, 2, mode='constant')
        plt.imshow(h_g.transpose(), origin='lower',
                   extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
                   cmap=plt.get_cmap('Blues'), aspect='auto')
        # Plot top tiers.
        top_tiers = [x for y, x in sorted(zip(model_done[1], model_done[0]))]
        print top_tiers[:10]
        top_t = [1, 2, 3, 4, 5, 6, 7, 8, 9]
        for i, txt in enumerate(top_t):
            ax20.scatter(zip(*top_tiers)[0][i], zip(*top_tiers)[1][i], s=1)
            ax20.annotate(txt, (zip(*top_tiers)[0][i], zip(*top_tiers)[1][i]),
                size=9)

    # Generate output file for each data file.
    fig.tight_layout()
    pl_fmt, pl_dpi = g.pl_params[1:3]
    plt.savefig(join(output_subdir, str(clust_name) + '.' + pl_fmt), dpi=pl_dpi)

    # Close to release memory.
    plt.clf()
    plt.close()
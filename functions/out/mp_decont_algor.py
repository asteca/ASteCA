# -*- coding: utf-8 -*-
"""
Created on Tue Dic 16 12:00:00 2014

@author: gabriel
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import FormatStrFormatter
from functions.exp_function import exp_3p
from .._in import get_in_params as g


def pl_mp_histo(gs, red_return, memb_prob_avrg_sort):
    '''
    Histogram for the distribution of membership probabilities from the
    decontamination algorithm.
    '''
    # Reduced membership.
    min_prob = red_return[1]
    ax = plt.subplot(gs[4:6, 2:4])
    plt.xlim(0., 1.)
    plt.xlabel('membership probability', fontsize=12)
    plt.ylabel('N', fontsize=12)
    ax.minorticks_on()
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
    prob_data = [star[7] for star in memb_prob_avrg_sort]
    # Histogram of the data.
    n, bins, patches = plt.hist(prob_data, 25, normed=1)
    # Get bin centers.
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    # scale values to interval [0,1]
    col = bin_centers - min(bin_centers)
    col /= max(col)
    cm = plt.cm.get_cmap('RdYlBu_r')
    # Plot histo colored according to colormap.
    for c, p in zip(col, patches):
        plt.setp(p, 'facecolor', cm(c))
    # Plot minimum probability value used in best isochrone fit function, if
    # it was used.
    bf_flag = g.bf_params[0]
    if bf_flag:
        text = '$prob_{min}=%.2f$' % min_prob
        plt.text(0.05, 0.92, text, transform=ax.transAxes,
             bbox=dict(facecolor='white', alpha=0.85), fontsize=14)
        plt.axvline(x=min_prob, linestyle='--', color='green', lw=2.5)


def pl_chart_mps(gs, fig, x_name, y_name, coord, x_zmin, x_zmax, y_zmin,
    y_zmax, center_cl, clust_rad, field_dens, memb_prob_avrg_sort, stars_out):
    '''
    Finding chart of cluster region with decontamination algorithm
    applied and colors assigned according to the probabilities obtained.
    '''
    ax = plt.subplot(gs[4:6, 4:6])
    #Set plot limits, Use 'zoom' x,y ranges.
    if coord == 'deg':
        # If RA is used, invert axis.
        plt.xlim(x_zmax, x_zmin)
        # Set x axis to not use scientific notation.
        ax.xaxis.set_major_formatter(FormatStrFormatter('%0.1f'))
    else:
        plt.xlim(x_zmin, x_zmax)
    plt.ylim(y_zmin, y_zmax)
    #Set axis labels
    plt.xlabel('{} ({})'.format(x_name, coord), fontsize=12)
    plt.ylabel('{} ({})'.format(y_name, coord), fontsize=12)
    # Set minor ticks
    ax.minorticks_on()
    # Radius
    circle = plt.Circle((center_cl[0], center_cl[1]), clust_rad,
                        color='red', fill=False)
    fig.gca().add_artist(circle)
    plt.text(0.63, 0.93, 'Cluster region', transform=ax.transAxes,
        bbox=dict(facecolor='white', alpha=0.8), fontsize=12)
    # Color map, higher prob stars look redder.
    cm = plt.cm.get_cmap('RdYlBu_r')
    # Star sizes for dense and not dense regions.
    st_size = 20 if field_dens > 0.005 else 35
    m_p_m_temp = [[], [], []]
    for star in memb_prob_avrg_sort:
        m_p_m_temp[0].append(star[1])
        m_p_m_temp[1].append(star[2])
        m_p_m_temp[2].append(star[7])
    # Create new list with inverted values so higher prob stars are on top.
    m_p_m_temp_inv = [i[::-1] for i in m_p_m_temp]
    plt.scatter(m_p_m_temp_inv[0], m_p_m_temp_inv[1], marker='o',
                c=m_p_m_temp_inv[2], s=st_size, edgecolors='black',
                cmap=cm, lw=0.5)
    out_clust_rad = [[], []]
    for star in stars_out:
        if x_zmin <= star[1] <= x_zmax and y_zmin <= star[2] <= y_zmax:
            dist = np.sqrt((center_cl[0] - star[1]) ** 2 +
            (center_cl[1] - star[2]) ** 2)
            # Only plot stars outside the cluster's radius.
            if dist >= clust_rad:
                out_clust_rad[0].append(star[1])
                out_clust_rad[1].append(star[2])
    plt.scatter(out_clust_rad[0], out_clust_rad[1], marker='o',
                s=st_size, edgecolors='black', facecolors='none', lw=0.5)


def pl_mps_phot_diag(gs, fig, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd,
    x_ax, y_ax, memb_prob_avrg_sort, shift_isoch, col_data, err_lst):
    '''
    Star's membership probabilities on cluster's photom diagram.
    '''
    ax = plt.subplot(gs[4:6, 6:8])
    #Set plot limits
    plt.xlim(x_min_cmd, x_max_cmd)
    plt.ylim(y_min_cmd, y_max_cmd)
    #Set axis labels
    plt.xlabel('$' + x_ax + '$', fontsize=18)
    plt.ylabel('$' + y_ax + '$', fontsize=18)
    tot_clust = len(memb_prob_avrg_sort)
    text = '$r \leq r_{cl}\,|\,N=%d$' % tot_clust
    plt.text(0.5, 0.93, text, transform=ax.transAxes,
             bbox=dict(facecolor='white', alpha=0.5), fontsize=16)
    # Set minor ticks
    ax.minorticks_on()
    ax.xaxis.set_major_locator(MultipleLocator(1.0))
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
    # This reversed colormap means higher prob stars will look redder.
    cm = plt.cm.get_cmap('RdYlBu_r')
    m_p_m_temp = [[], [], []]
    for star in memb_prob_avrg_sort:
        m_p_m_temp[0].append(star[5])
        m_p_m_temp[1].append(star[3])
        m_p_m_temp[2].append(star[7])
    # Create new list with inverted values so higher prob stars are on top.
    m_p_m_temp_inv = [i[::-1] for i in m_p_m_temp]
    v_min_mp, v_max_mp = round(min(m_p_m_temp[2]), 2), \
    round(max(m_p_m_temp[2]), 2)
    if v_min_mp != v_max_mp:
        col_select, c_iso = m_p_m_temp_inv[2], 'g'
    else:
        col_select, c_iso = '#4682b4', 'r'
    # Plot isochrone if best fit process was used.
    bf_flag = g.bf_params[0]
    if bf_flag:
        plt.plot(shift_isoch[0], shift_isoch[1], c=c_iso, lw=1.2)
    sca = plt.scatter(m_p_m_temp_inv[0], m_p_m_temp_inv[1], marker='o',
                c=col_select, s=40, cmap=cm, lw=0.5, vmin=v_min_mp,
                vmax=v_max_mp)
    # If list is not empty.
    if m_p_m_temp_inv[1]:
        # Plot error bars at several values.
        mag_y = np.arange(int(min(m_p_m_temp_inv[1]) + 0.5),
                          int(max(m_p_m_temp_inv[1]) + 0.5) + 0.1)
        x_val = [min(x_max_cmd, max(col_data) + 0.2) - 0.4] * len(mag_y)
        # Read average fitted values for exponential error fit.
        popt_mag, popt_col1 = err_lst[:2]
        plt.errorbar(x_val, mag_y, yerr=exp_3p(mag_y, *popt_mag),
                     xerr=exp_3p(mag_y, *popt_col1), fmt='k.', lw=0.8,
                     ms=0., zorder=4)
    # Plot colorbar (see bottom of make_plots file).
    trans = ax.transAxes + fig.transFigure.inverted()
    plot_colorbar = False
    if v_min_mp != v_max_mp:
        plot_colorbar = True

    return plot_colorbar, v_min_mp, v_max_mp, sca, trans


def plot(N, *args):
    '''
    Handle each plot separately.
    '''

    plt_map = {
        1: [pl_mp_histo, 'MPs histogram'],
        2: [pl_chart_mps, 'frame with MPs coloring']
    }

    fxn = plt_map.get(N, None)[0]
    if fxn is None:
        raise ValueError("  ERROR: there is no plot {}.".format(N))

    try:
        fxn(*args)
    except:
        print("  WARNING: error when plotting {}.".format(plt_map.get(N)[1]))
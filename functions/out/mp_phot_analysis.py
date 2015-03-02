# -*- coding: utf-8 -*-
"""
Created on Tue Dic 16 12:00:00 2014

@author: gabriel
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from .._in import get_in_params as g
from functions.exp_function import exp_3p


def phot_diag_st_size(x):
    '''
    Calculate optimal size for stars in photometric diagram.
    '''
    a, b, c, d = 2.99, -2.81, 563.36, 15.02
    return ((a - d) / (1 + ((x / c) ** b))) + d


def pl_phot_err(gs, fig, up_low, x_ax, y_ax, mag_data, err_plot, err_flags,
    cl_region, stars_in_rjct, stars_out, stars_out_rjct):
    '''
    Photometric error rejection.
    '''
    # Define parameters for upper and lower plots.
    if up_low == 'up':
        ax, ax_y, j = plt.subplot(gs[2, 2:4]), y_ax, 4
    else:
        ax, ax_y, j = plt.subplot(gs[3, 2:4]), x_ax, 6

    # Error parameters.
    er_mode, e_max, be, be_e, N_sig = g.er_params
    err_all_fallback, err_max_fallback = err_flags

    #Set plot limits
    x_min, x_max = min(mag_data) - 0.5, max(mag_data) + 0.5
    plt.xlim(x_min, x_max)
    plt.ylim(-0.005, e_max + (e_max / 5.))
    #Set axis labels
    plt.ylabel('$\sigma_{' + ax_y + '}$', fontsize=18)
    plt.xlabel('$' + y_ax + '$', fontsize=18)
    # Set minor ticks
    ax.minorticks_on()
    # Plot e_max line.
    ax.hlines(y=e_max, xmin=x_min, xmax=x_max, color='k',
        linestyles='dashed', zorder=2)
    # Plot rectangle.
    bright_end = min(mag_data) + be
    ax.vlines(x=bright_end + 0.05, ymin=-0.005, ymax=be_e, color='k',
               linestyles='dashed', zorder=2)
    ax.vlines(x=min(mag_data) - 0.05, ymin=-0.005, ymax=be_e, color='k',
               linestyles='dashed', zorder=2)
    ax.hlines(y=be_e, xmin=min(mag_data), xmax=bright_end, color='k',
               linestyles='dashed', zorder=2)
    # If any method could be used.
    if err_all_fallback is False and err_max_fallback is False:
        # Plot curve(s) according to the method used.
        if er_mode == 'eyefit':
            # Unpack params.
            popt_mag, pol_mag, popt_col, pol_col, mag_val_left, \
            mag_val_right, col_val_left, col_val_right = err_plot

            if up_low == 'up':
                val_left, popt_left = mag_val_left, popt_mag
                val_right, pol_right = mag_val_right, pol_mag
            else:
                val_left, popt_left = col_val_left, popt_col
                val_right, pol_right = col_val_right, pol_col

            # Plot left side of upper envelope (exponential).
            ax.plot(val_left, exp_3p(val_left, *popt_left),
                'k--', lw=2., zorder=3)
            # Plot right side of upper envelope (polynomial).
            ax.plot(val_right, np.polyval(pol_right, (val_right)),
                'k--', lw=2., zorder=3)
        elif er_mode == 'lowexp':
            mag_x = np.linspace(bright_end, max(mag_data), 50)
            # Unpack params.
            popt_mag, popt_col1 = err_plot
            # Plot exponential curve.
            ax.plot(mag_x, exp_3p(mag_x, *popt_mag), 'k-', zorder=3)

    # Plot rejected stars.
    if len(stars_out_rjct) > 0:
        # Only attempt to pot if any star is stored in the list.
        plt.scatter(zip(*stars_out_rjct)[3], zip(*stars_out_rjct)[j],
            marker='x', c='teal', s=15, zorder=1)
    if len(stars_in_rjct) > 0:
        plt.scatter(zip(*stars_in_rjct)[3], zip(*stars_in_rjct)[j],
            marker='x', c='teal', s=15, zorder=1)
    # Plot accepted stars.
    plt.scatter(zip(*stars_out)[3], zip(*stars_out)[j], marker='o', c='b',
        s=5, zorder=2, lw=0.3, label='$r > r_{cl}$')
    plt.scatter(zip(*cl_region)[3], zip(*cl_region)[j], marker='o', c='r',
        s=10, zorder=2, lw=0.6, label='$r \leq r_{cl}}$')
    if up_low == 'up':
        # Legends.
        leg = plt.legend(fancybox=True, loc='upper left', scatterpoints=1,
            fontsize=16, markerscale=2.5, prop={'size': 13})
        # Set the alpha value of the legend.
        leg.get_frame().set_alpha(0.7)


def pl_fl_diag(gs, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd, x_ax, y_ax,
    stars_out_rjct, field_regions):
    '''
    Field stars diagram (stars outside cluster's radius)
    '''
    ax = plt.subplot(gs[2:4, 4:6])
    #Set plot limits
    plt.xlim(x_min_cmd, x_max_cmd)
    plt.ylim(y_min_cmd, y_max_cmd)
    #Set axis labels
    plt.xlabel('$' + x_ax + '$', fontsize=18)
    plt.ylabel('$' + y_ax + '$', fontsize=18)
    # Set minor ticks
    ax.minorticks_on()
    # Only draw units on axis (ie: 1, 2, 3)
    ax.xaxis.set_major_locator(MultipleLocator(1.0))
    # Set grid
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
    # Plot rejected stars.
    # Outside of the cluster region.
    stars_rjct_temp = [[], []]
    for star in stars_out_rjct:
        stars_rjct_temp[0].append(star[5])
        stars_rjct_temp[1].append(star[3])
    plt.scatter(stars_rjct_temp[0], stars_rjct_temp[1], marker='x',
        c='teal', s=15, zorder=1)
    # Plot stars within the field regions defined, if at least one exists.
    if field_regions:
        stars_acpt_temp = [[], []]
        for fr in field_regions:
            for star in fr:
                stars_acpt_temp[0].append(star[5])
                stars_acpt_temp[1].append(star[3])
        plt.text(0.56, 0.93, '$\star_{field}\,|\,N=%d$' %
            len(stars_acpt_temp[0]), transform=ax.transAxes,
            bbox=dict(facecolor='white', alpha=0.5), fontsize=16)
        sz_pt = phot_diag_st_size(len(stars_acpt_temp[0]))
        plt.scatter(stars_acpt_temp[0], stars_acpt_temp[1], marker='o', c='b',
                    s=sz_pt, lw=0.3, zorder=2)


def pl_cl_diag(gs, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd, x_ax, y_ax,
    stars_in_rjct, cl_region, n_memb):
    '''
    Cluster's stars diagram (stars inside cluster's radius)
    '''
    ax = plt.subplot(gs[2:4, 6:8])
    #Set plot limits
    plt.xlim(x_min_cmd, x_max_cmd)
    plt.ylim(y_min_cmd, y_max_cmd)
    #Set axis labels
    plt.xlabel('$' + x_ax + '$', fontsize=18)
    plt.ylabel('$' + y_ax + '$', fontsize=18)
    # Set minor ticks
    ax.minorticks_on()
    # Only draw units on axis (ie: 1, 2, 3)
    ax.xaxis.set_major_locator(MultipleLocator(1.0))
    # Set grid
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
    # Calculate total number of stars whitin cluster's radius.
    tot_stars = len(stars_in_rjct) + len(cl_region)
    text1 = '$r \leq r_{{cl}}\,|\,N={}$'.format(tot_stars)
    text2 = r'$n_{{memb}} \approx {}$'.format(n_memb)
    text = text1 + '\n' + text2
    plt.text(0.55, 0.87, text, transform=ax.transAxes,
             bbox=dict(facecolor='white', alpha=0.5), fontsize=16)
    # Plot stars in CMD.
    if len(stars_in_rjct) > 0:
        # Only attempt to pot if any star is stored in the list.
        plt.scatter(zip(*stars_in_rjct)[5], zip(*stars_in_rjct)[3],
            marker='x', c='teal', s=12, zorder=1)
    sz_pt = phot_diag_st_size(tot_stars)
    plt.scatter(zip(*cl_region)[5], zip(*cl_region)[3], marker='o', c='r',
                s=sz_pt, lw=0.3, zorder=2)


def pl_lum_func(gs, mag_data, y_ax, x_cl, y_cl, flag_area_stronger, x_fl,
    y_fl, completeness):
    '''
    LF of stars in cluster region and outside.
    '''
    ax = plt.subplot(gs[2:4, 8:10])
    #Set plot limits
    x_min, x_max = min(mag_data) - 0.5, max(mag_data) + 0.5
    plt.xlim(x_max, x_min)
    ax.minorticks_on()
    # Only draw units on axis (ie: 1, 2, 3)
    ax.xaxis.set_major_locator(MultipleLocator(2.0))
    # Set grid
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
    #Set axis labels
    plt.xlabel('$' + y_ax + '$', fontsize=18)
    plt.ylabel('$N^{\star}/A_{cl}$', fontsize=18)
    # Cluster region LF (contaminated).
    plt.step(x_cl, y_cl, where='post', color='r', ls='--', lw=1.5,
        label='$LF_{cl+fl} \;(r \leq r_{cl})$', zorder=2)
    # Check if field regiones were defined.
    if flag_area_stronger is not True:
        # Average field regions LF.
        plt.step(x_fl, y_fl, where='post', color='b', ls='--', lw=1.5,
            label='$LF_{fl} \;(r  > r_{cl})$', zorder=3)
        # Cluster region LF - average field regions LF.
        plt.step(x_cl, y_cl - y_fl, where='post', color='g', lw=1.7,
            label='$LF_{cl}$', zorder=4)
        # Force y axis min to 0.
        max_y = max(max(y_cl), max(y_fl))
    else:
        # Force y axis min to 0.
        max_y = max(y_cl)
    plt.ylim(0., max_y + 0.05 * max_y)
    # Completeness maximum value.
    # completeness = [max_mag, bin_edges, max_indx, comp_perc]
    bin_edges, max_indx = completeness[1], completeness[2]
    mag_peak = bin_edges[max_indx]
    text = '$' + y_ax + r',_{compl}\,\approx\,%0.1f$' % mag_peak
    ax.vlines(x=mag_peak, ymin=0., ymax=plt.ylim()[1], color='k',
        lw=1.5, linestyles='dashed', label=text, zorder=1)
    # Legends.
    leg = plt.legend(fancybox=True, loc='upper right', numpoints=1,
                       fontsize=13)
    # Set the alpha value of the legend.
    leg.get_frame().set_alpha(0.7)


def pl_integ_mag(gs, integr_return, y_ax, x_ax0, flag_area_stronger):
    '''
    Integrated magnitudes.
    '''
    if integr_return:
        # Unpack values.
        m_ord = g.axes_params[2]
        cl_reg_mag1, fl_reg_mag1, integ_mag1, cl_reg_mag2, fl_reg_mag2, \
        integ_mag2 = integr_return
        # Make plot
        ax = plt.subplot(gs[2:4, 10:12])
        # If field lists are not empty.
        if fl_reg_mag1[0].any() and fl_reg_mag2[0].any():
            x_min = min(min(cl_reg_mag1[0]), min(fl_reg_mag1[0]),
                min(cl_reg_mag2[0]), min(fl_reg_mag2[0])) - 0.2
            x_max = max(max(cl_reg_mag1[0]), max(fl_reg_mag1[0]),
                max(cl_reg_mag2[0]), max(fl_reg_mag2[0])) + 0.2
            y_min = max(max(cl_reg_mag1[1]), max(fl_reg_mag1[1]),
                max(cl_reg_mag2[1]), max(fl_reg_mag2[1])) + 0.2
            y_max = min(min(cl_reg_mag1[1]), min(fl_reg_mag1[1]),
                min(cl_reg_mag2[1]), min(fl_reg_mag2[1])) - 0.2
        else:
            x_min, x_max = min(min(cl_reg_mag1[0]),
                min(cl_reg_mag2[0])) - 0.2, max(max(cl_reg_mag1[0]),
                max(cl_reg_mag2[0])) + 0.2
            y_min, y_max = max(max(cl_reg_mag1[1]),
                max(cl_reg_mag2[1])) + 0.2, min(min(cl_reg_mag1[1]),
                min(cl_reg_mag2[1])) - 0.2
        plt.xlim(x_min, x_max)
        plt.ylim(y_min, y_max)
        ax.set_xlabel('$mag$', fontsize=18)
        ax.set_ylabel('$mag^*$', fontsize=18)
        ax.minorticks_on()
        ax.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
        text1 = '$' + y_ax + '^{*}_{cl+fl}$'
        text2 = '$' + x_ax0 + '^{*}_{cl+fl}$'
        # Cluster + field integrated magnitude curve.
        plt.plot(cl_reg_mag1[0], cl_reg_mag1[1], 'r-', lw=1., label=text1)
        # Cluster integrated magnitude.
        plt.plot(cl_reg_mag2[0], cl_reg_mag2[1], 'r:', lw=2., label=text2)
        # Check if field regiones were defined.
        if not flag_area_stronger:
            text3 = '$' + y_ax + '^{*}_{fl}$'
            text4 = '$' + x_ax0 + '^{*}_{fl}$'
            # Field average integrated magnitude curve.
            plt.plot(fl_reg_mag1[0], fl_reg_mag1[1], 'b-', lw=1.,
                label=text3)
            # Field average integrated magnitude.
            plt.plot(fl_reg_mag2[0], fl_reg_mag2[1], 'b:', lw=2.,
                label=text4)
        # Check how the second magnitude whould be formed.
        if m_ord == 21:
            sig, text0 = 1., x_ax0 + '^{*} -' + y_ax
        elif m_ord == 12:
            sig, text0 = -1., y_ax + '^{*} -' + x_ax0
        int_col = sig * (integ_mag2 - integ_mag1)
        text = '$(' + text0 + '^{*} )_{cl} = %0.2f$' % int_col
        plt.text(0.22, 0.15, text, transform=ax.transAxes,
             bbox=dict(facecolor='white', alpha=0.75), fontsize=13)
        lines, labels = ax.get_legend_handles_labels()
        leg = ax.legend(lines, labels, loc='lower right', numpoints=1,
            fontsize=13)
        leg.get_frame().set_alpha(0.75)


def pl_p_vals(gs, flag_pval_test, pval_test_params):
    '''
    Distribution of KDE p_values.
    '''
    if flag_pval_test:
        # Extract parameters from list.
        prob_cl_kde, kde_cl_1d, kde_f_1d, x_kde, y_over = pval_test_params
        ax = plt.subplot(gs[4:6, 0:2])
        plt.xlim(-0.5, 1.5)
        plt.ylim(0, 1.02)
        plt.xlabel('p-values', fontsize=12)
        plt.ylabel('Density (normalized)', fontsize=12)
        ax.minorticks_on()
        ax.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
        # Grid to background.
        ax.set_axisbelow(True)
        # Plot field vs field KDE.
        if kde_f_1d.any():
            max_kde = max(max(kde_f_1d), max(kde_cl_1d))
            plt.plot(x_kde, kde_f_1d / max_kde, c='b', ls='-', lw=1.,
                label='$KDE_{fl}$')
        else:
            max_kde = max(kde_cl_1d)
        # Plot cluster vs field KDE.
        plt.plot(x_kde, kde_cl_1d / max_kde, c='r', ls='-', lw=1.,
            label='$KDE_{cl}$')
        # Fill overlap.
        if y_over:
            plt.fill_between(x_kde, np.asarray(y_over) / max_kde, 0,
                color='grey', alpha='0.5')
        text = '$P_{cl}^{KDE} = %0.2f$' % round(prob_cl_kde, 2)
        plt.text(0.05, 0.92, text, transform=ax.transAxes,
             bbox=dict(facecolor='white', alpha=0.6), fontsize=12)
        # Legend.
        handles, labels = ax.get_legend_handles_labels()
        leg = ax.legend(handles, labels, loc='upper right', numpoints=1,
                          fontsize=12)
        leg.get_frame().set_alpha(0.6)


def plot(N, *args):
    '''
    Handle each plot separately.
    '''

    plt_map = {
        1: [pl_phot_err, 'upper error rejection function'],
        2: [pl_phot_err, 'lower error rejection function'],
        3: [pl_fl_diag, 'field regions photometric diagram'],
        4: [pl_cl_diag, 'cluster region photometric diagram'],
        5: [pl_lum_func, 'luminosity function'],
        6: [pl_integ_mag, 'integrated magnitudes'],
        7: [pl_p_vals, 'KDE p-values distribution']
    }

    fxn = plt_map.get(N, None)[0]
    if fxn is None:
        raise ValueError("  ERROR: there is no plot {}.".format(N))

    try:
        fxn(*args)
    except:
        print("  WARNING: error when plotting {}.".format(plt_map.get(N)[1]))

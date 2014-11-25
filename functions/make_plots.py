"""
@author: gabriel
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle
import matplotlib.offsetbox as offsetbox
from matplotlib.ticker import MultipleLocator
from itertools import cycle
from scipy.ndimage.filters import gaussian_filter
from matplotlib.patches import Ellipse
from os.path import join
import warnings
# Custom functions.
from functions.exp_function import exp_func
import error_round as err_r


def star_size(mag_data):
    '''
    Convert magnitudes into intensities and define sizes of stars in
    finding chart.
    '''
    return 0.1 + 100. * 10 ** ((np.array(mag_data) - min(mag_data)) / -2.5)


def line(x, slope, intercept):
    '''
    Linar function.
    '''
    y = slope * x + intercept
    return y


def two_params(x, cd, rc, bg):
    '''
    Two parameters King profile fit.
    '''
    return bg + cd / (1 + (np.asarray(x) / rc) ** 2)


def three_params(x, rt, cd, rc, bg):
    '''
    Three parameters King profile fit.
    '''
    return cd * (1 / np.sqrt(1 + (np.asarray(x) / rc) ** 2) -
        1 / np.sqrt(1 + (rt / rc) ** 2)) ** 2 + bg


def reject_outliers(data, m=6.5):
    '''
    Reject outliers from array.
    http://stackoverflow.com/a/16562028/1391441
    '''
    d = np.abs(data - np.median(data))
    mdev = np.median(d)
    s = d / mdev if mdev else 0.
    return data[s < m]


def make_plots(output_subdir, clust_name, x_data, y_data, gd_params,
    bin_width, center_params, rdp_params, field_dens, radius_params,
    cont_index, mag_data, col1_data, err_plot, err_flags, kp_params,
    cl_region, stars_out, stars_in_rjct, stars_out_rjct, integr_return, n_memb,
    flag_area_stronger, field_regions, flag_pval_test,
    pval_test_params, memb_prob_avrg_sort, lum_func, completeness, ip_list,
    da_params, bf_params, red_return, err_lst, bf_return, ga_params, er_params,
    axes_params, pl_params):
    '''
    Make all plots.
    '''

    # Define names for CMD axes.
    y_axis = 0
    y_ax, x_ax0, m_ord = axes_params[0:3]
    if m_ord == 21:
        x_ax = '(' + x_ax0 + '-' + y_ax + ')'
    elif m_ord == 12:
        x_ax = '(' + y_ax + '-' + x_ax0 + ')'

    # Unpack coordinates and photometric data.
    #x_data, y_data = id_coords[1:]
    phot_x = col1_data
    phot_y = mag_data

    # Define system of coordinates used.
    px_deg = gd_params[-1]
    coord_lst = ['px', 'x', 'y'] if px_deg == 'px' else ['deg', 'ra', 'dec']
    coord, x_name, y_name = coord_lst

    # Define plot limits for *all* CMD diagrams.
    phot_x_s, phot_y_s = reject_outliers(phot_x), reject_outliers(phot_y)
    x_max_cmd, x_min_cmd = max(phot_x_s) + 0.5, min(phot_x_s) - 0.5
    y_min_cmd, y_max_cmd = max(phot_y_s) + 0.5, min(phot_y_s) - 0.5
    # If photometric axis y is a magnitude, make sure the brightest stars
    # are plotted.
    y_max_cmd = (min(phot_y) - 1.) if y_axis == 0 else y_max_cmd

    # Unpack params.
    # Parameters from get_center function.
    cent_bin, kde_cent, e_cent, approx_cents, st_dev_lst, hist_2d_g, \
    kde_pl = center_params[:7]
    center_cl = kde_cent
    # RDP params.
    radii, ring_density, poisson_error = rdp_params[:3]
    # Parameters from get_radius function.
    clust_rad, e_rad = radius_params[:2]
    # King prof params.
    rc, e_rc, rt, e_rt, n_c_k, cd, flag_2pk_conver, flag_3pk_conver = kp_params
    # Error parameters.
    er_mode, e_max, be, be_e, N_sig = er_params
    err_all_fallback, err_max_fallback = err_flags
    # Luminosity functions.
    x_cl, y_cl, x_fl, y_fl = lum_func
    # Reduced membership.
    min_prob = red_return[1]
    # Best isochrone fit params.
    bf_flag, best_fit_algor, N_b = bf_params
    # Genetic algorithm params.
    n_pop, n_gen, fdif, p_cross, cr_sel, p_mut, n_el, n_ei, n_es = ga_params
    # Best fitting process results.
    isoch_fit_params, isoch_fit_errors, shift_isoch, synth_clst = bf_return

    # Plot all outputs
    # figsize(x1, y1), GridSpec(y2, x2) --> To have square plots: x1/x2 =
    # y1/y2 = 2.5
    fig = plt.figure(figsize=(20, 40))  # create the top-level container
    gs1 = gridspec.GridSpec(16, 8)  # create a GridSpec object

    # 2D gaussian convolved histogram.
    ax0 = plt.subplot(gs1[0:2, 0:2])
    plt.xlabel('{} (bins)'.format(x_name), fontsize=12)
    plt.ylabel('{} (bins)'.format(y_name), fontsize=12)
    ax0.minorticks_on()
    plt.axvline(x=cent_bin[0], linestyle='--', color='white')
    plt.axhline(y=cent_bin[1], linestyle='--', color='white')
    # Radius
    circle = plt.Circle((cent_bin[0], cent_bin[1]),
        clust_rad / bin_width, color='w', fill=False)
    fig.gca().add_artist(circle)
    # Add text boxs.
    bin_w_r = err_r.round_to_y(bin_width)
    text = 'Bin $\simeq$ {0:g} {1}'.format(bin_w_r, coord)
    ob = offsetbox.AnchoredText(text, loc=1, prop=dict(size=10))
    ob.patch.set(boxstyle='square,pad=-0.2', alpha=0.85)
    ax0.add_artist(ob)
    plt.imshow(hist_2d_g.transpose(), origin='lower')

    # 2D Gaussian histograms' centers using different standard deviations.
    ax1 = plt.subplot(gs1[0:2, 2:4])
    # Get max and min values in x,y
    x_min, x_max = min(x_data), max(x_data)
    y_min, y_max = min(y_data), max(y_data)
    #Set plot limits
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    plt.xlabel('{} ({})'.format(x_name, coord), fontsize=12)
    plt.ylabel('{} ({})'.format(y_name, coord), fontsize=12)
    ax1.minorticks_on()
    ## Add lines through median values with std deviations.
    cent_median, cent_std_dev = np.mean(np.array(approx_cents), axis=0), \
    np.std(np.array(approx_cents), axis=0)
    plt.axvline(x=cent_median[0], linestyle='-', color='k')
    plt.axvline(x=cent_median[0] + cent_std_dev[0], linestyle='--',
        color='k')
    plt.axvline(x=cent_median[0] - cent_std_dev[0], linestyle='--',
        color='k')
    plt.axhline(y=cent_median[1], linestyle='-', color='k')
    plt.axhline(y=cent_median[1] + cent_std_dev[1], linestyle='--',
        color='k')
    plt.axhline(y=cent_median[1] - cent_std_dev[1], linestyle='--',
        color='k')
    # Add stats box.
    cent_med_r, cent_std_dev_r = err_r.round_sig_fig(cent_median, cent_std_dev)
    text1 = (r"$(\tilde{{{0}}},\tilde{{{1}}}) = ({3:g},\;{4:g})\,"
    "{2}$").format(x_name, y_name, coord, *cent_med_r)
    text2 = ("$(\sigma_{{{0}}},\sigma_{{{1}}}) = ({3:g},\;{4:g})\,"
    "{2}$").format(x_name, y_name, coord, *cent_std_dev_r)
    text = text1 + '\n' + text2
    plt.text(0.05, 0.9, text, transform=ax1.transAxes,
        bbox=dict(facecolor='white', alpha=0.8), fontsize=11)
    # Plot centers.
    cols = cycle(['red', 'blue', 'green', 'black', 'cyan'])
    for i, center in enumerate(approx_cents):
        boxes = plt.gca()
        length = (bin_width * st_dev_lst[i]) * 2.
        boxes.add_patch(Rectangle(
            (center[0] - (length / 2.), center[1] - (length / 2.)),
            length, length,
            facecolor='none', edgecolor=next(cols), ls='solid',
            lw=1.5, zorder=(len(st_dev_lst) - i),
            label='St dev: %.1f' % st_dev_lst[i]))
    # get handles
    handles, labels = ax1.get_legend_handles_labels()
    # use them in the legend
    leg1 = ax1.legend(handles, labels, loc='lower right', numpoints=1,
        fontsize=7)
    # Set the alpha value of the legend.
    leg1.get_frame().set_alpha(0.5)

    # x,y finding chart of full frame
    ax4 = plt.subplot(gs1[2:4, 0:2])
    #Set plot limits
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    #Set axis labels
    plt.xlabel('{} ({})'.format(x_name, coord), fontsize=12)
    plt.ylabel('{} ({})'.format(y_name, coord), fontsize=12)
    # Set minor ticks
    ax4.minorticks_on()
    # Plot r_cl.
    circle = plt.Circle((center_cl[0], center_cl[1]), clust_rad, color='r',
        fill=False, lw=2.5)
    fig.gca().add_artist(circle)
    if flag_3pk_conver is True:
        # Plot tidal radius.
        circle = plt.Circle((center_cl[0], center_cl[1]), rt, color='g',
            fill=False, lw=2.5)
        fig.gca().add_artist(circle)
        # Plot core radius.
        if rc > 0:
            circle = plt.Circle((center_cl[0], center_cl[1]), rc,
                color='g', fill=False, ls='dashed', lw=2.5)
            fig.gca().add_artist(circle)
    elif flag_2pk_conver is True:
        # Plot core radius.
        if rc > 0:
            circle = plt.Circle((center_cl[0], center_cl[1]), rc,
                color='g', fill=False, ls='dashed', lw=2.5)
            fig.gca().add_artist(circle)
    # Add text box
    center_cl_r, e_cent_r = err_r.round_sig_fig(center_cl, e_cent)
    text1 = '${0}_{{cent}} = {1:g} \pm {2:g}\,{3}$'.format(x_name,
        center_cl_r[0], e_cent_r[0], coord)
    text2 = '${0}_{{cent}} = {1:g} \pm {2:g}\,{3}$'.format(y_name,
        center_cl_r[1], e_cent_r[1], coord)
    text = text1 + '\n' + text2
    plt.text(0.05, 0.9, text, transform=ax4.transAxes,
        bbox=dict(facecolor='white', alpha=0.85), fontsize=11)
    # Plot stars.
    st_sizes_arr = star_size(mag_data)
    plt.scatter(x_data, y_data, marker='o', c='black', s=st_sizes_arr)

    # Radial density plot.
    ax5 = plt.subplot(gs1[2:4, 2:6])
    # Get max and min values in x,y
    x_min, x_max = max(min(radii) - (max(radii) / 20.), 0), \
    max(radii) + (max(radii) / 20.)
    delta_total = (max(ring_density) - field_dens)
    delta_backg = 0.2 * delta_total
    y_min, y_max = max((field_dens - delta_backg) - (max(ring_density) -
    min(ring_density)) / 10, 0), max(ring_density) + (max(ring_density) -
    min(ring_density)) / 10
    # Set plot limits
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    # Set axes labels
    plt.xlabel('radius ({})'.format(coord), fontsize=12)
    plt.ylabel("stars/{}$^{{2}}$".format(coord), fontsize=12)
    # Set grid
    ax5.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
    # Cluster's name.
    text = str(clust_name)
    plt.text(0.4, 0.9, text, transform=ax5.transAxes, fontsize=14)
    # Round radii values.
    rads_r, e_rads_r = err_r.round_sig_fig([rc, rt, clust_rad],
        [e_rc, e_rt, e_rad])
    # Legend texts
    kp_text = '3P' if flag_3pk_conver else '2P'
    texts = ['RDP ($\sim${0:g} {1})'.format(bin_w_r, coord),
            '$d_{{field}}$ = {:.1E} $st/{}^{{2}}$'.format(field_dens, coord),
            '{} King profile'.format(kp_text),
            'r$_c$ = {0:g} $\pm$ {1:g} {2}'.format(rads_r[0], e_rads_r[0],
                coord),
            'r$_t$ = {0:g} $\pm$ {1:g} {2}'.format(rads_r[1], e_rads_r[1],
                coord),
            'r$_{{cl}}$ = {0:g} $\pm$ {1:g} {2}'.format(rads_r[2], e_rads_r[2],
                coord)]
    # Plot density profile with the smallest bin size
    ax5.plot(radii, ring_density, 'ko-', zorder=3, label=texts[0])
    # Plot poisson error bars
    plt.errorbar(radii, ring_density, yerr=poisson_error, fmt='ko',
                 zorder=1)
    # Plot background level.
    ax5.hlines(y=field_dens, xmin=0, xmax=max(radii),
               label=texts[1], color='b', zorder=5)
    # Approx middle of the graph.
    arr_y_up = (y_max - y_min) / 2.3 + y_min
    # Length and width of arrow head.
    head_w, head_l = x_max * 0.023, (y_max - y_min) * 0.045
    # Length of arrow.
    arr_y_dwn = -1. * abs(arr_y_up - field_dens) * 0.76
    # Plot 3-P King profile.
    if flag_3pk_conver is True:
        # Plot curve.
        ax5.plot(radii, three_params(radii, rt, cd, rc, field_dens), 'g--',
            label=texts[2], lw=2., zorder=3)
        # Plot r_t radius as an arrow. vline is there to show the label.
        ax5.vlines(x=rt, ymin=0., ymax=0., label=texts[4], color='g')
        ax5.arrow(rt, arr_y_up, 0., arr_y_dwn, fc="g", ec="g",
                  head_width=head_w, head_length=head_l, zorder=5)
        # Plot r_c as a dashed line.
        ax5.vlines(x=rc, ymin=0, ymax=three_params(rc, rt, cd, rc, field_dens),
            label=texts[3], color='g', linestyles=':', lw=4., zorder=4)
    # Plot 2-P King profile if 3-P was not found.
    elif flag_2pk_conver is True:
        # Plot curve.
        ax5.plot(radii, two_params(radii, cd, rc, field_dens), 'g--',
            label=texts[2], lw=2., zorder=3)
        # Plot r_c as a dashed line.
        ax5.vlines(x=rc, ymin=0, ymax=two_params(rc, cd, rc, field_dens),
                   label=texts[3], color='g', linestyles=':', lw=4., zorder=4)
    # Plot radius.
    ax5.vlines(x=clust_rad, ymin=0, ymax=0., label=texts[5], color='r')
    ax5.arrow(clust_rad, arr_y_up, 0., arr_y_dwn, fc="r",
              ec="r", head_width=head_w, head_length=head_l, zorder=5)
    # Plot radius error zone.
    if e_rad > 0.:
        plt.axvspan((clust_rad - e_rad), (clust_rad + e_rad), facecolor='grey',
            alpha=0.5)
    # get handles
    handles, labels = ax5.get_legend_handles_labels()
    # use them in the legend
    ax5.legend(handles, labels, loc='upper right', numpoints=2, fontsize=12)
    ax5.minorticks_on()

    # Zoom on x,y finding chart
    ax6 = plt.subplot(gs1[2:4, 6:8])
    #Set plot limits
    x_min, x_max = min(x_data), max(x_data)
    y_min, y_max = min(y_data), max(y_data)
    # If possible, zoom in.
    x_zmin, x_zmax = max(x_min, (center_cl[0] - 1.5 * clust_rad)), \
    min(x_max, (center_cl[0] + 1.5 * clust_rad))
    y_zmin, y_zmax = max(y_min, (center_cl[1] - 1.5 * clust_rad)), \
    min(y_max, (center_cl[1] + 1.5 * clust_rad))
    # Prevent axis stretching.
    if (x_zmax - x_zmin) != (y_zmax - y_zmin):
        lst = [(x_zmax - x_zmin), (y_zmax - y_zmin)]
        val, idx = min((val, idx) for (idx, val) in enumerate(lst))
        if idx == 0:
            x_zmax = x_zmin + lst[1]
        else:
            y_zmax = y_zmin + lst[0]
    plt.xlim(x_zmin, x_zmax)
    plt.ylim(y_zmin, y_zmax)
    #Set axis labels
    plt.xlabel('{} ({})'.format(x_name, coord), fontsize=12)
    plt.ylabel('{} ({})'.format(y_name, coord), fontsize=12)
    # Set minor ticks
    ax6.minorticks_on()
    # Add circle
    circle = plt.Circle((center_cl[0], center_cl[1]), clust_rad, color='r',
        fill=False, lw=1.5)
    fig.gca().add_artist(circle)
    text1 = 'Cluster (zoom)\n'
    text2 = 'CI = %0.2f' % (cont_index)
    text = text1 + text2
    plt.text(0.62, 0.9, text, transform=ax6.transAxes,
             bbox=dict(facecolor='white', alpha=0.85), fontsize=12)
    # Plot contour levels.
    # 'semi' and 'manual' modes do not produce any of these parameters.
    # Skip KDE plot if one of those modes was used.
    if kde_pl:
        ext_range, x, y, k_pos = kde_pl
        kde = np.reshape(k_pos.T, x.shape)
        plt.imshow(np.rot90(kde), cmap=plt.cm.YlOrBr, extent=ext_range)
        plt.contour(x, y, kde, 10, colors='k', linewidths=0.6)
    # Plot stars.
    plt.scatter(x_data, y_data, marker='o', c='black', s=st_sizes_arr, zorder=4)
    #Plot center.
    plt.scatter(center_cl[0], center_cl[1], color='w', s=40, lw=0.8,
        marker='x', zorder=5)

    # Cluster and field regions defined.
    ax7 = plt.subplot(gs1[4:6, 0:2])
    # Get max and min values in x,y
    x_min, x_max = min(x_data), max(x_data)
    y_min, y_max = min(y_data), max(y_data)
    #Set plot limits
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    #Set axis labels
    plt.xlabel('{} ({})'.format(x_name, coord), fontsize=12)
    plt.ylabel('{} ({})'.format(y_name, coord), fontsize=12)
    # Set minor ticks
    ax7.minorticks_on()
    ax7.grid(b=True, which='both', color='gray', linestyle='--', lw=0.5)
    # Radius
    circle = plt.Circle((center_cl[0], center_cl[1]), clust_rad,
                        color='k', fill=False)
    fig.gca().add_artist(circle)
    plt.text(0.4, 0.92, 'Cluster + %d Field regions' % (len(field_regions)),
             transform=ax7.transAxes,
             bbox=dict(facecolor='white', alpha=0.8), fontsize=12)
    # Plot cluster region.
    plt.scatter(zip(*cl_region)[1], zip(*cl_region)[2], marker='o', c='red',
                s=8, edgecolors='none')
    if not flag_area_stronger:
        # Plot field stars regions.
        col = cycle(['DimGray', 'ForestGreen', 'maroon', 'RoyalBlue'])
        for i, reg in enumerate(field_regions):
            stars_reg_temp = [[], []]
            for star in reg:
                # star[1] is the x coordinate and star[2] the y coordinate
                stars_reg_temp[0].append(star[1])
                stars_reg_temp[1].append(star[2])
            plt.scatter(stars_reg_temp[0], stars_reg_temp[1], marker='o',
                        c=next(col), s=8, edgecolors='none')

    # Field stars CMD (stars outside cluster's radius)
    ax8 = plt.subplot(gs1[4:6, 2:4])
    #Set plot limits
    plt.xlim(x_min_cmd, x_max_cmd)
    plt.ylim(y_min_cmd, y_max_cmd)
    #Set axis labels
    plt.xlabel('$' + x_ax + '$', fontsize=18)
    plt.ylabel('$' + y_ax + '$', fontsize=18)
    # Set minor ticks
    ax8.minorticks_on()
    # Only draw units on axis (ie: 1, 2, 3)
    ax8.xaxis.set_major_locator(MultipleLocator(1.0))
    # Set grid
    ax8.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
    # Plot rejected stars.
    # Outside of the cluster region.
    stars_rjct_temp = [[], []]
    for star in stars_out_rjct:
        stars_rjct_temp[0].append(star[5])
        stars_rjct_temp[1].append(star[3])
    plt.scatter(stars_rjct_temp[0], stars_rjct_temp[1], marker='x', c='teal',
                s=15, zorder=1)
    # Plot stars within the field regions defined.
    stars_acpt_temp = [[], []]
    for fr in field_regions:
        for star in fr:
            stars_acpt_temp[0].append(star[5])
            stars_acpt_temp[1].append(star[3])
    sz_pt = 0.2 if (len(stars_out_rjct) + len(stars_out)) > 5000 else 0.5
    plt.scatter(stars_acpt_temp[0], stars_acpt_temp[1], marker='o', c='k',
                s=sz_pt, zorder=2)
    plt.text(0.53, 0.93, '$r > r_{cl}\,|\,N=%d$' % len(stars_acpt_temp[0]),
        transform=ax8.transAxes, bbox=dict(facecolor='white', alpha=0.5),
        fontsize=16)

    # Cluster's stars CMD (stars inside cluster's radius)
    ax9 = plt.subplot(gs1[4:6, 4:6])
    #Set plot limits
    plt.xlim(x_min_cmd, x_max_cmd)
    plt.ylim(y_min_cmd, y_max_cmd)
    #Set axis labels
    plt.xlabel('$' + x_ax + '$', fontsize=18)
    plt.ylabel('$' + y_ax + '$', fontsize=18)
    # Set minor ticks
    ax9.minorticks_on()
    # Only draw units on axis (ie: 1, 2, 3)
    ax9.xaxis.set_major_locator(MultipleLocator(1.0))
    # Set grid
    ax9.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
    # Calculate total number of stars whitin cluster's radius.
    tot_stars = len(stars_in_rjct) + len(cl_region)
    text1 = '$r \leq r_{{cl}}\,|\,N={}$'.format(tot_stars)
    text2 = r'$n_{{memb}} \approx {}$'.format(n_memb)
    text = text1 + '\n' + text2
    plt.text(0.55, 0.87, text, transform=ax9.transAxes,
             bbox=dict(facecolor='white', alpha=0.5), fontsize=16)
    # Plot stars in CMD.
    if len(stars_in_rjct) > 0:
        # Only attempt to pot if any star is stored in the list.
        plt.scatter(zip(*stars_in_rjct)[5], zip(*stars_in_rjct)[3], marker='x',
            c='teal', s=12, zorder=1)
    sz_pt = 0.5 if (len(stars_in_rjct) + len(cl_region)) > 1000 else 1.
    plt.scatter(zip(*cl_region)[5], zip(*cl_region)[3], marker='o', c='k',
                s=sz_pt, zorder=2)

    # Magnitude error
    ax10 = plt.subplot(gs1[4, 6:8])
    #Set plot limits
    x_min, x_max = min(mag_data) - 0.5, max(mag_data) + 0.5
    plt.xlim(x_min, x_max)
    plt.ylim(-0.005, e_max + (e_max / 5.))
    #Set axis labels
    plt.ylabel('$\sigma_{' + y_ax + '}$', fontsize=18)
    plt.xlabel('$' + y_ax + '$', fontsize=18)
    # Set minor ticks
    ax10.minorticks_on()
    # Plot e_max line.
    ax10.hlines(y=e_max, xmin=x_min, xmax=x_max, color='r',
        linestyles='dashed', zorder=2)
    # Plot rectangle.
    bright_end = min(mag_data) + be
    ax10.vlines(x=bright_end + 0.05, ymin=-0.005, ymax=be_e, color='r',
               linestyles='dashed', zorder=2)
    ax10.vlines(x=min(mag_data) - 0.05, ymin=-0.005, ymax=be_e, color='r',
               linestyles='dashed', zorder=2)
    ax10.hlines(y=be_e, xmin=min(mag_data), xmax=bright_end, color='r',
               linestyles='dashed', zorder=2)
    # If any method could be used.
    if err_all_fallback is False and err_max_fallback is False:
        # Plot curve(s) according to the method used.
        if er_mode == 'eyefit':
            # Unpack params.
            popt_umag, pol_mag, popt_ucol1, pol_col1, mag_val_left, \
            mag_val_right, col1_val_left, col1_val_right = err_plot

            # Plot left side of upper envelope (exponential).
            ax10.plot(mag_val_left, exp_func(mag_val_left, *popt_umag), 'r--',
                lw=2., zorder=3)
            # Plot right side of upper envelope (polynomial).
            ax10.plot(mag_val_right, np.polyval(pol_mag, (mag_val_right)),
                'r--', lw=2., zorder=3)
        elif er_mode == 'lowexp':
            # Unpack params.
            popt_mag, popt_col1 = err_plot
            # Plot exponential curve.
            mag_x = np.linspace(bright_end, max(mag_data), 50)
            ax10.plot(mag_x, exp_func(mag_x, *popt_mag), 'r-', zorder=3)
    # Plot rejected stars.
    if len(stars_out_rjct) > 0:
        # Only attempt to pot if any star is stored in the list.
        plt.scatter(zip(*stars_out_rjct)[3], zip(*stars_out_rjct)[4],
            marker='x', c='teal', s=15, zorder=1)
    if len(stars_in_rjct) > 0:
        plt.scatter(zip(*stars_in_rjct)[3], zip(*stars_in_rjct)[3], marker='x',
            c='teal', s=15, zorder=1)
    # Plot accepted stars.
    plt.scatter(zip(*cl_region)[3], zip(*cl_region)[4], marker='o', c='k',
                s=1, zorder=2)
    plt.scatter(zip(*stars_out)[3], zip(*stars_out)[4], marker='o', c='k',
                s=1, zorder=2)

    # Color error
    ax11 = plt.subplot(gs1[5, 6:8])
    #Set plot limits
    x_min, x_max = min(mag_data) - 0.5, max(mag_data) + 0.5
    plt.xlim(x_min, x_max)
    plt.ylim(-0.005, e_max + (e_max / 5.))
    #Set axis labels
    plt.ylabel('$\sigma_{' + x_ax + '}$', fontsize=18)
    plt.xlabel('$' + y_ax + '$', fontsize=18)
    # Set minor ticks
    ax11.minorticks_on()
    # Plot e_max line.
    ax11.hlines(y=e_max, xmin=x_min, xmax=x_max, color='r',
               linestyles='dashed', zorder=2)
    # Plot rectangle.
    bright_end = min(mag_data) + be
    ax11.vlines(x=bright_end + 0.05, ymin=-0.005, ymax=be_e, color='r',
               linestyles='dashed', zorder=2)
    ax11.vlines(x=min(mag_data) - 0.05, ymin=-0.005, ymax=be_e, color='r',
               linestyles='dashed', zorder=2)
    ax11.hlines(y=be_e, xmin=min(mag_data), xmax=bright_end, color='r',
               linestyles='dashed', zorder=2)
    # If any method could be used.
    if err_all_fallback is False and err_max_fallback is False:
        # Plot curve(s) according to the method used.
        if er_mode == 'eyefit':
            # Unpack params.
            popt_mag, pol_mag, popt_col1, pol_col1, mag_val_left, \
            mag_val_right, col1_val_left, col1_val_right = err_plot

            # Plot left side: exponential envelope.
            ax11.plot(col1_val_left, exp_func(col1_val_left, *popt_col1), 'r--',
                lw=2., zorder=3)
            # Plot right side: polynomial envelope.
            ax11.plot(col1_val_right, np.polyval(pol_col1, (col1_val_right)),
                'r--', lw=2., zorder=3)
        elif er_mode == 'lowexp':
            # Unpack params.
            popt_mag, popt_col1 = err_plot
            # Plot exponential curve.
            ax11.plot(mag_x, exp_func(mag_x, *popt_col1), 'r-', zorder=3)
    # Plot rejected stars.
    if len(stars_out_rjct) > 0:
        # Only attempt to pot if any star is stored in the list.
        plt.scatter(zip(*stars_out_rjct)[3], zip(*stars_out_rjct)[6],
            marker='x', c='teal', s=15, zorder=1)
    if len(stars_in_rjct) > 0:
        plt.scatter(zip(*stars_in_rjct)[3], zip(*stars_in_rjct)[6], marker='x',
            c='teal', s=15, zorder=1)
    # Plot accepted stars.
    plt.scatter(zip(*cl_region)[3], zip(*cl_region)[6], marker='o', c='k',
                s=1, zorder=2)
    plt.scatter(zip(*stars_out)[3], zip(*stars_out)[6], marker='o', c='k',
                s=1, zorder=2)

    # LF of stars in cluster region and outside.
    ax12 = plt.subplot(gs1[6:8, 0:2])
    #Set plot limits
    x_min, x_max = min(mag_data) - 0.5, max(mag_data) + 0.5
    plt.xlim(x_max, x_min)
    ax12.minorticks_on()
    # Only draw units on axis (ie: 1, 2, 3)
    ax12.xaxis.set_major_locator(MultipleLocator(2.0))
    # Set grid
    ax12.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
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
    plt.ylim(0., plt.ylim()[1])
    # Completeness maximum value.
    # completeness = [max_mag, bin_edges, max_indx, comp_perc]
    bin_edges, max_indx = completeness[1], completeness[2]
    mag_peak = bin_edges[max_indx]
    text = '$' + y_ax + r',_{compl}\,\approx\,%0.1f$' % mag_peak
    ax12.vlines(x=mag_peak, ymin=0., ymax=plt.ylim()[1], color='k',
        lw=1.5, linestyles='dashed', label=text, zorder=1)
    # Legends.
    leg11 = plt.legend(fancybox=True, loc='upper right', numpoints=1,
                       fontsize=13)
    # Set the alpha value of the legend.
    leg11.get_frame().set_alpha(0.7)

    # Integrated magnitude and color.
    if integr_return:
        # Unpack values.
        cl_reg_mag1, fl_reg_mag1, integ_mag1, cl_reg_mag2, fl_reg_mag2, \
        integ_mag2 = integr_return
        # Make plot
        ax13 = plt.subplot(gs1[6:8, 2:4])
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
            x_min, x_max = min(min(cl_reg_mag1[0]), min(cl_reg_mag2[0])) - 0.2,\
            max(max(cl_reg_mag1[0]), max(cl_reg_mag2[0])) + 0.2
            y_min, y_max = max(max(cl_reg_mag1[1]), max(cl_reg_mag2[1])) + 0.2,\
            min(min(cl_reg_mag1[1]), min(cl_reg_mag2[1])) - 0.2
        plt.xlim(x_min, x_max)
        plt.ylim(y_min, y_max)
        ax13.set_xlabel('$mag$', fontsize=18)
        ax13.set_ylabel('$mag^*$', fontsize=18)
        ax13.minorticks_on()
        ax13.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
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
            plt.plot(fl_reg_mag1[0], fl_reg_mag1[1], 'b-', lw=1., label=text3)
            # Field average integrated magnitude.
            plt.plot(fl_reg_mag2[0], fl_reg_mag2[1], 'b:', lw=2., label=text4)
        # Check how the second magnitude whould be formed.
        if m_ord == 21:
            sig, text0 = 1., x_ax0 + '^{*} -' + y_ax
        elif m_ord == 12:
            sig, text0 = -1., y_ax + '^{*} -' + x_ax0
        int_col = sig * (integ_mag2 - integ_mag1)
        text = '$(' + text0 + '^{*} )_{cl} = %0.2f$' % int_col
        plt.text(0.25, 0.15, text, transform=ax13.transAxes,
             bbox=dict(facecolor='white', alpha=0.75), fontsize=13)
        lines, labels = ax13.get_legend_handles_labels()
        leg = ax13.legend(lines, labels, loc='lower right', numpoints=1,
            fontsize=13)
        leg.get_frame().set_alpha(0.75)

    # Distribution of p_values.
    if flag_pval_test and not flag_area_stronger:
        # Extract parameters from list.
        prob_cl_kde, p_vals_cl, p_vals_f, kde_cl_1d, kde_f_1d, x_kde, y_over\
        = pval_test_params
        ax14 = plt.subplot(gs1[6:8, 4:6])
        plt.xlim(-0.5, 1.5)
        plt.ylim(0, max(max(kde_f_1d), max(kde_cl_1d)) + 0.5)
        plt.xlabel('p-values', fontsize=12)
        plt.ylabel('Density', fontsize=12)
        ax14.minorticks_on()
        ax14.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
        # Grid to background.
        ax14.set_axisbelow(True)
        # Plot cluster vs field KDE.
        plt.plot(x_kde, kde_cl_1d, c='b', ls='-', lw=1., label='$KDE_{cl}$')
        # Plot field vs field KDE.
        plt.plot(x_kde, kde_f_1d, c='r', ls='-', lw=1., label='$KDE_{f}$')
        # Fill overlap.
        plt.fill_between(x_kde, y_over, 0, color='grey', alpha='0.5')
        text = '$P_{cl}^{KDE} = %0.2f$' % round(prob_cl_kde, 2)
        plt.text(0.05, 0.92, text, transform=ax14.transAxes,
             bbox=dict(facecolor='white', alpha=0.6), fontsize=12)
        # Legend.
        handles, labels = ax14.get_legend_handles_labels()
        leg = ax14.legend(handles, labels, loc='upper right', numpoints=1,
                          fontsize=12)
        leg.get_frame().set_alpha(0.6)

    plot_colorbar = False
    if da_params[0] != 'skip':
        # Histogram for the distribution of membership probabilities from the
        # decontamination algorithm.
        ax16 = plt.subplot(gs1[8:10, 0:2])
        plt.xlim(0., 1.)
        plt.xlabel('membership probability', fontsize=12)
        plt.ylabel('N', fontsize=12)
        ax16.minorticks_on()
        ax16.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
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
        if bf_flag:
            # Plot minimum probability value used in best isochrone fit
            # function.
            text = '$prob_{min}=%.2f$' % min_prob
            plt.text(0.05, 0.92, text, transform=ax16.transAxes,
                 bbox=dict(facecolor='white', alpha=0.85), fontsize=14)
            plt.axvline(x=min_prob, linestyle='--', color='green', lw=2.5)

        # Finding chart of cluster region with decontamination algorithm
        # applied and colors assigned according to the probabilities obtained.
        ax17 = plt.subplot(gs1[8:10, 2:4])
        #Set plot limits, Use 'zoom' x,y ranges.
        plt.xlim(x_zmin, x_zmax)
        plt.ylim(y_zmin, y_zmax)
        #Set axis labels
        plt.xlabel('x (px)', fontsize=12)
        plt.ylabel('y (px)', fontsize=12)
        # Set minor ticks
        ax17.minorticks_on()
        # Radius
        circle = plt.Circle((center_cl[0], center_cl[1]), clust_rad,
                            color='red', fill=False)
        fig.gca().add_artist(circle)
        plt.text(0.63, 0.93, 'Cluster region', transform=ax17.transAxes,
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

    if da_params[0] != 'skip' or bf_flag:
        # Star's membership probabilities on cluster's CMD.
        ax18 = plt.subplot(gs1[8:10, 4:6])
        #Set plot limits
        plt.xlim(x_min_cmd, x_max_cmd)
        plt.ylim(y_min_cmd, y_max_cmd)
        #Set axis labels
        plt.xlabel('$' + x_ax + '$', fontsize=18)
        plt.ylabel('$' + y_ax + '$', fontsize=18)
        tot_clust = len(memb_prob_avrg_sort)
        text = '$r \leq r_{cl}\,|\,N=%d$' % tot_clust
        plt.text(0.5, 0.93, text, transform=ax18.transAxes,
                 bbox=dict(facecolor='white', alpha=0.5), fontsize=16)
        # Set minor ticks
        ax18.minorticks_on()
        ax18.xaxis.set_major_locator(MultipleLocator(1.0))
        ax18.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
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
        if bf_flag:
            # Plot isochrone if best fit process was used.
            plt.plot(shift_isoch[0], shift_isoch[1], c=c_iso, lw=1.2)
        sca = plt.scatter(m_p_m_temp_inv[0], m_p_m_temp_inv[1], marker='o',
                    c=col_select, s=40, cmap=cm, lw=0.5, vmin=v_min_mp,
                    vmax=v_max_mp)
        # If list is not empty.
        if m_p_m_temp_inv[1]:
            # Plot error bars at several mag values.
            mag_y = np.arange(int(min(m_p_m_temp_inv[1]) + 0.5),
                              int(max(m_p_m_temp_inv[1]) + 0.5) + 0.1)
            x_val = [min(x_max_cmd, max(col1_data) + 0.2) - 0.4] * len(mag_y)
            # Read average fitted values for exponential error fit.
            popt_mag, popt_col1 = err_lst[:2]
            plt.errorbar(x_val, mag_y, yerr=exp_func(mag_y, *popt_mag),
                         xerr=exp_func(mag_y, *popt_col1), fmt='k.', lw=0.8,
                         ms=0., zorder=4)
            # Plot colorbar (see bottom of file).
            if v_min_mp != v_max_mp:
                plot_colorbar = True

    # Synthetic cluster.
    if bf_flag:
        ax19 = plt.subplot(gs1[8:10, 6:8])
        #Set plot limits
        plt.xlim(x_min_cmd, x_max_cmd)
        plt.ylim(y_min_cmd, y_max_cmd)
        #Set axis labels
        plt.xlabel('$' + x_ax + '$', fontsize=18)
        plt.ylabel('$' + y_ax + '$', fontsize=18)
        # Set minor ticks
        ax19.minorticks_on()
        ax19.xaxis.set_major_locator(MultipleLocator(1.0))
        ax19.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
        # Add text box
        # See if bootstrap process was applied.
        if N_b < 2:
            # Round cluster params using the g format.
            cp_r = ['{:g}'.format(_) for _ in isoch_fit_params[0]]
            cp_e = isoch_fit_errors
        else:
            # Round cluster parameters.
            cp_r, cp_e = err_r.round_sig_fig(isoch_fit_params[0],
                isoch_fit_errors)
        text1 = '$N = {}$\n'.format(len(synth_clst[0]))
        text2 = '$z = {} \pm {}$\n'.format(cp_r[0], cp_e[0])
        text3 = '$log(age) = {} \pm {}$\n'.format(cp_r[1], cp_e[1])
        text4 = '$E_{{(B-V)}} = {} \pm {}$\n'.format(cp_r[2], cp_e[2])
        text5 = '$(m-M)_o = {} \pm {}$\n'.format(cp_r[3], cp_e[3])
        text6 = '$M_{{\odot}} = {} \pm {}$\n'.format(cp_r[4], cp_e[4])
        text7 = '$b_{{frac}} = {} \pm {}$'.format(cp_r[5], cp_e[5])
        text = text1 + text2 + text3 + text4 + text5 + text6 + text7
        plt.text(0.5, 0.65, text, transform=ax19.transAxes,
                 bbox=dict(facecolor='white', alpha=0.6), fontsize=12)
        # Plot isochrone.
        plt.plot(shift_isoch[0], shift_isoch[1], 'r', lw=1.2)
        # Plot synth clust.
        plt.scatter(synth_clst[0], synth_clst[2], marker='o', s=40,
                    c='#4682b4', lw=0.5)

    # Best fitting process plots for GA.
    if bf_flag and best_fit_algor == 'genet':

        m, a, e, d, mass, binar_f = isoch_fit_params[0]
        if N_b >= 2:
            e_m, e_a, e_e, e_d, e_mass, e_bin = isoch_fit_errors
        else:
            e_m, e_a, e_e, e_d, e_mass, e_bin = [0.] * len(isoch_fit_errors)

        # Set ranges used by plots below.
        min_max_p = []
        for param in ip_list[1]:
            if max(param) != min(param):
                delta_p = (max(param) - min(param)) * 0.05
            else:
                # The first max is there for when max(param)=0 (for example
                # when 0 reddening is used)
                delta_p = max(max(param) / 100., 0.0005)
            min_max_p.append([min(param) - delta_p, max(param) + delta_p])
        m_min, m_max = min_max_p[0]
        a_min, a_max = min_max_p[1]
        e_min, e_max = min_max_p[2]
        d_min, d_max = min_max_p[3]
        mass_min, mass_max = min_max_p[4]
        bin_min, bin_max = min_max_p[5]

        lkl_old, new_bs_indx, model_done = isoch_fit_params[1], \
        isoch_fit_params[2], isoch_fit_params[3]

        # Age vs metallicity GA diagram.
        plt.subplot(gs1[10:12, 0:2])
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
            ax20 = plt.gca()
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

        # GA diagram.
        ax21 = plt.subplot(gs1[10:12, 2:6])
        plt.xlim(-0.5, n_gen + int(0.01 * n_gen))
        lkl_range = max(lkl_old[1]) - min(lkl_old[0])
        plt.ylim(min(lkl_old[0]) - 0.1 * lkl_range,
                 max(lkl_old[1]) + 0.1 * lkl_range)
        # Set minor ticks
        ax21.minorticks_on()
        ax21.tick_params(axis='y', which='major', labelsize=9)
        ax21.grid(b=True, which='major', color='gray', linestyle='--', lw=0.6)
        plt.xlabel('Generation', fontsize=12)
        plt.ylabel('Likelihood', fontsize=12)
        text1 = '$N_{total} = %.2e\,;\,N_{btst} = %d$' '\n' % \
        (len(model_done[0]), N_b)
        text2 = '$n_{gen}=%d\,;\,n_{pop}=%d$' '\n' % (n_gen, n_pop)
        text3 = '$f_{dif}=%0.2f\,;\,cr_{sel}=%s$' '\n' % (fdif, cr_sel)
        text4 = '$p_{cross}=%0.2f\,;\,p_{mut}=%0.2f$' '\n' % (p_cross, p_mut)
        text5 = '$n_{el}=%d\,;\,n_{ei}=%d\,;\,n_{es}=%d$' % (n_el, n_ei, n_es)
        text = text1 + text2 + text3 + text4 + text5
        plt.text(0.05, 0.73, text, transform=ax21.transAxes,
            bbox=dict(facecolor='white', alpha=0.75), fontsize=12)
        # Plot likelihood minimum and mean lines.
        ax21.plot(range(len(lkl_old[0])), lkl_old[0], lw=1., c='black',
                  label='$L_{min}$')
        ax21.plot(range(len(lkl_old[0])), lkl_old[1], lw=1., c='blue',
                  label='$L_{mean}$')
        # Plot line marking a new best solution found.
        for lin in new_bs_indx:
            lw_lin = 2. if lin > 0.05 * len(lkl_old[0]) else 0.5
            plt.axvline(x=lin, linestyle='--', lw=lw_lin, color='green')
        # Legend.
        handles, labels = ax21.get_legend_handles_labels()
        leg = ax21.legend(handles, labels, loc='upper right', numpoints=1,
                          fontsize=12)
        leg.get_frame().set_alpha(0.6)

        # Extinction vs distance modulus GA diagram.
        plt.subplot(gs1[10:12, 6:8])
        plt.xlim(e_min, e_max)
        plt.ylim(d_min, d_max)
        plt.xlabel('$E_{(B-V)}$', fontsize=16)
        plt.ylabel('$(m-M)_o$', fontsize=16)
        plt.minorticks_on()
        # Plot best fit point.
        plt.scatter(e, d, marker='o', c='b', s=30)
        # Check if errors in both dimensions are defined.
        if all([i > 0. for i in [e_e, e_d]]):
            # Plot ellipse error.
            ax21 = plt.gca()
            ellipse = Ellipse(xy=(e, d), width=2 * e_e, height=2 * e_d,
                                    edgecolor='b', fc='None', lw=1.)
            ax21.add_patch(ellipse)
        elif e_e < 0.:
            plt.errorbar(e, d, yerr=e_d, color='b')
        elif e_d < 0.:
            plt.errorbar(e, d, xerr=e_e, color='b')
        # Plot density map.
        hist, xedges, yedges = np.histogram2d(zip(*model_done[0])[2],
                                              zip(*model_done[0])[3], bins=100)
        # H_g is the 2D histogram with a gaussian filter applied
        h_g = gaussian_filter(hist, 2, mode='constant')
        plt.imshow(h_g.transpose(), origin='lower',
                   extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
                   cmap=plt.get_cmap('Reds'), aspect='auto')

        #
        # Metallicity/likelihood plot.
        ax22 = plt.subplot(gs1[12:14, 0:2])
        plt.ylim(max(0, min(lkl_old[0]) - 0.3 * min(lkl_old[0])),
            max(lkl_old[1]))
        plt.xlim(m_min, m_max)
        # Set minor ticks
        ax22.minorticks_on()
        ax22.tick_params(axis='y', which='major', labelsize=9)
        plt.ylabel('Likelihood', fontsize=12)
        plt.xlabel('$z$', fontsize=16)
        text = '$z = {} \pm {}$'.format(cp_r[0], cp_e[0])
        plt.text(0.1, 0.93, text, transform=ax22.transAxes,
            bbox=dict(facecolor='white', alpha=0.5), fontsize=12)
        hist, xedges, yedges = np.histogram2d(zip(*model_done[0])[0],
                                              model_done[1], bins=100)
        # H_g is the 2D histogram with a gaussian filter applied
        h_g = gaussian_filter(hist, 2, mode='constant')
        y_min_edge = max(0, min(lkl_old[0]) - 0.3 * min(lkl_old[0]))
        plt.imshow(h_g.transpose(), origin='lower',
                   extent=[xedges[0], xedges[-1], y_min_edge, yedges[-1]],
                   cmap=plt.get_cmap('gist_yarg'), aspect='auto')
        plt.axvline(x=m + e_m, linestyle='--', color='red')
        plt.axvline(x=m - e_m, linestyle='--', color='red')
        plt.axvline(x=m, linestyle='--', color='blue', zorder=3)

        #
        # Age/likelihood plot.
        ax23 = plt.subplot(gs1[12:14, 2:4])
        plt.ylim(max(0, min(lkl_old[0]) - 0.3 * min(lkl_old[0])),
            max(lkl_old[1]))
        plt.xlim(a_min, a_max)
        # Set minor ticks
        ax23.minorticks_on()
        ax23.tick_params(axis='y', which='major', labelsize=9)
        plt.ylabel('Likelihood', fontsize=12)
        plt.xlabel('$log(age)$', fontsize=16)
        text = '$log(age) = {} \pm {}$'.format(cp_r[1], cp_e[1])
        plt.text(0.1, 0.93, text, transform=ax23.transAxes,
            bbox=dict(facecolor='white', alpha=0.5), fontsize=12)
        hist, xedges, yedges = np.histogram2d(zip(*model_done[0])[1],
                                              model_done[1], bins=100)
        # H_g is the 2D histogram with a gaussian filter applied
        h_g = gaussian_filter(hist, 2, mode='constant')
        plt.imshow(h_g.transpose(), origin='lower',
                   extent=[xedges[0], xedges[-1], y_min_edge, yedges[-1]],
                   cmap=plt.get_cmap('gist_yarg'), aspect='auto')
        plt.axvline(x=a + e_a, linestyle='--', color='red')
        plt.axvline(x=a - e_a, linestyle='--', color='red')
        plt.axvline(x=a, linestyle='--', color='blue', zorder=3)

        #
        # Reddening/likelihood plot.
        ax24 = plt.subplot(gs1[12:14, 4:6])
        plt.ylim(max(0, min(lkl_old[0]) - 0.3 * min(lkl_old[0])),
            max(lkl_old[1]))
        plt.xlim(e_min, e_max)
        # Set minor ticks
        ax24.minorticks_on()
        ax24.tick_params(axis='y', which='major', labelsize=9)
        plt.ylabel('Likelihood', fontsize=12)
        plt.xlabel('$E_{(B-V)}$', fontsize=16)
        text = '$E_{{(B-V)}} = {} \pm {}$'.format(cp_r[2], cp_e[2])
        plt.text(0.1, 0.93, text, transform=ax24.transAxes,
            bbox=dict(facecolor='white', alpha=0.5), fontsize=12)
        hist, xedges, yedges = np.histogram2d(zip(*model_done[0])[2],
                                              model_done[1], bins=100)
        # H_g is the 2D histogram with a gaussian filter applied
        h_g = gaussian_filter(hist, 2, mode='constant')
        plt.imshow(h_g.transpose(), origin='lower',
                   extent=[xedges[0], xedges[-1], y_min_edge, yedges[-1]],
                   cmap=plt.get_cmap('gist_yarg'), aspect='auto')
        plt.axvline(x=e + e_e, linestyle='--', color='red')
        plt.axvline(x=e - e_e, linestyle='--', color='red')
        plt.axvline(x=e, linestyle='--', color='blue', zorder=3)

        #
        # Dist/likelihood plot.
        ax25 = plt.subplot(gs1[12:14, 6:8])
        plt.ylim(max(0, min(lkl_old[0]) - 0.3 * min(lkl_old[0])),
            max(lkl_old[1]))
        plt.xlim(d_min, d_max)
        # Set minor ticks
        ax25.minorticks_on()
        ax25.tick_params(axis='y', which='major', labelsize=9)
        plt.ylabel('Likelihood', fontsize=12)
        plt.xlabel('$(m-M)_o$', fontsize=16)
        text = '$(m-M)_o = {} \pm {}$'.format(cp_r[3], cp_e[3])
        plt.text(0.1, 0.93, text, transform=ax25.transAxes,
            bbox=dict(facecolor='white', alpha=0.5), fontsize=12)
        hist, xedges, yedges = np.histogram2d(zip(*model_done[0])[3],
                                              model_done[1], bins=100)
        # H_g is the 2D histogram with a gaussian filter applied
        h_g = gaussian_filter(hist, 2, mode='constant')
        plt.imshow(h_g.transpose(), origin='lower',
                   extent=[xedges[0], xedges[-1], y_min_edge, yedges[-1]],
                   cmap=plt.get_cmap('gist_yarg'), aspect='auto')
        plt.axvline(x=d + e_d, linestyle='--', color='red')
        plt.axvline(x=d - e_d, linestyle='--', color='red')
        plt.axvline(x=d, linestyle='--', color='blue', zorder=3)

        #
        # Mass/binary plot.
        ax26 = plt.subplot(gs1[14:16, 0:2])
        plt.xlim(mass_min, mass_max)
        plt.ylim(bin_min, bin_max)
        plt.xlabel('$M_{\odot}$', fontsize=16)
        plt.ylabel('$b_{frac}$', fontsize=16)
        plt.minorticks_on()
        # Plot best fit point.
        plt.scatter(mass, binar_f, marker='o', c='b', s=30)
        # Check if errors in both dimensions are defined.
        if all([i > 0. for i in [e_mass, e_bin]]):
            # Plot ellipse error.
            ax26 = plt.gca()
            ellipse = Ellipse(xy=(mass, binar_f), width=2 * e_mass,
                height=2 * e_bin, edgecolor='b', fc='None', lw=1.)
            ax26.add_patch(ellipse)
        elif e_mass < 0.:
            plt.errorbar(mass, binar_f, yerr=e_bin, color='b')
        elif e_bin < 0.:
            plt.errorbar(mass, binar_f, xerr=e_mass, color='b')
        # Plot density map.
        hist, xedges, yedges = np.histogram2d(zip(*model_done[0])[4],
                                              zip(*model_done[0])[5], bins=100)
        # H_g is the 2D histogram with a gaussian filter applied
        h_g = gaussian_filter(hist, 2, mode='constant')
        plt.imshow(h_g.transpose(), origin='lower',
                   extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
                   cmap=plt.get_cmap('Reds'), aspect='auto')

        #
        # Mass/likelihood plot.
        ax27 = plt.subplot(gs1[14:16, 2:4])
        plt.ylim(max(0, min(lkl_old[0]) - 0.3 * min(lkl_old[0])),
            max(lkl_old[1]))
        plt.xlim(mass_min, mass_max)
        # Set minor ticks
        ax27.minorticks_on()
        ax27.tick_params(axis='y', which='major', labelsize=9)
        plt.ylabel('Likelihood', fontsize=12)
        plt.xlabel('$M_{\odot}$', fontsize=16)
        text = '$M_{{\odot}} = {} \pm {}$'.format(cp_r[4], cp_e[4])
        plt.text(0.1, 0.93, text, transform=ax27.transAxes,
            bbox=dict(facecolor='white', alpha=0.5), fontsize=12)
        hist, xedges, yedges = np.histogram2d(zip(*model_done[0])[4],
                                              model_done[1], bins=100)
        # H_g is the 2D histogram with a gaussian filter applied
        h_g = gaussian_filter(hist, 2, mode='constant')
        plt.imshow(h_g.transpose(), origin='lower',
                   extent=[xedges[0], xedges[-1], y_min_edge, yedges[-1]],
                   cmap=plt.get_cmap('gist_yarg'), aspect='auto')
        plt.axvline(x=mass + e_mass, linestyle='--', color='red')
        plt.axvline(x=mass - e_mass, linestyle='--', color='red')
        plt.axvline(x=mass, linestyle='--', color='blue', zorder=3)

        #
        # Binary/likelihood plot.
        ax28 = plt.subplot(gs1[14:16, 4:6])
        plt.ylim(max(0, min(lkl_old[0]) - 0.3 * min(lkl_old[0])),
            max(lkl_old[1]))
        plt.xlim(bin_min, bin_max)
        # Set minor ticks
        ax28.minorticks_on()
        ax28.tick_params(axis='y', which='major', labelsize=9)
        plt.ylabel('Likelihood', fontsize=12)
        plt.xlabel('$b_{frac}$', fontsize=16)
        text = '$b_{{frac}} = {} \pm {}$'.format(cp_r[5], cp_e[5])
        plt.text(0.1, 0.93, text, transform=ax28.transAxes,
            bbox=dict(facecolor='white', alpha=0.5), fontsize=12)
        hist, xedges, yedges = np.histogram2d(zip(*model_done[0])[5],
                                              model_done[1], bins=100)
        # H_g is the 2D histogram with a gaussian filter applied
        h_g = gaussian_filter(hist, 2, mode='constant')
        plt.imshow(h_g.transpose(), origin='lower',
                   extent=[xedges[0], xedges[-1], y_min_edge, yedges[-1]],
                   cmap=plt.get_cmap('gist_yarg'), aspect='auto')
        plt.axvline(x=binar_f + e_bin, linestyle='--', color='red')
        plt.axvline(x=binar_f - e_bin, linestyle='--', color='red')
        plt.axvline(x=binar_f, linestyle='--', color='blue', zorder=3)

    # Ignore warning issued by colorbar plotted in CMD with membership
    # probabilities.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fig.tight_layout()

    # Plot colorbar down here so tight_layout won't move it around.
    if plot_colorbar is True:
        import matplotlib
        # Position and dimensions relative to the axes.
        x0, y0, width, height = [0.6, 0.85, 0.2, 0.04]
        # Transform them to get the ABSOLUTE POSITION AND DIMENSIONS
        Bbox = matplotlib.transforms.Bbox.from_bounds(x0, y0, width, height)
        trans = ax18.transAxes + fig.transFigure.inverted()
        l, b, w, h = matplotlib.transforms.TransformedBbox(Bbox, trans).bounds
        # Create the axes and the colorbar.
        cbaxes = fig.add_axes([l, b, w, h])
        cbar = plt.colorbar(sca, cax=cbaxes, ticks=[v_min_mp, v_max_mp],
            orientation='horizontal')
        cbar.ax.tick_params(labelsize=9)

    # Generate output file for each data file.
    pl_fmt, pl_dpi = pl_params[1:3]
    plt.savefig(join(output_subdir, str(clust_name) + '.' + pl_fmt), dpi=pl_dpi)

    # Close to release memory.
    plt.clf()
    plt.close()
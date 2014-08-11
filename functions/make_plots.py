"""
@author: gabriel
"""

from functions.exp_function import exp_func
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle
from matplotlib.ticker import MultipleLocator
from itertools import cycle
from scipy.ndimage.filters import gaussian_filter
from matplotlib.patches import Ellipse
from os.path import join
import warnings


def make_plots(output_subdir, clust_name, x_data, y_data, center_params,
    rdp_params, field_dens, radius_params,
    cont_index, mag_data, col1_data, err_plot, err_flags, kp_params,
    stars_in, stars_out, stars_in_rjct, stars_out_rjct, integr_return, n_c,
    flag_area_stronger, cluster_region, field_region, flag_pval_test,
    pval_test_params, memb_prob_avrg_sort, lum_func, completeness, bf_params,
    red_return, err_lst, bf_return, ga_params, er_params, axes_params,
    ps_params, pl_params):
    '''
    Make all plots.
    '''

    def star_size(x, a, c, area):
        '''
        Function to obtain the optimal star size for the scatter plot.
        '''
        return sum(a * np.exp(x * mag_data ** c)) / area - 0.001

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

    # Name for axes.
    x_ax, y_ax = axes_params[0], axes_params[1]

    # Define plot limits for *all* CMD diagrams.
    x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd = axes_params[2]
    col1_min, col1_max = max(x_min_cmd, min(col1_data) - 0.2),\
    min(x_max_cmd, max(col1_data) + 0.2)
    mag_min, mag_max = min(y_max_cmd, max(mag_data) + 0.5),\
    max(y_min_cmd, min(mag_data) - 0.5)

    # Unpack params.
    # Selected system params.
    cmd_select = ps_params[1]
    m_rs, a_rs, e_rs, d_rs = ps_params[3:]
    # Parameters from get_center function.
    bin_list, h_filter, bin_center, centers_kde, cent_stats, kde_pl = \
    center_params[0], center_params[3], center_params[4], center_params[5], \
    center_params[6], center_params[7]
    center_cl = [center_params[5][0][0], center_params[5][0][1]]
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
    # Integrated magnitude distribution.
    cl_reg_mag, fl_reg_mag, integ_mag, cl_reg_col, fl_reg_col, integ_col =\
    integr_return
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
    fig = plt.figure(figsize=(20, 35))  # create the top-level container
    gs1 = gridspec.GridSpec(14, 8)  # create a GridSpec object
    #gs1.update(wspace=.09, hspace=.0)

    # 2D not-weighted gaussian convolved histogram, smallest bin width.
    ax0 = plt.subplot(gs1[0:2, 0:2])
    plt.xlabel('x (bins)', fontsize=12)
    plt.ylabel('y (bins)', fontsize=12)
    ax0.minorticks_on()
    plt.axvline(x=bin_center[0], linestyle='--', color='white')
    plt.axhline(y=bin_center[1], linestyle='--', color='white')
    # Radius
    circle = plt.Circle((bin_center[0], bin_center[1]),
        clust_rad / bin_list[0], color='w', fill=False)
    fig.gca().add_artist(circle)
    # Add text boxs.
    text = 'Bin: %.1f px' % (bin_list[0])
    plt.text(0.7, 0.94, text, transform=ax0.transAxes,
             bbox=dict(facecolor='white', alpha=0.8), fontsize=10)
    plt.imshow(h_filter.transpose(), origin='lower', aspect='auto')

    # 2D not-weighted histograms' centers.
    ax1 = plt.subplot(gs1[0:2, 2:4])
    # Get max and min values in x,y
    x_min, x_max = min(x_data), max(x_data)
    y_min, y_max = min(y_data), max(y_data)
    #Set plot limits
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    plt.xlabel('x (px)', fontsize=12)
    plt.ylabel('y (px)', fontsize=12)
    ax1.minorticks_on()
    # Add lines through meadian values with std deviations.
    plt.axvline(x=cent_stats[1][0], linestyle='-', color='k')
    plt.axvline(x=cent_stats[1][0] + cent_stats[2][0], linestyle='--',
        color='k')
    plt.axvline(x=cent_stats[1][0] - cent_stats[2][0], linestyle='--',
        color='k')
    plt.axhline(y=cent_stats[1][1], linestyle='-', color='k')
    plt.axhline(y=cent_stats[1][1] + cent_stats[2][1], linestyle='--',
        color='k')
    plt.axhline(y=cent_stats[1][1] - cent_stats[2][1], linestyle='--',
        color='k')
    # Add stats box.
    text1 = r'$(\tilde{x},\, \tilde{y}) = (%.1f, %.1f)\,px$' '\n' % \
    (cent_stats[1][0], cent_stats[1][1])
    text2 = '$(\sigma_x,\, \sigma_y) = (%.1f, %.1f)\,px$' % \
    (cent_stats[2][0], cent_stats[2][1])
    text = text1 + text2
    plt.text(0.05, 0.9, text, transform=ax1.transAxes,
        bbox=dict(facecolor='white', alpha=0.8), fontsize=11)
    cols = ['red', 'blue', 'green', 'black']
    for i in range(len(bin_list)):
        boxes = plt.gca()
        boxes.add_patch(Rectangle(((centers_kde[i][0] - bin_list[i]),
            (centers_kde[i][1] - bin_list[i])), bin_list[i] * 2.,
            bin_list[i] * 2., facecolor='none', edgecolor=cols[i], ls='solid',
            lw=1.5, zorder=(len(bin_list) - i),
            label='Bin: %.1f px' % bin_list[i]))
    # get handles
    handles, labels = ax1.get_legend_handles_labels()
    # use them in the legend
    leg1 = ax1.legend(handles, labels, loc='upper right', numpoints=1,
        fontsize=7)
    # Set the alpha value of the legend.
    leg1.get_frame().set_alpha(0.5)

    # x,y finding chart of full frame
    ax4 = plt.subplot(gs1[2:4, 0:2])
    #Set plot limits
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    #Set axis labels
    plt.xlabel('x (px)', fontsize=12)
    plt.ylabel('y (px)', fontsize=12)
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
    e_cent = cent_stats[0]
    text1 = '$x_{cent} = %.1f \pm %.1f px$' '\n' % (center_cl[0], e_cent)
    text2 = '$y_{cent} = %.1f \pm %.1f px$' % (center_cl[1], e_cent)
    text = text1 + text2
    plt.text(0.05, 0.9, text, transform=ax4.transAxes,
        bbox=dict(facecolor='white', alpha=0.85), fontsize=11)
    # Plot stars.
    # Solve for optimal star size.
    #from scipy.optimize import fsolve
    #a, c = 75., 2.5
    #area = (max(x_data) - min(x_data)) * (max(y_data) - min(y_data))
    #b = fsolve(star_size, 0.01, args=(a, c, area))
    #st_sizes_arr = a * np.exp(b * mag_data ** c)
    st_sizes_arr = 0.1 + 100. * 10 ** ((np.array(mag_data) -
        min(mag_data)) / -2.5)
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
    plt.xlabel('radius (px)', fontsize=12)
    plt.ylabel("stars/px$^{2}$", fontsize=12)
    # Set grid
    ax5.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
    # Cluster's name.
    text = str(clust_name)
    plt.text(0.4, 0.9, text, transform=ax5.transAxes, fontsize=14)
    # Legend texts
    kp_text = '3P' if flag_3pk_conver else '2P'
    texts = ['RDP (%0.1f px)' % bin_list[0],
            '$d_{field}$ = %.1E $st/px^{2}$' % field_dens,
            '%s King profile' % kp_text,
            'r$_c$ = %0.1f $\pm$ %0.1f px' % (rc, e_rc),
            'r$_t$ = %0.1f $\pm$ %0.1f px' % (rt, e_rt),
            'r$_{cl}$ = %0.1f $\pm$ %0.1f px' % (clust_rad, e_rad)]
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
    x_min, x_max = max(x_min, (center_cl[0] - 1.5 * clust_rad)), \
    min(x_max, (center_cl[0] + 1.5 * clust_rad))
    y_min, y_max = max(y_min, (center_cl[1] - 1.5 * clust_rad)), \
    min(y_max, (center_cl[1] + 1.5 * clust_rad))
    # Prevent axis stretching.
    if (x_max - x_min) != (y_max - y_min):
        lst = [(x_max - x_min), (y_max - y_min)]
        val, idx = min((val, idx) for (idx, val) in enumerate(lst))
        if idx == 0:
            x_max = x_min + lst[1]
        else:
            y_max = y_min + lst[0]
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    #Set axis labels
    plt.xlabel('x (px)', fontsize=12)
    plt.ylabel('y (px)', fontsize=12)
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
    # Plor contour levels.
    # Get KDE for CMD intrinsic position of most probable members.
    #x, y = np.mgrid[x_min:x_max:100j, y_min:y_max:100j]
    #positions = np.vstack([x.ravel(), y.ravel()])
    #x_zoom, y_zoom = [], []
    #for indx, star_x in enumerate(x_data):
        #if x_min < star_x < x_max and y_min < y_data[indx] < y_max:
            #x_zoom.append(star_x)
            #y_zoom.append(y_data[indx])
    #values = np.vstack([x_zoom, y_zoom])
    ## The results are HEAVILY dependant on the bandwidth used here.
    #kernel = stats.gaussian_kde(values)

    #k_pos = kernel(positions)
    ## Print x,y coordinates of max value.
    #new_cent = positions.T[np.argmax(k_pos)]
    #print new_cent
    ext_range, x, y, k_pos = kde_pl
    kde = np.reshape(k_pos.T, x.shape)
    plt.imshow(np.rot90(kde), cmap=plt.cm.YlOrBr, extent=ext_range)
    plt.contour(x, y, kde, 10, colors='k', linewidths=0.6)
    # Plot stars.
    plt.scatter(x_data, y_data, marker='o', c='black', s=st_sizes_arr, zorder=4)
    #Plot center.
    plt.scatter(center_cl[0], center_cl[1], color='w', s=40, lw=0.8,
        marker='x', zorder=5)
    #plt.scatter(new_cent[0], new_cent[1], color='g', s=40, lw=0.8,
        #marker='x', zorder=5)

    # Cluster and field regions defined.
    ax7 = plt.subplot(gs1[4:6, 0:2])
    # Get max and min values in x,y
    x_min, x_max = min(x_data), max(x_data)
    y_min, y_max = min(y_data), max(y_data)
    #Set plot limits
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    #Set axis labels
    plt.xlabel('x (px)', fontsize=12)
    plt.ylabel('y (px)', fontsize=12)
    # Set minor ticks
    ax7.minorticks_on()
    ax7.grid(b=True, which='both', color='gray', linestyle='--', lw=0.5)
    # Radius
    circle = plt.Circle((center_cl[0], center_cl[1]), clust_rad,
                        color='k', fill=False)
    fig.gca().add_artist(circle)
    plt.text(0.4, 0.92, 'Cluster + %d Field regions' % (len(field_region)),
             transform=ax7.transAxes,
             bbox=dict(facecolor='white', alpha=0.8), fontsize=12)
    # Plot cluster region.
    clust_reg_temp = [[], []]
    for star in cluster_region:
        dist = np.sqrt((center_cl[0] - star[1]) ** 2 +
        (center_cl[1] - star[2]) ** 2)
        # Only plot stars inside the cluster's radius.
        if dist <= clust_rad:
            clust_reg_temp[0].append(star[1])
            clust_reg_temp[1].append(star[2])
    plt.scatter(clust_reg_temp[0], clust_reg_temp[1], marker='o', c='red',
                s=8, edgecolors='none')
    if not flag_area_stronger:
        # Plot field stars regions.
        col = cycle(['DimGray', 'ForestGreen', 'maroon', 'RoyalBlue'])
        for i, reg in enumerate(field_region):
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
    plt.xlim(col1_min, col1_max)
    plt.ylim(mag_min, mag_max)
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
    for fr in field_region:
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
    plt.xlim(col1_min, col1_max)
    plt.ylim(mag_min, mag_max)
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
    tot_stars = len(stars_in_rjct) + len(stars_in)
    plt.text(0.55, 0.93, '$r \leq r_{cl}\,|\,N=%d$' % tot_stars,
             transform=ax9.transAxes,
             bbox=dict(facecolor='white', alpha=0.5), fontsize=16)
    # Plot stars.
    stars_rjct_temp = [[], []]
    for star in stars_in_rjct:
        stars_rjct_temp[0].append(star[5])
        stars_rjct_temp[1].append(star[3])
    plt.scatter(stars_rjct_temp[0], stars_rjct_temp[1], marker='x', c='teal',
                s=12, zorder=1)
    stars_acpt_temp = [[], []]
    for star in stars_in:
        stars_acpt_temp[0].append(star[5])
        stars_acpt_temp[1].append(star[3])
    sz_pt = 0.5 if (len(stars_in_rjct) + len(stars_in)) > 1000 else 1.
    plt.scatter(stars_acpt_temp[0], stars_acpt_temp[1], marker='o', c='k',
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
    # Plot stars.
    stars_rjct_temp = [[], []]
    for star in stars_out_rjct:
        stars_rjct_temp[0].append(star[3])
        stars_rjct_temp[1].append(star[4])
    for star in stars_in_rjct:
        stars_rjct_temp[0].append(star[3])
        stars_rjct_temp[1].append(star[4])
    plt.scatter(stars_rjct_temp[0], stars_rjct_temp[1], marker='x', c='teal',
                s=15, zorder=1)
    stars_acpt_temp = [[], []]
    for star in stars_out:
        stars_acpt_temp[0].append(star[3])
        stars_acpt_temp[1].append(star[4])
    for star in stars_in:
        stars_acpt_temp[0].append(star[3])
        stars_acpt_temp[1].append(star[4])
    plt.scatter(stars_acpt_temp[0], stars_acpt_temp[1], marker='o', c='k',
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
    # Plot stars.
    stars_rjct_temp = [[], []]
    for star in stars_out_rjct:
        stars_rjct_temp[0].append(star[3])
        stars_rjct_temp[1].append(star[6])
    for star in stars_in_rjct:
        stars_rjct_temp[0].append(star[3])
        stars_rjct_temp[1].append(star[6])
    plt.scatter(stars_rjct_temp[0], stars_rjct_temp[1], marker='x', c='teal',
                s=15, zorder=1)
    stars_acpt_temp = [[], []]
    for star in stars_out:
        stars_acpt_temp[0].append(star[3])
        stars_acpt_temp[1].append(star[6])
    for star in stars_in:
        stars_acpt_temp[0].append(star[3])
        stars_acpt_temp[1].append(star[6])
    plt.scatter(stars_acpt_temp[0], stars_acpt_temp[1], marker='o', c='k',
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
    ax12.vlines(x=mag_peak, ymin=0., ymax=plt.ylim()[1], color='k',
        lw=1.5, linestyles='dashed',
        label='$' + y_ax + r'_{compl}\,\approx\,%0.1f$' % mag_peak, zorder=1)
    # Legends.
    leg11 = plt.legend(fancybox=True, loc='upper right', numpoints=1,
                       fontsize=13)
    # Set the alpha value of the legend.
    leg11.get_frame().set_alpha(0.7)

    # Integrated magnitude and color.
    ax13 = plt.subplot(gs1[6:8, 2:4])
    # If field lists are not empty.
    if fl_reg_mag[0].any() and fl_reg_col[0].any():
        x_min = min(min(cl_reg_mag[0]), min(fl_reg_mag[0]),
            min(cl_reg_col[0]), min(fl_reg_col[0])) - 0.2
        x_max = max(max(cl_reg_mag[0]), max(fl_reg_mag[0]),
            max(cl_reg_col[0]), max(fl_reg_col[0])) + 0.2
        y_min = max(max(cl_reg_mag[1]), max(fl_reg_mag[1]),
            max(cl_reg_col[1]), max(fl_reg_col[1])) + 0.2
        y_max = min(min(cl_reg_mag[1]), min(fl_reg_mag[1]),
            min(cl_reg_col[1]), min(fl_reg_col[1])) - 0.2
    else:
        x_min, x_max = min(min(cl_reg_mag[0]), min(cl_reg_col[0])) - 0.2,\
        max(max(cl_reg_mag[0]), max(cl_reg_col[0])) + 0.2
        y_min, y_max = max(max(cl_reg_mag[1]), max(cl_reg_col[1])) + 0.2,\
        min(min(cl_reg_mag[1]), min(cl_reg_col[1])) - 0.2
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    ax13.set_xlabel('$mag$', fontsize=18)
    ax13.set_ylabel('$mag^*$', fontsize=18)
    ax13.minorticks_on()
    ax13.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
    # System used.
    if cmd_select == 1:
        x_ax0 = 'B'
    elif cmd_select == 2:
        x_ax0 = 'I'
    elif cmd_select == 3:
        x_ax0 = 'U'
    elif cmd_select == 4:
        x_ax0 = 'C'
    text1 = '$' + y_ax + '^{*}_{cl+fl}$'
    text2 = '$' + x_ax0 + '^{*}_{cl+fl}$'
    # Cluster + field integrated magnitude curve.
    plt.plot(cl_reg_mag[0], cl_reg_mag[1], 'r-', lw=1., label=text1)
    # Cluster integrated magnitude.
    plt.plot(cl_reg_col[0], cl_reg_col[1], 'r:', lw=2., label=text2)
    # Check if field regiones were defined.
    if not flag_area_stronger:
        text3 = '$' + y_ax + '^{*}_{fl}$'
        text4 = '$' + x_ax0 + '^{*}_{fl}$'
        # Field average integrated magnitude curve.
        plt.plot(fl_reg_mag[0], fl_reg_mag[1], 'b-', lw=1., label=text3)
        # Field average integrated magnitude.
        plt.plot(fl_reg_col[0], fl_reg_col[1], 'b:', lw=2., label=text4)
    text = '$(' + x_ax0 + '^{*} -' + y_ax + '^{*} )_{cl} = %0.2f$' % \
    (integ_col - integ_mag)
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

    # Norm fit for decontamination algorithm probability values.
    plot_colorbar = False
    ax16 = plt.subplot(gs1[8:10, 0:2])
    plt.xlim(0., 1.)
    plt.xlabel('membership probability', fontsize=12)
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
    # applied.
    # Used for the finding chart with colors assigned according to the
    # probabilities obtained.
    # Check if decont algorithm was applied.
    ax17 = plt.subplot(gs1[8:10, 2:4])
    # Get max and min values in x,y
    x_min, x_max = 10000., -10000
    y_min, y_max = 10000., -10000
    for star in cluster_region:
        x_min, x_max = min(star[1], x_min), max(star[1], x_max)
        y_min, y_max = min(star[2], y_min), max(star[2], y_max)
    #Set plot limits
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
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
    star_size = 20 if field_dens > 0.005 else 35
    m_p_m_temp = [[], [], []]
    for star in memb_prob_avrg_sort:
        m_p_m_temp[0].append(star[1])
        m_p_m_temp[1].append(star[2])
        m_p_m_temp[2].append(star[7])
    # Create new list with inverted values so higher prob stars are on top.
    m_p_m_temp_inv = [i[::-1] for i in m_p_m_temp]
    plt.scatter(m_p_m_temp_inv[0], m_p_m_temp_inv[1], marker='o',
                c=m_p_m_temp_inv[2], s=star_size, edgecolors='black',
                cmap=cm, lw=0.5)
    out_clust_rad = [[], []]
    for star in cluster_region:
        dist = np.sqrt((center_cl[0] - star[1]) ** 2 +
        (center_cl[1] - star[2]) ** 2)
        # Only plot stars outside the cluster's radius.
        if dist >= clust_rad:
            out_clust_rad[0].append(star[1])
            out_clust_rad[1].append(star[2])
    plt.scatter(out_clust_rad[0], out_clust_rad[1], marker='o',
                s=star_size, edgecolors='black', facecolors='none', lw=0.5)

    # Star's membership probabilities on cluster's CMD.
    ax18 = plt.subplot(gs1[8:10, 4:6])
    #Set plot limits
    plt.xlim(col1_min, col1_max)
    plt.ylim(mag_min, mag_max)
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
    if bf_flag:
        # Plot isochrone if best fit process was used.
        plt.plot(shift_isoch[0], shift_isoch[1], 'g', lw=1.2)
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
    sca = plt.scatter(m_p_m_temp_inv[0], m_p_m_temp_inv[1], marker='o',
                c=m_p_m_temp_inv[2], s=40, cmap=cm, lw=0.5, vmin=v_min_mp,
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
        plt.xlim(col1_min, col1_max)
        plt.ylim(mag_min, mag_max)
        #Set axis labels
        plt.xlabel('$' + x_ax + '$', fontsize=18)
        plt.ylabel('$' + y_ax + '$', fontsize=18)
        # Set minor ticks
        ax19.minorticks_on()
        ax19.xaxis.set_major_locator(MultipleLocator(1.0))
        ax19.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
        # Add text box
        m, a, e, d = isoch_fit_params[0]
        e_m, e_a, e_e, e_d = isoch_fit_errors
        text1 = '$z = %0.4f \pm %0.4f$' '\n' % (m, e_m)
        text2 = '$log(age) = %0.2f \pm %0.2f$' '\n' % (a, e_a)
        text3 = '$E_{(B-V)} = %0.2f \pm %0.2f$' '\n' % (e, e_e)
        text4 = '$(m-M)_o = %0.2f \pm %0.2f$' % (d, e_d)
        text = text1 + text2 + text3 + text4
        plt.text(0.5, 0.77, text, transform=ax19.transAxes,
                 bbox=dict(facecolor='white', alpha=0.6), fontsize=12)
        # Plot isochrone.
        plt.plot(shift_isoch[0], shift_isoch[1], 'r', lw=1.2)
        # Plot synth clust.
        plt.scatter(synth_clst[0], synth_clst[2], marker='o', s=30,
                    c='#4682b4', lw=0.5)

    # Best fitting process plots for GA.
    if bf_flag and best_fit_algor == 'genet':

        # Set ranges used by plots below.
        m_min, m_max = m_rs[:2]
        a_min, a_max = a_rs[:2]
        e_min, e_max = e_rs[:2]
        d_min, d_max = d_rs[:2]
        if m_min == m_max:
            m_min, m_max = m_min - 0.1 * m_min, m_max + 0.1 * m_min
        if a_min == a_max:
            a_min, a_max = a_min - 0.1 * a_min, a_max + 0.1 * a_min
        if e_min == e_max:
            e_min, e_max = e_min - 0.1 * e_min, e_max + 0.1 * e_min
        if d_min == d_max:
            d_min, d_max = d_min - 0.1 * d_min, d_max + 0.1 * d_min

        # Age vs metallicity GA diagram.
        isoch_done = isoch_fit_params[3]
        plt.subplot(gs1[10:12, 0:2])
        # Axis limits.
        plt.xlim(m_min, m_max)
        plt.ylim(a_min, a_max)
        plt.xlabel('$z$', fontsize=16)
        plt.ylabel('$log(age)$', fontsize=16)
        plt.minorticks_on()
        # Plot best fit point.
        plt.scatter(m, a, marker='o', c='r', s=30)
        # Plot ellipse error.
        ax20 = plt.gca()
        ellipse = Ellipse(xy=(m, a), width=2 * e_m, height=2 * e_a,
                                edgecolor='r', fc='None', lw=1.)
        ax20.add_patch(ellipse)
        # Plot density map.
        hist, xedges, yedges = np.histogram2d(zip(*isoch_done[0])[0],
                                              zip(*isoch_done[0])[1], bins=100)
        # H_g is the 2D histogram with a gaussian filter applied
        h_g = gaussian_filter(hist, 2, mode='constant')
        plt.imshow(h_g.transpose(), origin='lower',
                   extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
                   cmap=plt.get_cmap('Blues'), aspect='auto')

        # GA diagram.
        lkl_old, ext_imm_indx, isoch_done = isoch_fit_params[1], \
        isoch_fit_params[2], isoch_fit_params[3]
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
        (len(isoch_done[0]), N_b)
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
        # Plot extinction/immigration lines.
        for lin in ext_imm_indx:
            plt.axvline(x=lin, linestyle='--', lw=2., color='green')
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
        # Plot ellipse error.
        ax21 = plt.gca()
        ellipse = Ellipse(xy=(e, d), width=2 * e_e, height=2 * e_d,
                                edgecolor='b', fc='None', lw=1.)
        ax21.add_patch(ellipse)
        # Plot density map.
        hist, xedges, yedges = np.histogram2d(zip(*isoch_done[0])[2],
                                              zip(*isoch_done[0])[3], bins=100)
        # H_g is the 2D histogram with a gaussian filter applied
        h_g = gaussian_filter(hist, 2, mode='constant')
        plt.imshow(h_g.transpose(), origin='lower',
                   extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
                   cmap=plt.get_cmap('Reds'), aspect='auto')

        ax22 = plt.subplot(gs1[12:14, 0:2])
        plt.ylim(max(0, min(lkl_old[0]) - 0.3 * min(lkl_old[0])),
            max(lkl_old[1]))
        plt.xlim(m_min, m_max)
        # Set minor ticks
        ax22.minorticks_on()
        ax22.tick_params(axis='y', which='major', labelsize=9)
        plt.ylabel('Likelihood', fontsize=12)
        plt.xlabel('$z$', fontsize=16)
        text = '$z = %0.4f \pm %0.4f$' % (m, e_m)
        plt.text(0.1, 0.93, text, transform=ax22.transAxes,
            bbox=dict(facecolor='white', alpha=0.5), fontsize=12)
        hist, xedges, yedges = np.histogram2d(zip(*isoch_done[0])[0],
                                              isoch_done[1], bins=100)
        # H_g is the 2D histogram with a gaussian filter applied
        h_g = gaussian_filter(hist, 2, mode='constant')
        y_min_edge = max(0, min(lkl_old[0]) - 0.3 * min(lkl_old[0]))
        plt.imshow(h_g.transpose(), origin='lower',
                   extent=[xedges[0], xedges[-1], y_min_edge, yedges[-1]],
                   cmap=plt.get_cmap('gist_yarg'), aspect='auto')
        plt.axvline(x=m, linestyle='--', color='blue')
        plt.axvline(x=m + e_m, linestyle='--', color='red')
        plt.axvline(x=m - e_m, linestyle='--', color='red')

        ax23 = plt.subplot(gs1[12:14, 2:4])
        plt.ylim(max(0, min(lkl_old[0]) - 0.3 * min(lkl_old[0])),
            max(lkl_old[1]))
        plt.xlim(a_min, a_max)
        # Set minor ticks
        ax23.minorticks_on()
        ax23.tick_params(axis='y', which='major', labelsize=9)
        plt.ylabel('Likelihood', fontsize=12)
        plt.xlabel('$log(age/yr)$', fontsize=16)
        text = '$log(age/yr) = %0.2f \pm %0.2f$' % (a, e_a)
        plt.text(0.1, 0.93, text, transform=ax23.transAxes,
            bbox=dict(facecolor='white', alpha=0.5), fontsize=12)
        hist, xedges, yedges = np.histogram2d(zip(*isoch_done[0])[1],
                                              isoch_done[1], bins=100)
        # H_g is the 2D histogram with a gaussian filter applied
        h_g = gaussian_filter(hist, 2, mode='constant')
        plt.imshow(h_g.transpose(), origin='lower',
                   extent=[xedges[0], xedges[-1], y_min_edge, yedges[-1]],
                   cmap=plt.get_cmap('gist_yarg'), aspect='auto')
        plt.axvline(x=a, linestyle='--', color='blue')
        plt.axvline(x=a + e_a, linestyle='--', color='red')
        plt.axvline(x=a - e_a, linestyle='--', color='red')

        ax24 = plt.subplot(gs1[12:14, 4:6])
        plt.ylim(max(0, min(lkl_old[0]) - 0.3 * min(lkl_old[0])),
            max(lkl_old[1]))
        plt.xlim(e_min, e_max)
        # Set minor ticks
        ax24.minorticks_on()
        ax24.tick_params(axis='y', which='major', labelsize=9)
        plt.ylabel('Likelihood', fontsize=12)
        plt.xlabel('$E_{(B-V)}$', fontsize=16)
        text = '$E_{(B-V)} = %0.2f \pm %0.2f$' % (e, e_e)
        plt.text(0.1, 0.93, text, transform=ax24.transAxes,
            bbox=dict(facecolor='white', alpha=0.5), fontsize=12)
        hist, xedges, yedges = np.histogram2d(zip(*isoch_done[0])[2],
                                              isoch_done[1], bins=100)
        # H_g is the 2D histogram with a gaussian filter applied
        h_g = gaussian_filter(hist, 2, mode='constant')
        plt.imshow(h_g.transpose(), origin='lower',
                   extent=[xedges[0], xedges[-1], y_min_edge, yedges[-1]],
                   cmap=plt.get_cmap('gist_yarg'), aspect='auto')
        plt.axvline(x=e, linestyle='--', color='blue')
        plt.axvline(x=e + e_e, linestyle='--', color='red')
        plt.axvline(x=e - e_e, linestyle='--', color='red')

        ax25 = plt.subplot(gs1[12:14, 6:8])
        plt.ylim(max(0, min(lkl_old[0]) - 0.3 * min(lkl_old[0])),
            max(lkl_old[1]))
        plt.xlim(d_min, d_max)
        # Set minor ticks
        ax25.minorticks_on()
        ax25.tick_params(axis='y', which='major', labelsize=9)
        plt.ylabel('Likelihood', fontsize=12)
        plt.xlabel('$(m-M)_o$', fontsize=16)
        text = '$(m-M)_o = %0.2f \pm %0.2f$' % (d, e_d)
        plt.text(0.1, 0.93, text, transform=ax25.transAxes,
            bbox=dict(facecolor='white', alpha=0.5), fontsize=12)
        hist, xedges, yedges = np.histogram2d(zip(*isoch_done[0])[3],
                                              isoch_done[1], bins=100)
        # H_g is the 2D histogram with a gaussian filter applied
        h_g = gaussian_filter(hist, 2, mode='constant')
        plt.imshow(h_g.transpose(), origin='lower',
                   extent=[xedges[0], xedges[-1], y_min_edge, yedges[-1]],
                   cmap=plt.get_cmap('gist_yarg'), aspect='auto')
        plt.axvline(x=d, linestyle='--', color='blue')
        plt.axvline(x=d + e_d, linestyle='--', color='red')
        plt.axvline(x=d - e_d, linestyle='--', color='red')

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
    pl_fmt = pl_params[1]
    pl_dpi = pl_params[2]
    plt.savefig(join(output_subdir, str(clust_name) + '.' + pl_fmt), dpi=pl_dpi)

    # Close to release memory.
    plt.clf()
    plt.close()
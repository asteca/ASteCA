
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.offsetbox as offsetbox
from itertools import cycle
from matplotlib.patches import Rectangle
from ..errors import error_round
from ..structure import king_prof_funcs as kpf


def pl_centers(gs, x_name, y_name, coord, x_min, x_max, y_min, y_max,
               asp_ratio, approx_cents, bin_width, st_dev_lst):
    '''
    2D Gaussian histograms' centers using different standard deviations.
    '''
    ax = plt.subplot(gs[0:2, 0:2])
    ax.set_aspect(aspect=asp_ratio)
    # Set plot limits.
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    if coord == 'deg':
        # If RA is used, invert axis.
        ax.invert_xaxis()
    plt.xlabel('{} ({})'.format(x_name, coord), fontsize=12)
    plt.ylabel('{} ({})'.format(y_name, coord), fontsize=12)
    ax.minorticks_on()
    # Add lines through median values with std deviations.
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
    cent_med_r, cent_std_dev_r = error_round.round_sig_fig(
        cent_median, cent_std_dev)
    text1 = (r"$(\tilde{{{0}}},\,\tilde{{{1}}}) = ({3:g},\,{4:g})\,"
             "{2}$").format(x_name, y_name, coord, *cent_med_r)
    text2 = ("$(\sigma_{{{0}}},\,\sigma_{{{1}}}) = ({3:g},\,{4:g})\,"
             "{2}$").format(x_name, y_name, coord, *cent_std_dev_r)
    text = text1 + '\n' + text2
    ob = offsetbox.AnchoredText(text, loc=2, prop=dict(size=11))
    ob.patch.set(alpha=0.85)
    ax.add_artist(ob)
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
    handles, labels = ax.get_legend_handles_labels()
    # use them in the legend
    leg = ax.legend(handles, labels, loc='lower right', numpoints=1,
                    fontsize=7)
    # Set the alpha value of the legend.
    leg.get_frame().set_alpha(0.5)


def pl_hist_g(gs, fig, asp_ratio, x_name, y_name, coord, cent_bin, clust_rad,
              bin_width, hist_2d_g):
    '''
    2D Gaussian convolved histogram.
    '''
    bin_w_r = error_round.round_to_y(bin_width)

    ax = plt.subplot(gs[0:2, 2:4])
    plt.xlabel('{} (bins)'.format(x_name), fontsize=12)
    plt.ylabel('{} (bins)'.format(y_name), fontsize=12)
    ax.minorticks_on()
    plt.axvline(x=cent_bin[0], linestyle='--', color='white')
    plt.axhline(y=cent_bin[1], linestyle='--', color='white')
    # Radius
    circle = plt.Circle((cent_bin[0], cent_bin[1]), clust_rad / bin_width,
                        color='w', fill=False)
    fig.gca().add_artist(circle)
    # Add text box.
    text = 'Bin $\simeq$ {0:g} {1}'.format(bin_w_r, coord)
    ob = offsetbox.AnchoredText(text, pad=0.2, loc=1, prop=dict(size=10))
    ob.patch.set(alpha=0.85)
    ax.add_artist(ob)
    plt.imshow(hist_2d_g.transpose(), origin='lower')
    # If RA is used, invert axis.
    if coord == 'deg':
        ax.invert_xaxis()
    ax.set_aspect(aspect=asp_ratio)


def pl_rad_dens(gs, mode, radii, rdp_points, field_dens, coord, clust_name,
                clust_rad, e_rad, poisson_error, bin_width, core_rad, e_core,
                tidal_rad, e_tidal, K_cent_dens, flag_2pk_conver,
                flag_3pk_conver):
    '''
    Radial density plot.
    '''
    bin_w_r = error_round.round_to_y(bin_width)

    ax = plt.subplot(gs[0:2, 4:8])
    # Get max and min values in x,y
    x_min, x_max = max(min(radii) - (max(radii) / 20.), 0), \
        max(radii) + (max(radii) / 20.)
    delta_total = (max(rdp_points) - field_dens)
    delta_backg = 0.2 * delta_total
    y_min = max((field_dens - delta_backg) - (max(rdp_points) -
                min(rdp_points)) / 10, 0)
    y_max = max(rdp_points) + (max(rdp_points) - min(rdp_points)) / 10
    # Set plot limits
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    ax.minorticks_on()
    # Set axes labels
    plt.xlabel('radius ({})'.format(coord), fontsize=12)
    plt.ylabel("stars/{}$^{{2}}$".format(coord), fontsize=12)
    # Set grid
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
    # Cluster's name and mode used to process it.
    plt.title(str(clust_name) + ' (' + mode + ')', fontsize=14)
    # Round radii values.
    rads_r, e_rads_r = error_round.round_sig_fig(
        [core_rad, tidal_rad, clust_rad], [e_core, e_tidal, e_rad])
    # Legend texts
    kp_text = '3P' if flag_3pk_conver else '2P'
    texts = [
        'RDP ($\sim${0:g} {1})'.format(bin_w_r, coord),
        '$d_{{field}}$ = {:.1E} $st/{}^{{2}}$'.format(field_dens, coord),
        '{} King profile'.format(kp_text),
        'r$_c$ = {0:g} $\pm$ {1:g} {2}'.format(rads_r[0], e_rads_r[0], coord),
        'r$_t$ = {0:g} $\pm$ {1:g} {2}'.format(rads_r[1], e_rads_r[1], coord),
        'r$_{{cl}}$ = {0:g} $\pm$ {1:g} {2}'.format(rads_r[2], e_rads_r[2],
                                                    coord)
    ]
    # Plot density profile with the smallest bin size
    ax.plot(radii, rdp_points, 'ko-', zorder=3, label=texts[0])
    # Plot Poisson error bars
    plt.errorbar(radii, rdp_points, yerr=poisson_error, fmt='ko',
                 zorder=1)
    # Plot background level.
    ax.hlines(y=field_dens, xmin=0, xmax=max(radii), label=texts[1],
              color='b', zorder=5)
    # Approx middle of the graph.
    arr_y_up = (y_max - y_min) / 2.3 + y_min
    # Length and width of arrow head.
    head_w, head_l = x_max * 0.023, (y_max - y_min) * 0.045
    # Length of arrow.
    arr_y_dwn = -1. * abs(arr_y_up - field_dens) * 0.76
    # Plot 3-P King profile.
    if flag_3pk_conver:
        # Plot curve.
        ax.plot(radii, kpf.three_params(
            radii, tidal_rad, K_cent_dens, core_rad, field_dens),
            'g--', label=texts[2], lw=2., zorder=3)
        # Plot r_t radius as an arrow. vline is there to show the label.
        ax.vlines(x=tidal_rad, ymin=0., ymax=0., label=texts[4], color='g')
        ax.arrow(tidal_rad, arr_y_up, 0., arr_y_dwn, fc="g", ec="g",
                 head_width=head_w, head_length=head_l, zorder=5)
        # Plot r_c as a dashed line.
        ax.vlines(x=core_rad, ymin=0, ymax=kpf.three_params(
            core_rad, tidal_rad, K_cent_dens, core_rad, field_dens),
            label=texts[3], color='g', linestyles=':', lw=4., zorder=4)
    # Plot 2-P King profile if 3-P was not found.
    elif flag_2pk_conver:
        # Plot curve.
        ax.plot(radii, kpf.two_params(
            radii, K_cent_dens, core_rad, field_dens), 'g--', label=texts[2],
            lw=2., zorder=3)
        # Plot r_c as a dashed line.
        ax.vlines(x=core_rad, ymin=0, ymax=kpf.two_params(
            core_rad, K_cent_dens, core_rad, field_dens), label=texts[3],
            color='g', linestyles=':', lw=4., zorder=4)
    # Plot radius.
    ax.vlines(x=clust_rad, ymin=0, ymax=0., label=texts[5], color='r')
    ax.arrow(clust_rad, arr_y_up, 0., arr_y_dwn, fc="r", ec="r",
             head_width=head_w, head_length=head_l, zorder=5)
    # Plot radius error zone.
    if e_rad > 0.:
        plt.axvspan((clust_rad - e_rad), (clust_rad + e_rad),
                    facecolor='grey', alpha=0.5)
    # get handles
    handles, labels = ax.get_legend_handles_labels()
    # use them in the legend
    ax.legend(handles, labels, loc='upper right', numpoints=2, fontsize=11)


def pl_full_frame(
        gs, fig, x_name, y_name, coord, x_min, x_max, y_min, y_max, asp_ratio,
        clust_cent, clust_rad, e_cent, x, y, st_sizes_arr, core_rad, e_core,
        tidal_rad, e_tidal, K_conct_par, flag_2pk_conver, flag_3pk_conver):
    '''
    x,y finding chart of full frame
    '''
    ax = plt.subplot(gs[0:2, 8:10])
    ax.set_aspect(aspect=asp_ratio)
    # Set plot limits
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    # If RA is used, invert axis.
    if coord == 'deg':
        ax.invert_xaxis()
    # Set axis labels
    plt.xlabel('{} ({})'.format(x_name, coord), fontsize=12)
    plt.ylabel('{} ({})'.format(y_name, coord), fontsize=12)
    # Set minor ticks
    ax.minorticks_on()
    # Plot r_cl.
    circle = plt.Circle((clust_cent[0], clust_cent[1]), clust_rad, color='r',
                        fill=False, lw=1.5)
    fig.gca().add_artist(circle)
    if flag_3pk_conver is True:
        # Plot tidal radius.
        circle = plt.Circle(
            (clust_cent[0], clust_cent[1]), tidal_rad, color='g', fill=False,
            lw=1.5)
        fig.gca().add_artist(circle)
        # Plot core radius.
        if core_rad > 0:
            circle = plt.Circle(
                (clust_cent[0], clust_cent[1]), core_rad, color='g',
                fill=False, ls='dashed', lw=1.)
            fig.gca().add_artist(circle)
    elif flag_2pk_conver is True:
        # Plot core radius.
        if core_rad > 0:
            circle = plt.Circle(
                (clust_cent[0], clust_cent[1]), core_rad, color='g',
                fill=False, ls='dashed', lw=1.)
            fig.gca().add_artist(circle)
    # Add text box
    center_cl_r, e_cent_r = error_round.round_sig_fig(clust_cent, e_cent)
    text1 = '${0}_{{cent}} = {1:g} \pm {2:g}\,{3}$'.format(
        x_name, center_cl_r[0], e_cent_r[0], coord)
    text2 = '${0}_{{cent}} = {1:g} \pm {2:g}\,{3}$'.format(
        y_name, center_cl_r[1], e_cent_r[1], coord)
    text = text1 + '\n' + text2
    ob = offsetbox.AnchoredText(text, pad=0.2, loc=2, prop=dict(size=11))
    ob.patch.set(alpha=0.85)
    ax.add_artist(ob)
    # Plot stars.
    plt.scatter(x, y, marker='o', c='black', s=st_sizes_arr)


def pl_zoom_frame(gs, fig, x_name, y_name, coord, x_zmin, x_zmax, y_zmin,
                  y_zmax, cont_index, kde_plot, x_data, y_data, st_sizes_arr,
                  center_cl, clust_rad):
    '''
    Zoom on x,y finding chart.
    '''
    ax = plt.subplot(gs[0:2, 10:12])
    # Force square plot.
    ax.set_aspect('equal')
    # Set plot limits.
    plt.xlim(x_zmin, x_zmax)
    plt.ylim(y_zmin, y_zmax)
    if coord == 'deg':
        # If RA is used, invert axis.
        ax.invert_xaxis()
    # Set axis labels
    plt.xlabel('{} ({})'.format(x_name, coord), fontsize=12)
    plt.ylabel('{} ({})'.format(y_name, coord), fontsize=12)
    # Set minor ticks
    ax.minorticks_on()
    # Plot radius.
    circle = plt.Circle((center_cl[0], center_cl[1]), clust_rad, color='r',
                        fill=False, lw=1.5, zorder=5)
    fig.gca().add_artist(circle)
    # Add text box.
    text = 'Cluster zoom\nCI = %0.2f' % (cont_index)
    ob = offsetbox.AnchoredText(text, pad=0.2, loc=1, prop=dict(size=12))
    ob.patch.set(alpha=0.85)
    ax.add_artist(ob)
    # Plot contour levels if it was obtained.
    if kde_plot:
        ext_range, x, y, k_pos = kde_plot
        # Number of countour lines depends on how large the plotted area is
        # compared with the area where the posotional KDE was obtained.
        frac_xy = (x_zmax - x_zmin) / (ext_range[1] - ext_range[0])
        if frac_xy > 2.:
            c_lines = 10
        elif 1. < frac_xy <= 2.:
            c_lines = 15
        elif 0.5 < frac_xy <= 1.:
            c_lines = 20
        elif 0.2 < frac_xy <= 0.5:
            c_lines = 25
        elif frac_xy <= 0.2:
            c_lines = 30
        kde = np.reshape(k_pos.T, x.shape)
        plt.imshow(np.rot90(kde), cmap=plt.cm.YlOrBr, extent=ext_range)
        plt.contour(x, y, kde, c_lines, colors='b', linewidths=0.6)
    # Plot stars.
    plt.scatter(x_data, y_data, marker='o', c='black', s=st_sizes_arr,
                zorder=1)
    # Plot center.
    plt.scatter(center_cl[0], center_cl[1], color='w', s=40, lw=0.8,
                marker='x', zorder=5)


def pl_cl_fl_regions(gs, fig, x_name, y_name, coord, x_min, x_max, y_min,
                     y_max, asp_ratio, center_cl, clust_rad, field_regions,
                     cl_region, flag_no_fl_regs):
    '''
    Cluster and field regions defined.
    '''
    ax = plt.subplot(gs[2:4, 0:2])
    ax.set_aspect(aspect=asp_ratio)
    # Set plot limits
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    # If RA is used, invert axis.
    if coord == 'deg':
        ax.invert_xaxis()
    # Set axis labels
    plt.xlabel('{} ({})'.format(x_name, coord), fontsize=12)
    plt.ylabel('{} ({})'.format(y_name, coord), fontsize=12)
    # Set minor ticks
    ax.minorticks_on()
    ax.grid(b=True, which='both', color='gray', linestyle='--', lw=0.5)
    # Radius
    circle = plt.Circle((center_cl[0], center_cl[1]), clust_rad,
                        color='k', fill=False)
    fig.gca().add_artist(circle)
    # Add text box.
    text = 'Cluster + %d Field regions' % (len(field_regions))
    ob = offsetbox.AnchoredText(text, pad=0.2, loc=1, prop=dict(size=12))
    ob.patch.set(alpha=0.85)
    ax.add_artist(ob)
    # Plot cluster region.
    plt.scatter(zip(*cl_region)[1], zip(*cl_region)[2], marker='o', c='red',
                s=8, edgecolors='none')
    if not flag_no_fl_regs:
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


def plot(N, *args):
    '''
    Handle each plot separately.
    '''

    plt_map = {
        0: [pl_hist_g, 'density positional chart'],
        1: [pl_centers, 'obtained centers in positional chart'],
        2: [pl_full_frame, 'full frame'],
        3: [pl_rad_dens, 'radial density function'],
        4: [pl_zoom_frame, 'zoomed frame'],
        5: [pl_cl_fl_regions, 'cluster and field regions defined']
    }

    fxn = plt_map.get(N, None)[0]
    if fxn is None:
        raise ValueError("  ERROR: there is no plot {}.".format(N))

    try:
        fxn(*args)
    except:
        print("  WARNING: error when plotting {}.".format(plt_map.get(N)[1]))
        # import traceback
        # print traceback.format_exc()

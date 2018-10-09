
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.offsetbox as offsetbox
from itertools import cycle
from ..structure import king_prof_funcs as kpf


def pl_center(gs, fig, asp_ratio, x_name, y_name, coord, bin_cent, clust_rad,
              bin_width, hist_2d_g):
    '''
    2D Gaussian convolved histogram.
    '''

    ax = plt.subplot(gs[0:2, 0:2])
    plt.xlabel('{} (bins)'.format(x_name), fontsize=12)
    plt.ylabel('{} (bins)'.format(y_name), fontsize=12)
    ax.minorticks_on()
    plt.axvline(x=bin_cent[0], linestyle='--', color='green')
    plt.axhline(y=bin_cent[1], linestyle='--', color='green')
    # Radius
    circle = plt.Circle(
        (bin_cent[0], bin_cent[1]), clust_rad / bin_width, color='green',
        fill=False)
    fig.gca().add_artist(circle)
    # Add text box.
    text = 'Bin $\simeq$ {0:g} {1}'.format(round(bin_width, 1), coord)
    ob = offsetbox.AnchoredText(text, pad=0.2, loc=1, prop=dict(size=10))
    ob.patch.set(alpha=0.85)
    ax.add_artist(ob)
    plt.imshow(hist_2d_g.transpose(), origin='lower',
               cmap=plt.get_cmap('RdYlBu_r'))
    plt.contour(hist_2d_g.transpose(), 5, colors='#551a8b', linewidths=0.5)
    # If RA is used, invert axis.
    if coord == 'deg':
        ax.invert_xaxis()
    ax.set_aspect(aspect=asp_ratio)


def pl_rad_dens(gs, mode, radii, rdp_points, field_dens, coord, clust_name,
                clust_rad, e_rad, poisson_error, bin_width, core_rad,
                e_core, tidal_rad, e_tidal, K_cent_dens, flag_2pk_conver,
                flag_3pk_conver):
    '''
    Radial density plot.
    '''
    ax = plt.subplot(gs[0:2, 2:6])
    # Get max and min values in x,y
    x_min, x_max = max(min(radii) - (max(radii) / 20.), 0), \
        min(max(radii) + (max(radii) / 20.), 3. * clust_rad)
    delta_total = (max(rdp_points) - field_dens)
    delta_backg = 0.2 * delta_total
    y_min = max((field_dens - delta_backg) -
                (max(rdp_points) - min(rdp_points)) / 10, 0)
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
    plt.title(str(clust_name) + ' (' + mode + ')', fontsize=9)
    # Legend texts
    kp_text = '3P' if flag_3pk_conver else '2P'
    texts = [
        'RDP ($\sim${:.0f} {})'.format(bin_width, coord),
        '$d_{{field}}$ = {:.1E} $st/{}^{{2}}$'.format(field_dens, coord),
        '{} King profile'.format(kp_text),
        'r$_c$ = {0:.0f} $\pm$ {1:.0f} {2}'.format(core_rad, e_core, coord),
        'r$_t$ = {0:.0f} $\pm$ {1:.0f} {2}'.format(tidal_rad, e_tidal, coord),
        'r$_{{cl}}$ = {0:.0f} $\pm$ {1:.0f} {2}'.format(
            clust_rad, e_rad, coord)
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
        kde_cent, clust_rad, x, y, st_sizes_arr, core_rad, e_core,
        tidal_rad, e_tidal, K_conct_par, flag_2pk_conver, flag_3pk_conver):
    '''
    x,y finding chart of full frame
    '''
    ax = plt.subplot(gs[2:4, 0:2])
    ax.set_aspect(aspect=asp_ratio)
    ax.set_title(
        r"$N_{{stars}}$={} (incomp frame)".format(len(x)), fontsize=9)
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
    circle = plt.Circle((kde_cent[0], kde_cent[1]), clust_rad, color='r',
                        fill=False, lw=1.5)
    fig.gca().add_artist(circle)
    if flag_3pk_conver is True:
        # Plot tidal radius.
        circle = plt.Circle(
            (kde_cent[0], kde_cent[1]), tidal_rad, color='g', fill=False,
            lw=1.5)
        fig.gca().add_artist(circle)
        # Plot core radius.
        if core_rad > 0:
            circle = plt.Circle(
                (kde_cent[0], kde_cent[1]), core_rad, color='g',
                fill=False, ls='dashed', lw=1.)
            fig.gca().add_artist(circle)
    elif flag_2pk_conver is True:
        # Plot core radius.
        if core_rad > 0:
            circle = plt.Circle(
                (kde_cent[0], kde_cent[1]), core_rad, color='g',
                fill=False, ls='dashed', lw=1.)
            fig.gca().add_artist(circle)
    # Add text box
    text1 = '${0}_{{cent}} = {1:.0f}\,{2}$'.format(x_name, kde_cent[0], coord)
    text2 = '${0}_{{cent}} = {1:.0f}\,{2}$'.format(y_name, kde_cent[1], coord)
    text = text1 + '\n' + text2
    ob = offsetbox.AnchoredText(text, pad=0.2, loc=2, prop=dict(size=11))
    ob.patch.set(alpha=0.85)
    ax.add_artist(ob)
    # Plot stars.
    plt.scatter(x, y, marker='o', c='black', s=st_sizes_arr)


def pl_zoom_frame(gs, fig, x_name, y_name, coord, x_zmin, x_zmax, y_zmin,
                  y_zmax, cont_index, kde_plot, x_data, y_data, st_sizes_arr,
                  kde_cent, clust_rad):
    '''
    Zoom on x,y finding chart.
    '''
    ax = plt.subplot(gs[2:4, 2:4])
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
    circle = plt.Circle((kde_cent[0], kde_cent[1]), clust_rad, color='r',
                        fill=False, lw=1.5, zorder=5)
    fig.gca().add_artist(circle)
    # Add text box.
    text = "CI = {:.2f}".format(cont_index)
    ob = offsetbox.AnchoredText(text, pad=0.2, loc=1, prop=dict(size=12))
    ob.patch.set(alpha=0.85)
    ax.add_artist(ob)
    # Plot contour levels if it was obtained.
    if kde_plot:
        ext_range, x, y, k_pos = kde_plot
        # Number of contour lines depends on how large the plotted area is
        # compared with the area where the positional KDE was obtained.
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
    plt.scatter(kde_cent[0], kde_cent[1], color='w', s=40, lw=0.8,
                marker='x', zorder=5)


def pl_cl_fl_regions(gs, fig, x_name, y_name, coord, x_min, x_max, y_min,
                     y_max, asp_ratio, kde_cent, clust_rad, field_regions_i,
                     field_regions_rjct_i, cl_region_i, cl_region_rjct_i,
                     flag_no_fl_regs_i):
    '''
    Cluster and field regions defined.
    '''
    ax = plt.subplot(gs[2:4, 4:6])
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
    circle = plt.Circle((kde_cent[0], kde_cent[1]), clust_rad,
                        color='k', fill=False)
    fig.gca().add_artist(circle)
    # Add text box.
    text = 'Cluster + %d Field regions' % (len(field_regions_i))
    ob = offsetbox.AnchoredText(text, pad=0.2, loc=1, prop=dict(size=12))
    ob.patch.set(alpha=0.85)
    ax.add_artist(ob)
    # Plot cluster region.
    if len(cl_region_rjct_i) > 0:
        plt.scatter(
            list(zip(*cl_region_rjct_i))[1], list(zip(*cl_region_rjct_i))[2],
            marker='x', c='orange', s=5, lw=.5, edgecolors='none')
    plt.scatter(list(zip(*cl_region_i))[1], list(zip(*cl_region_i))[2],
                marker='o', c='red', s=8, edgecolors='none')

    N_flrg = 0
    if not flag_no_fl_regs_i:
        col = cycle(['DimGray', 'ForestGreen', 'maroon', 'RoyalBlue'])
        # Stars inside the field regions with accepted errors.
        for i, reg in enumerate(field_regions_i):
            fl_reg = list(zip(*reg))
            N_flrg += len(fl_reg[0])
            plt.scatter(fl_reg[1], fl_reg[2], marker='o',
                        c=next(col), s=8, edgecolors='none')
        # Stars inside the field regions with rejected errors.
        for i, reg in enumerate(field_regions_rjct_i):
            fl_reg = list(zip(*reg))
            N_flrg += len(fl_reg[0])
            plt.scatter(fl_reg[1], fl_reg[2], marker='x',
                        c=next(col), s=5, lw=.5, edgecolors='none')

    ax.set_title(r"$N_{{stars}}$={} (incomp frame)".format(
        len(cl_region_i) + len(cl_region_rjct_i) + N_flrg), fontsize=9)


def plot(N, *args):
    '''
    Handle each plot separately.
    '''

    plt_map = {
        0: [pl_center, 'density positional chart'],
        1: [pl_full_frame, 'full frame'],
        2: [pl_rad_dens, 'radial density function'],
        3: [pl_zoom_frame, 'zoomed frame'],
        4: [pl_cl_fl_regions, 'cluster and field regions defined']
    }

    fxn = plt_map.get(N, None)[0]
    if fxn is None:
        raise ValueError("  ERROR: there is no plot {}.".format(N))

    try:
        fxn(*args)
    except Exception:
        print("  WARNING: error when plotting {}.".format(plt_map.get(N)[1]))
        import traceback
        print(traceback.format_exc())

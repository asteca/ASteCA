
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.offsetbox as offsetbox
from mpl_toolkits.axes_grid1 import make_axes_locatable
from itertools import cycle
from ..structure import king_prof_funcs as kpf


def pl_center(gs, fig, asp_ratio, x_name, y_name, coord, bw_list, kde_cent,
              frame_kde_cent, kde_dens_max, kde_dens_min, clust_rad):
    '''
    2D Gaussian convolved histogram.
    '''

    ax = plt.subplot(gs[0:2, 0:2])
    frmt = '{:.4f}' if coord == 'deg' else '{:.0f}'
    ax.set_title((r'$KDE_{{bdw}}$ =' + frmt + ' [{}]').format(
        bw_list[1], coord), fontsize=9)
    plt.xlabel('{} ({})'.format(x_name, coord), fontsize=12)
    plt.ylabel('{} ({})'.format(y_name, coord), fontsize=12)
    ax.minorticks_on()
    plt.axvline(x=kde_cent[0], linestyle='--', color='green')
    plt.axhline(y=kde_cent[1], linestyle='--', color='green')
    # Radius
    circle = plt.Circle(
        (kde_cent[0], kde_cent[1]), clust_rad, color='green', fill=False)
    fig.gca().add_artist(circle)

    ext_range, x_grid, y_grid, k_pos = frame_kde_cent
    kde = np.reshape(k_pos.T, x_grid.shape)
    im = plt.imshow(
        np.rot90(kde), cmap=plt.get_cmap('RdYlBu_r'), extent=ext_range)
    plt.contour(x_grid, y_grid, kde, colors='#551a8b', linewidths=0.5)
    # If RA is used, invert axis.
    if coord == 'deg':
        ax.invert_xaxis()

    # Colorbar on the right side of ax.
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.05)
    cbar = plt.colorbar(im, cax=cax)
    cbar.set_ticks([np.min(kde), np.ptp(kde) * .5, np.max(kde)])
    scale = 3600. if coord == 'deg' else 1.

    midpt = ((kde_dens_max + kde_dens_min) * .5) / scale
    frmt = '{:.2E}' if midpt > 100. or midpt < .1 else '{:.0f}'
    cbar.ax.set_yticklabels([
        frmt.format(kde_dens_min / scale),
        frmt.format(midpt),
        frmt.format(kde_dens_max / scale)], rotation=90)
    cbar.ax.tick_params(labelsize=9)
    # Align bottom and middle labels. Don't include last label (one at the
    # top) as we don't want to change its alignment.
    for i, label in enumerate(cbar.ax.get_yticklabels()[:-1]):
        if i == 0:
            label.set_va("bottom")
        else:
            label.set_va("center")

    ax.set_aspect(aspect=asp_ratio)


def pl_rad_dens(
    gs, mode, radii, rdp_points, field_dens, coord, clust_name, clust_rad,
        e_rad, poisson_error, bin_width, core_rad, e_core, tidal_rad, e_tidal,
        K_cent_dens, flag_2pk_conver, flag_3pk_conver):
    '''
    Radial density plot.
    '''
    # Convert from deg to arcmin if (ra,dec) were used.
    if coord == 'deg':
        radii = np.array(radii) * 60.
        clust_rad, e_rad = clust_rad * 60., e_rad * 60.
        bin_width = bin_width * 60.
        core_rad, e_core, tidal_rad, e_tidal = core_rad * 60., e_core * 60.,\
            tidal_rad * 60., e_tidal * 60.
        field_dens, K_cent_dens = field_dens / 3600., K_cent_dens / 3600.
        rdp_points = np.array(rdp_points) / 3600.
        poisson_error = np.array(poisson_error) / 3600.
        coord2 = 'arcmin'
    else:
        coord2 = 'px'

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
    plt.xlabel(r'radius $[{}]$'.format(coord2), fontsize=12)
    plt.ylabel(r"N $[st/{}^{{2}}]$".format(coord2), fontsize=12)
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=.5)
    # Cluster's name and mode used to process it.
    plt.title(str(clust_name) + ' (' + mode + ')', fontsize=12)
    # Legend texts
    kp_text = '3P' if flag_3pk_conver else '2P'
    r_frmt = '{:.0f}' if coord2 == 'px' else '{:.1f}'
    texts = [
        (r'RDP ($bin_{{w}}\approx$' + r_frmt + r' $[{}]$)').format(
            bin_width, coord2),
        r"$d_{{fld}}=$ {:.1E} $[st/{}^{{2}}]$".format(field_dens, coord2),
        '{} King profile'.format(kp_text),
        (r"$r_c=$" + r_frmt + r"$\pm$" + r_frmt + r" $[{}]$").format(
            core_rad, e_core, coord2),
        (r"$r_t=$" + r_frmt + r"$\pm$" + r_frmt + r" $[{}]$").format(
            tidal_rad, e_tidal, coord2),
        (r"$r_{{cl}}=$" + r_frmt + r"$\pm$" + r_frmt + r' $[{}]$').format(
            clust_rad, e_rad, coord2)
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
                 head_width=head_w, head_length=head_l, zorder=4)
        # Plot r_c as a dashed line.
        ax.vlines(x=core_rad, ymin=0, ymax=kpf.three_params(
            core_rad, tidal_rad, K_cent_dens, core_rad, field_dens),
            label=texts[3], color='g', linestyles=':', lw=4., zorder=5)
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
    ax.legend(handles, labels, loc='upper right', numpoints=2, fontsize=10)


def pl_full_frame(
    gs, fig, x_name, y_name, coord, x_min, x_max, y_min, y_max, asp_ratio,
        kde_cent, clust_rad, x, y, st_sizes_arr, core_rad, e_core, tidal_rad,
        e_tidal, K_conct_par, flag_2pk_conver, flag_3pk_conver):
    '''
    x,y finding chart of full frame
    '''
    ax = plt.subplot(gs[2:4, 0:2])
    ax.set_aspect(aspect=asp_ratio)
    ax.set_title(
        r"$N_{{stars}}$={} (phot incomp)".format(len(x)), fontsize=9)
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
    r_frmt = '{:.0f}' if coord == 'px' else '{:.3f}'
    t1 = (r'${}_{{cent}} =$' + r_frmt + '$\,{}$').format(
        x_name, kde_cent[0], coord)
    t2 = (r'${}_{{cent}} =$' + r_frmt + '$\,{}$').format(
        y_name, kde_cent[1], coord)
    text = t1 + '\n' + t2
    ob = offsetbox.AnchoredText(text, pad=0.2, loc=2, prop=dict(size=11))
    ob.patch.set(alpha=0.85)
    ax.add_artist(ob)
    # Plot stars.
    plt.scatter(x, y, marker='o', c='black', s=st_sizes_arr)


def pl_zoom_frame(
    gs, fig, x_name, y_name, coord, x_zmin, x_zmax, y_zmin, y_zmax, cont_index,
        kde_plot, x_data, y_data, st_sizes_arr, kde_cent, clust_rad):
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
    N_r = [1 for pt in list(zip(x_data, y_data)) if circle.contains_point(pt)]
    ax.set_title(r"$N_{{r<r_{{cl}}}}$={} [$CI = {:.2f}$] (phot incomp)".format(
        sum(N_r), cont_index), fontsize=9)
    fig.gca().add_artist(circle)

    # Plot contour levels if it was obtained.
    if kde_plot:
        ext_range, x, y, k_pos = kde_plot
        # Number of contour lines depends on how large the cluster area is
        # compared with the area where the positional KDE was obtained.
        frac_xy = (ext_range[1] - ext_range[0]) / (2. * clust_rad)
        if frac_xy <= 1.5:
            c_lines = 5
        elif 1.5 < frac_xy <= 2.:
            c_lines = 10
        elif 2. < frac_xy <= 2.5:
            c_lines = 15
        elif 2.5 < frac_xy <= 3.:
            c_lines = 20
        elif 3 < frac_xy <= 3.5:
            c_lines = 25
        elif 3.5 < frac_xy:
            c_lines = 30

        kde = np.reshape(k_pos.T, x.shape)
        plt.contour(x, y, kde, c_lines, colors='b', linewidths=0.6)
    # Plot stars.
    plt.scatter(x_data, y_data, marker='o', c='black', s=st_sizes_arr,
                zorder=1)
    # Plot center.
    plt.scatter(kde_cent[0], kde_cent[1], color='r', s=40, lw=0.8,
                marker='x', zorder=5)


def pl_cl_fl_regions(
    gs, fig, x_name, y_name, coord, x_min, x_max, y_min, y_max, asp_ratio,
        kde_cent, clust_rad, field_regions_i, field_regions_rjct_i,
        cl_region_i, cl_region_rjct_i, flag_no_fl_regs_i):
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
    ax.grid(b=True, which='both', color='gray', linestyle='--', lw=.5)
    # Radius
    circle = plt.Circle((kde_cent[0], kde_cent[1]), clust_rad,
                        color='k', fill=False)
    fig.gca().add_artist(circle)

    # Plot cluster region.
    if len(cl_region_rjct_i) > 0:
        plt.scatter(
            list(zip(*cl_region_rjct_i))[1], list(zip(*cl_region_rjct_i))[2],
            marker='o', c='red', s=8, edgecolors='w', lw=.2)
    plt.scatter(list(zip(*cl_region_i))[1], list(zip(*cl_region_i))[2],
                marker='o', c='red', s=8, edgecolors='w', lw=.2)

    N_flrg = 0
    if not flag_no_fl_regs_i:
        col0 = cycle(['DimGray', 'ForestGreen', 'maroon', 'RoyalBlue'])
        col1 = cycle(['DimGray', 'ForestGreen', 'maroon', 'RoyalBlue'])
        # Stars inside the field regions with accepted errors.
        for i, reg in enumerate(field_regions_i):
            fl_reg = list(zip(*reg))
            N_flrg += len(fl_reg[0])
            plt.scatter(fl_reg[1], fl_reg[2], marker='o',
                        c=next(col0), s=8, edgecolors='w', lw=.2)
        # Stars inside the field regions with rejected errors.
        for i, reg in enumerate(field_regions_rjct_i):
            if reg:
                fl_reg = list(zip(*reg))
                N_flrg += len(fl_reg[0])
                plt.scatter(fl_reg[1], fl_reg[2], marker='o',
                            c=next(col1), s=8, edgecolors='w', lw=.2)

    ax.set_title(r"$N_{{stars}}$={} (phot incomp); $N_{{fregs}}$={}".format(
        len(cl_region_i) + len(cl_region_rjct_i) + N_flrg,
        len(field_regions_i)), fontsize=9)


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

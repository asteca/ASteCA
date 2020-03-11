
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.offsetbox as offsetbox
from mpl_toolkits.axes_grid1 import make_axes_locatable
from . prep_plots import xylabelsize, xytickssize, titlesize, cbartickssize,\
    legendsize


def pl_densmap(
    gs, fig, asp_ratio, x_name, y_name, coord, bw_list, kde_cent,
        frame_kde_cent, fr_dens, clust_rad):
    """
    Coordinates 2D KDE.
    """

    x_name, y_name = r'$\alpha_{2000}$', r'$\delta_{2000}$'

    ax = plt.subplot(gs[0:2, 0:2])
    frmt = '{:.4f}' if coord == 'deg' else '{:.0f}'
    # ax.set_title((r'$KDE_{{bdw}}$ =' + frmt + ' [{}]').format(
    #     bw_list[1], coord), fontsize=titlesize)
    plt.xlabel('{} ({})'.format(x_name, coord), fontsize=xylabelsize)
    plt.ylabel('{} ({})'.format(y_name, coord), fontsize=xylabelsize)

    ax.minorticks_on()
    ax.tick_params(axis='both', which='major', labelsize=8) #xytickssize)

    # plt.axvline(x=kde_cent[0], linestyle='--', color='green')
    # plt.axhline(y=kde_cent[1], linestyle='--', color='green')

    # Radius
    circle = plt.Circle(
        (kde_cent[0], kde_cent[1]), clust_rad, color='green', fill=False)
    ax.add_artist(circle)

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
    # cbar.set_ticks([np.min(kde), np.ptp(kde) * .5, np.max(kde)])
    cbar.set_ticks([np.min(kde), np.max(kde)])
    scale = 3600. if coord == 'deg' else 1.

    kde_dens_min, kde_dens_max = fr_dens.min(), fr_dens.max()
    midpt = ((kde_dens_max + kde_dens_min) * .5) / scale
    frmt = '{:.2E}' if midpt > 100. or midpt < .1 else '{:.0f}'
    cbar.ax.set_yticklabels([
        frmt.format(kde_dens_min / scale),
        # frmt.format(midpt),
        frmt.format(kde_dens_max / scale)], rotation=0) #90)
    cbar.ax.tick_params(labelsize=cbartickssize)
    # # Align bottom and middle labels. Don't include last label (one at the
    # # top) as we don't want to change its alignment.
    # for i, label in enumerate(cbar.ax.get_yticklabels()[:-1]):
    #     if i == 0:
    #         label.set_va("bottom")
    #     else:
    #         label.set_va("center")

    # Add text box
    kde_cent = (186.085049, -61.8505119)
    r_frmt = '{:.0f}' if coord == 'px' else '{:.5f}'
    x_cent = kde_cent[0]
    x_name, y_name = r'\alpha', r'\delta'
    t1 = (r'${}_{{c}} =$' + r_frmt).format(x_name, x_cent)
    t2 = (r'${}_{{c}} =$' + r_frmt).format(y_name, kde_cent[1])
    text = '   NGC4349\n' + t1 + '\n' + t2
    ob = offsetbox.AnchoredText(
        text, pad=0.2, loc=2, prop=dict(size=legendsize))
    ob.patch.set(alpha=0.85)
    ax.add_artist(ob)

    import matplotlib.patches as patches
    x_outer_min, x_outer_max = -0.09, 0.109
    y_outer_min, y_outer_max = -0.085, 0.115
    x_inner_min, x_inner_max = -0.07, 0.085
    y_inner_min, y_inner_max = -0.065, 0.09
    w_in, h_in = x_inner_max - x_inner_min, y_inner_max - y_inner_min
    rect = patches.Rectangle(
        (x_inner_min, y_inner_min), w_in, h_in, linewidth=2., linestyle='--',
        edgecolor='k', facecolor='none')
    ax.add_patch(rect)
    w_ou, h_ou = x_outer_max - x_outer_min, y_outer_max - y_outer_min
    rect = patches.Rectangle(
        (x_outer_min, y_outer_min), w_ou, h_ou, linewidth=2., linestyle='--',
        edgecolor='k', facecolor='none')
    ax.add_patch(rect)

    ax.set_aspect(aspect=asp_ratio)


def pl_knn_dens(
    gs, fig, asp_ratio, x_min, x_max, y_min, y_max, x_name, y_name, coord,
        NN_dd, xy_filtered, fr_dens, NN_dist, kde_cent, clust_rad):
    """
    """
    ax = plt.subplot(gs[0:2, 2:4])
    ax.set_aspect(aspect=asp_ratio)
    ax.tick_params(axis='both', which='major', labelsize=xytickssize)
    # Set plot limits
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    # If RA is used, invert axis.
    if coord == 'deg':
        ax.invert_xaxis()
    ax.set_title(
        r'$kNN={}\;(d\leq d_{{p=25\%}})$'.format(NN_dd), fontsize=titlesize)

    plt.xlabel('{} ({})'.format(x_name, coord), fontsize=xylabelsize)
    plt.ylabel('{} ({})'.format(y_name, coord), fontsize=xylabelsize)
    ax.minorticks_on()
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=.5)

    perc = np.percentile(NN_dist, 25)
    msk = NN_dist < perc
    xy, NN_d = xy_filtered[msk], NN_dist[msk]
    for i, (x, y) in enumerate(xy):
        circle = plt.Circle(
            (x, y), NN_d[i], color='k', lw=.5, alpha=.5, fill=False)
        ax.add_artist(circle)

    # Star with the smallest associated density.
    idx = np.argmin(fr_dens)
    circle = plt.Circle(
        (xy_filtered[idx][0], xy_filtered[idx][1]), NN_dist[idx], color='b',
        lw=1.5, fill=False, zorder=4)
    fig.gca().add_artist(circle)
    plt.plot([], [], color='b', lw=2., label=r"$dens_{min}$")
    # Star with the largest associated density.
    idx = np.argmax(fr_dens)
    circle = plt.Circle(
        (xy_filtered[idx][0], xy_filtered[idx][1]), NN_dist[idx], color='r',
        lw=1., fill=False, zorder=5)
    fig.gca().add_artist(circle)
    plt.plot([], [], color='r', lw=2., label=r"$dens_{max}$")

    # Assigned center.
    plt.scatter(kde_cent[0], kde_cent[1], color='g', s=40, lw=1.5,
                marker='x', zorder=5)
    # Radius
    circle = plt.Circle(
        (kde_cent[0], kde_cent[1]), clust_rad, color='red', fill=False)
    ax.add_artist(circle)

    leg = plt.legend(
        fancybox=True, fontsize=legendsize, handlelength=1., loc='best')
    leg.get_frame().set_alpha(0.7)


def pl_full_frame(
    gs, fig, project, x_offset, y_offset, x_name, y_name, coord, x_min, x_max,
        y_min, y_max, asp_ratio, kde_cent, x, y, st_sizes_arr, clust_rad):
    """
    x,y finding chart of full frame
    """
    ax = plt.subplot(gs[0:2, 4:6])
    ax.minorticks_on()
    ax.set_aspect(aspect=asp_ratio)
    ax.set_title(
        r"$N_{{stars}}$={} (phot incomp)".format(len(x)), fontsize=titlesize)
    # Set plot limits
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    # If RA is used, invert axis.
    if coord == 'deg':
        ax.invert_xaxis()
    # Set axis labels
    plt.xlabel('{} ({})'.format(x_name, coord), fontsize=xylabelsize)
    plt.ylabel('{} ({})'.format(y_name, coord), fontsize=xylabelsize)
    ax.tick_params(axis='both', which='major', labelsize=xytickssize)

    N_max = 50000  # HARDCODED
    if len(x) > N_max:
        print("  WARNING: too many stars. Plotting {} random samples.".format(
            N_max))
        ids = np.random.choice(np.arange(len(x)), N_max, replace=False)
        x, y = x[ids], y[ids]

    # Plot stars.
    plt.scatter(x, y, marker='o', c='black', s=st_sizes_arr)
    # plt.axvline(x=kde_cent[0], linestyle='--', color='green')
    # plt.axhline(y=kde_cent[1], linestyle='--', color='green')
    plt.scatter(*kde_cent, marker='x', c='red', s=50)
    # Radius
    circle = plt.Circle(
        (kde_cent[0], kde_cent[1]), clust_rad, color='red', fill=False)
    ax.add_artist(circle)

    # Add text box
    r_frmt = '{:.0f}' if coord == 'px' else '{:.5f}'
    if coord == 'deg' and project:
        x_cent = (kde_cent[0] / np.cos(np.deg2rad(kde_cent[1] + y_offset))) +\
            x_offset
    else:
        x_cent = kde_cent[0]
    t1 = (r'${}_{{c}} =$' + r_frmt + r'$\,{}$').format(x_name, x_cent, coord)
    t2 = (r'${}_{{c}} =$' + r_frmt + r'$\,{}$').format(
        y_name, kde_cent[1] + y_offset, coord)
    text = t1 + '\n' + t2
    ob = offsetbox.AnchoredText(
        text, pad=0.2, loc=2, prop=dict(size=legendsize))
    ob.patch.set(alpha=0.85)
    ax.add_artist(ob)


def pl_field_dens(
    gs, coord, fdens_method, fr_dist, fr_dens, fdens_min_d, fdens_lst,
        fdens_std_lst, field_dens_d, field_dens, field_dens_std):
    """
    Field density values for different percentiles.
    """
    # Convert from deg to arcmin if (ra,dec) were used.
    if coord == 'deg':
        fr_dist, fdens_min_d, field_dens_d = np.array(fr_dist) * 60.,\
            np.array(fdens_min_d) * 60., field_dens_d * 60.
        fr_dens, fdens_lst, fdens_std_lst = [
            np.array(_) / 3600. for _ in (fr_dens, fdens_lst, fdens_std_lst)]
        field_dens = field_dens / 3600.
        coord2 = 'arcmin'
    else:
        coord2 = 'px'

    delta_y = np.ptp(fr_dens) * .1
    ymin = max(0., min(fr_dens) - delta_y)
    ymax = max(fr_dens) + delta_y

    ax = plt.subplot(gs[2:4, 0:2])
    ax.set_title(("Method: '{}'").format(fdens_method), fontsize=titlesize)
    plt.ylim(ymin, ymax)
    ax.minorticks_on()
    plt.xlabel(
        r'Distance to center $[{}]$'.format(coord2), fontsize=xylabelsize)
    plt.ylabel(r"$\rho$ $[st/{}^{{2}}]$".format(coord2), fontsize=xylabelsize)
    ax.tick_params(axis='both', which='major', labelsize=xytickssize)
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=.5)

    plt.scatter(fr_dist, fr_dens, c='k', s=5, alpha=.2, zorder=1)
    plt.errorbar(
        fdens_min_d, fdens_lst, yerr=fdens_std_lst, fmt='b', ms=25,
        ecolor='r', lw=1.2)

    t1 = r"$d_{{field}}=$ {:.1E} $[st/{}^{{2}}]$".format(
        field_dens, coord2)

    # Check if a manual value was used
    if not np.isnan(field_dens_d):
        plt.scatter(
            field_dens_d, field_dens, marker='o', s=25, c='g', label=t1,
            zorder=5)
    else:
        ax.hlines(
            field_dens, xmin=fdens_min_d[0], xmax=fdens_min_d[-1], color='g',
            label=t1)

    leg = plt.legend(fancybox=True, fontsize=legendsize, loc='upper right')
    leg.get_frame().set_alpha(0.7)


def pl_mag_cent(gs, coord, y_ax, integ_dists, integ_mags):
    """
    """
    # Convert from deg to arcmin if (ra,dec) were used.
    if coord == 'deg':
        integ_dists = np.array(integ_dists) * 60.
        coord2 = 'arcmin'
    else:
        coord2 = 'px'

    ax = plt.subplot(gs[2:4, 2:4])
    ax.minorticks_on()
    ax.tick_params(axis='both', which='major', labelsize=xytickssize)
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=.5)
    ax.set_title('Integrated magnitude vs distance', fontsize=titlesize)

    # Nstp = max(1, int(round(len(rdp_radii) / 10.)))
    # rdp_radii, rdp_mags = rdp_radii[::Nstp], rdp_mags[::Nstp]

    # w = (rdp_radii[-1] - rdp_radii[0]) / 15.
    # plt.boxplot(rdp_mags, positions=rdp_radii, widths=w)
    # plt.xlim(rdp_radii[0] - 1.5 * w, rdp_radii[-1] + 1.5 * w)
    # ax.set_xticks(rdp_radii)
    # if coord2 == "arcmin":
    #     tck = ["{:.1f}".format(_) for _ in rdp_radii]
    # else:
    #     tck = ["{:.0f}".format(_) for _ in rdp_radii]
    # ax.set_xticklabels(tck)

    plt.plot(integ_dists, integ_mags)

    plt.ylabel('$' + y_ax + '$' + r" $^{*}$", fontsize=xylabelsize)
    plt.xlabel("Distance to center [{}]".format(coord2), fontsize=xylabelsize)
    ax.invert_yaxis()


def pl_rdp_rings(
    gs, fig, asp_ratio, x_min, x_max, y_min, y_max, x_name, y_name, coord,
        kde_cent, rdp_radii):
    """
    """
    ax = plt.subplot(gs[2:4, 4:6])
    ax.minorticks_on()
    ax.set_aspect(aspect=asp_ratio)
    # Set plot limits
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    plt.xlabel('{} ({})'.format(x_name, coord), fontsize=xylabelsize)
    plt.ylabel('{} ({})'.format(y_name, coord), fontsize=xylabelsize)
    ax.tick_params(axis='both', which='major', labelsize=xytickssize)
    # If RA is used, invert axis.
    if coord == 'deg':
        ax.invert_xaxis()
    ax.set_title(
        r'$N_{{rings}}={}$'.format(len(rdp_radii)), fontsize=titlesize)

    for rad in rdp_radii:
        circle = plt.Circle(
            (kde_cent[0], kde_cent[1]), rad, color='g', fill=False, ls=':',
            zorder=5)
        ax.add_artist(circle)


def plot(N, *args):
    """
    Handle each plot separately.
    """

    plt_map = {
        0: [pl_densmap, 'density positional chart'],
        1: [pl_knn_dens, 'kNN per-star densities'],
        2: [pl_full_frame, 'full frame'],
        3: [pl_field_dens, 'Field density'],
        4: [pl_mag_cent, 'magnitude vs dist to center'],
        5: [pl_rdp_rings, 'RDP rings']
    }

    fxn = plt_map.get(N, None)[0]
    if fxn is None:
        raise ValueError("  ERROR: there is no plot {}.".format(N))

    try:
        fxn(*args)
    except Exception:
        print("  WARNING: error when plotting {}".format(plt_map.get(N)[1]))
        import traceback
        print(traceback.format_exc())


import numpy as np
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
import matplotlib.offsetbox as offsetbox
from mpl_toolkits.axes_grid1 import make_axes_locatable


def pl_full_frame(
    gs, fig, x_name, y_name, coord, x_min, x_max, y_min, y_max, asp_ratio,
        kde_cent, x, y, st_sizes_arr, clust_rad):
    """
    x,y finding chart of full frame
    """
    ax = plt.subplot(gs[0:4, 0:4])
    ax.set_aspect(aspect=asp_ratio)
    ax.set_title(r"$N_{{stars}}$={}".format(len(x)))
    # Set plot limits
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    ax.invert_xaxis()
    # Set axis labels
    plt.xlabel('{} ({})'.format(x_name, coord))
    plt.ylabel('{} ({})'.format(y_name, coord))

    N_max = 50000  # HARDCODED
    if len(x) > N_max:
        print("  WARNING: too many stars. Plotting {} random samples.".format(
            N_max))
        ids = np.random.choice(np.arange(len(x)), N_max, replace=False)
        x, y, st_sizes_arr = x[ids], y[ids], st_sizes_arr[ids]

    # Plot stars.
    plt.scatter(x, y, marker='o', c='black', s=st_sizes_arr * 1.5)
    # plt.axvline(x=kde_cent[0], linestyle='--', color='green')
    # plt.axhline(y=kde_cent[1], linestyle='--', color='green')
    plt.scatter(*kde_cent, marker='x', c='red', s=50)
    # Radius
    circle = plt.Circle(
        (kde_cent[0], kde_cent[1]), clust_rad, color='red', fill=False)
    ax.add_artist(circle)

    # Add text box
    r_frmt = '{:.5f}'
    t1 = (r'${}_{{c}} =$' + r_frmt + r'$\,{}$').format(
        x_name, kde_cent[0], coord)
    t2 = (r'${}_{{c}} =$' + r_frmt + r'$\,{}$').format(
        y_name, kde_cent[1], coord)
    text = t1 + '\n' + t2
    # prop={'fontsize': 15})
    ob = offsetbox.AnchoredText(text, pad=0.2, loc=2)
    ob.patch.set(alpha=0.85)
    ax.add_artist(ob)


def pl_densmap(
    gs, fig, asp_ratio, x_name, y_name, coord, bw_list, kde_cent,
        frame_kde_cent, fr_dens, clust_rad):
    """
    Coordinates 2D KDE.
    """

    ax = plt.subplot(gs[0:2, 4:6])
    ax.set_title((r'$KDE_{{bdw}}$ ={:.4f} [{}]').format(
        bw_list[1], coord))
    plt.xlabel('{} ({})'.format(x_name, coord))
    plt.ylabel('{} ({})'.format(y_name, coord))

    plt.axvline(x=kde_cent[0], linestyle='--', lw=1., color='white')
    plt.axhline(y=kde_cent[1], linestyle='--', lw=1., color='white')
    # Radius
    circle = plt.Circle(
        (kde_cent[0], kde_cent[1]), clust_rad, color='r', fill=False)
    ax.add_artist(circle)

    ext_range, x_grid, y_grid, k_pos = frame_kde_cent
    kde = np.reshape(k_pos.T, x_grid.shape)
    im = plt.imshow(
        np.rot90(kde), cmap=plt.get_cmap('RdYlBu_r'), extent=ext_range)
    plt.contour(x_grid, y_grid, kde, colors='#551a8b', linewidths=0.5)
    ax.invert_xaxis()

    # Colorbar on the right side of ax.
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.05)
    cbar = plt.colorbar(im, cax=cax)
    cbar.set_ticks([
        np.min(kde), (np.max(kde) + np.min(kde)) * .5, np.max(kde)])
    cbar.ax.minorticks_off()
    scale = 3600.

    kde_dens_min, kde_dens_max = fr_dens.min(), fr_dens.max()
    midpt = ((kde_dens_max + kde_dens_min) * .5) / scale
    frmt = '{:.2E}' if midpt > 100. else '{:.1f}'
    cbar.ax.set_yticklabels([
        frmt.format(kde_dens_min / scale),
        frmt.format(midpt),
        frmt.format(kde_dens_max / scale)], rotation=90)
    cbar.ax.tick_params()

    # Align bottom and middle labels. Don't include last label (one at the
    # top) as we don't want to change its alignment.
    for i, label in enumerate(cbar.ax.get_yticklabels()[:-1]):
        if i == 0:
            label.set_va("bottom")
        else:
            label.set_va("center")

    ax.set_aspect(aspect=asp_ratio)


def pl_knn_dens(
    gs, fig, plot_style, asp_ratio, x_min, x_max, y_min, y_max, x_name,
    y_name, coord, NN_dd, xy_filtered, fr_dens, NN_dist, kde_cent,
        clust_rad):
    """
    """
    ax = plt.subplot(gs[2:4, 4:6])
    ax.set_aspect(aspect=asp_ratio)
    # Set plot limits
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    ax.invert_xaxis()
    ax.set_title(r'$kNN={}\;|\;max(100; (d\leq d_{{p=25\%}}))$'.format(NN_dd))

    plt.xlabel('{} ({})'.format(x_name, coord))
    plt.ylabel('{} ({})'.format(y_name, coord))
    if plot_style == 'asteca':
        ax.grid()

    perc = np.percentile(NN_dist, 25)
    msk = NN_dist < perc
    xy, NN_d = xy_filtered[msk], NN_dist[msk]
    # 100 maximum
    xy, NN_d = xy[:100], NN_d[:100]
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
    # Add coords to legend
    r_frmt = '{:.0f}' if coord == 'px' else '{:.5f}'
    x_max_dens = xy_filtered[idx][0]
    y_max_dens = xy_filtered[idx][1]
    x_max_dens = (x_max_dens / np.cos(np.deg2rad(y_max_dens)))
    t1 = ('(' + r_frmt + ',\n' + r_frmt + ')').format(x_max_dens, y_max_dens)
    plt.plot([], [], color='r', lw=2., label=r"$dens_{max}$" + "\n" + t1)

    # Assigned center.
    plt.scatter(kde_cent[0], kde_cent[1], color='g', s=40, lw=1.5,
                marker='x', zorder=5)
    # Radius
    circle = plt.Circle(
        (kde_cent[0], kde_cent[1]), clust_rad, color='red', fill=False)
    ax.add_artist(circle)

    leg = plt.legend(fancybox=True, handlelength=1., loc='best')
    leg.get_frame().set_alpha(0.7)


def pl_field_dens(
    gs, plot_style, coord, fdens_method, fr_dist, fr_dens, fdens_min_d,
        fdens_lst, fdens_std_lst, field_dens):
    """
    Field density values for different percentiles.
    """
    # Convert from deg to arcmin
    fr_dist, fdens_min_d = np.array(fr_dist) * 60., np.array(fdens_min_d) * 60.
    fr_dens, fdens_lst, fdens_std_lst = [
        np.array(_) / 3600. for _ in (fr_dens, fdens_lst, fdens_std_lst)]
    field_dens = field_dens / 3600.
    coord2 = 'arcmin'

    idx = np.argmin(abs(field_dens - np.array(fdens_lst)))
    field_dens_d = fdens_min_d[idx]

    delta_y = np.ptp(fr_dens) * .1
    ymin = max(0., min(fr_dens) - delta_y)
    ymax = max(fr_dens) + delta_y

    ax = plt.subplot(gs[4:6, 0:4])
    # ax.set_title(("Method: '{}'").format(fdens_method))
    plt.ylim(ymin, ymax)
    plt.xlabel(r'Distance to center $[{}]$'.format(coord2))
    plt.ylabel(r"$\rho$ $[st/{}^{{2}}]$".format(coord2))
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    if plot_style == 'asteca':
        ax.grid()

    plt.scatter(fr_dist, fr_dens, c='k', s=5, alpha=.2, zorder=1)
    plt.errorbar(
        fdens_min_d, fdens_lst, yerr=fdens_std_lst, fmt='b', ms=25,
        ecolor='r', lw=1.2)

    t1 = r"$d_{{field}}=$ {:.3E} $[st/{}^{{2}}]$".format(
        field_dens, coord2)

    # Check if a manual value was used
    if fdens_method != 'a':
        ax.hlines(
            field_dens, xmin=fdens_min_d[0], xmax=fdens_min_d[-1], color='g',
            label=t1, zorder=5)
    else:
        plt.scatter(
            field_dens_d, field_dens, marker='o', s=25, c='g', label=t1,
            zorder=5)

    # from ..structure.king_profile import KingProf as kpf
    # kpf_xvals = np.linspace(0, KP_Bys_rt[1] * 60., 100)
    # kpf_yvals = (KP_cent_dens/3600.) * kpf(
    #     kpf_xvals, KP_Bys_rc[3] * 60., KP_Bys_rt[3] * 60.) + field_dens
    # ax.plot(kpf_xvals, kpf_yvals, 'g--', lw=2., zorder=3)

    leg = plt.legend(fancybox=True, loc='upper right')
    leg.get_frame().set_alpha(0.7)


def pl_centdist_vs_mag(
    gs, fig, plot_style, y_ax, coord, xi, yi, magi, kde_cent, clust_rad,
        integ_dists, integ_mags):
    """
    """
    # Distances of all stars to the center of the cluster.
    xy_cent_dist = cdist([kde_cent], np.array((xi, yi)).T)[0]

    # Convert from deg to arcmin
    clust_rad, xy_cent_dist = clust_rad * 60., xy_cent_dist * 60.
    integ_dists = np.array(integ_dists) * 60.
    coord2 = 'arcmin'

    ax = plt.subplot(gs[4:6, 4:6])
    if plot_style == 'asteca':
        ax.grid()
    # plt.xlim(xy_cent_dist.min(), min(xy_cent_dist.max(), clust_rad * 10.))

    plt.xlabel(r'Distance to center $[{}]$'.format(coord2))
    plt.ylabel('$' + y_ax + '$')

    msk = xy_cent_dist <= clust_rad
    plt.scatter(
        xy_cent_dist[~msk], magi[~msk], s=10, c='grey', alpha=.75,
        edgecolor='w', lw=.3, label=r"$r>r_{cl}$")
    plt.scatter(
        xy_cent_dist[msk], magi[msk], s=25, c='k',
        edgecolor='w', lw=.3, label=r"$r\leq r_{cl}$")
    ax.invert_yaxis()
    plt.legend(fancybox=True, loc='lower right')

    ax2 = ax.twinx()
    plt.plot(integ_dists, integ_mags, c='r', lw=3., label="Integ mag")

    plt.ylabel('$' + y_ax + '$' + r" $^{*}$")
    ax2.invert_yaxis()

    leg = plt.legend(fancybox=True, loc='upper center')
    leg.get_frame().set_alpha(0.7)


def plot(N, *args):
    """
    Handle each plot separately.
    """

    plt_map = {
        0: [pl_full_frame, 'full frame'],
        1: [pl_densmap, 'density positional chart'],
        2: [pl_knn_dens, 'kNN per-star densities'],
        3: [pl_field_dens, 'Field density'],
        4: [pl_centdist_vs_mag, 'Distance to center vs magnitude']
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

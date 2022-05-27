
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Ellipse, Rectangle
import numpy as np
from astropy.stats import sigma_clipped_stats


def pms_VPD_all(gs, plot_style, xlabel, coord, y_ax, allfr_PMs):
    """
    PMs for all the frame.
    """
    # This is set to '3' in pms_analysis.kde_2d()
    PM_KDE_std = 3

    ax = plt.subplot(gs[0:2, 0:2])
    ax.set_title(
        r'All stars in frame + KDE ($\mu\pm\sigma={:.2f}$)'.format(PM_KDE_std))
    if plot_style == 'asteca':
        ax.grid()
    plt.xlabel(xlabel)
    plt.ylabel(r"$\mu_{{\delta}} \, \mathrm{[mas/yr]}$")

    cmap = plt.cm.get_cmap('viridis')
    mmag = allfr_PMs['mmag']
    msk = np.argsort(mmag)[::-1]
    pmRA_all, pmDE_all, mag = allfr_PMs['pmRA'][msk], allfr_PMs['pmDE'][msk],\
        mmag[msk]
    plt.scatter(
        pmRA_all, pmDE_all, s=8, c=mag, cmap=cmap,
        edgecolor='k', lw=.2, zorder=3)
    cbar = plt.colorbar(pad=.01, fraction=.02, aspect=50)
    cbar.ax.minorticks_off()
    cbar.ax.invert_yaxis()
    cbar.set_label(y_ax)

    ra_mean, ra_median, ra_std = sigma_clipped_stats(pmRA_all)
    de_mean, de_median, de_std = sigma_clipped_stats(pmDE_all)
    plt.xlim(ra_median - PM_KDE_std * ra_std, ra_median + PM_KDE_std * ra_std)
    plt.ylim(de_median - PM_KDE_std * de_std, de_median + PM_KDE_std * de_std)


def pms_VPD_KDE_all(gs, xlabel, coord, y_ax, allr_KDE_PMs, xydelta, xyrang):
    """
    KDE for the PMS of all the stars in the frame.
    """
    PMx, PMy, PMz = allr_KDE_PMs['x'], allr_KDE_PMs['y'], allr_KDE_PMs['z']

    ax = plt.subplot(gs[0:2, 2:4])
    ax.set_title("All stars KDE (zoomed)")
    plt.xlabel(xlabel)
    plt.ylabel(r"$\mu_{{\delta}} \, \mathrm{[mas/yr]}$")

    # Zoom +/-1 sigma around the max density value.
    xmax, ymax = allr_KDE_PMs['zmax_x'], allr_KDE_PMs['zmax_y']
    zmax_x, zmax_y = allr_KDE_PMs['x'][xmax][ymax],\
        allr_KDE_PMs['y'][xmax][ymax]

    ax.contourf(PMx, PMy, PMz, cmap='Blues')
    plt.contour(
        PMx, PMy, PMz, 10, colors='k', linewidths=.5, zorder=5)
    plt.scatter(
        zmax_x, zmax_y, marker='x', c='r', s=20,
        label="[{:.2f}, {:.2f}]".format(zmax_x, zmax_y), zorder=10)

    # (bottom, left)
    xb, yl = zmax_x - xydelta, zmax_y - xydelta
    width, height = 2. * xydelta, 2. * xydelta
    rect = Rectangle(
        (xb, yl), width, height, linewidth=1, edgecolor='r',
        facecolor='none', zorder=10)
    ax.add_patch(rect)

    plt.legend()
    plt.xlim(zmax_x - xyrang, zmax_x + xyrang)
    plt.ylim(zmax_y - xyrang, zmax_y + xyrang)


def pms_coords_all(
    fig, gs, plot_style, coord, x_min, x_max, y_min, y_max, asp_ratio, x_name,
        y_name, kde_cent, clust_rad, allfr_PMs, allr_KDE_PMs, xydelta):
    """
    """
    ax = plt.subplot(gs[0:2, 4:6])
    ax.set_title(
        "In square w,h: {:.3f} [mas/yr]".format(2. * xydelta))
    ax.tick_params(axis='both', which='major')
    if plot_style == 'asteca':
        ax.grid()
    plt.xlabel('{} ({})'.format(x_name, coord))
    plt.ylabel('{} ({})'.format(y_name, coord))

    # Mask for stars inside the rectangle in the above plot.
    xmax, ymax = allr_KDE_PMs['zmax_x'], allr_KDE_PMs['zmax_y']
    zmax_x, zmax_y = allr_KDE_PMs['x'][xmax][ymax],\
        allr_KDE_PMs['y'][xmax][ymax]
    msk_rect = (allfr_PMs['pmRA'] < (zmax_x + xydelta)) &\
        (allfr_PMs['pmRA'] > (zmax_x - xydelta)) &\
        (allfr_PMs['pmDE'] < (zmax_y + xydelta)) &\
        (allfr_PMs['pmDE'] > (zmax_y - xydelta))

    pmRA, pmDE, xRA, yDE = allfr_PMs['pmRA'][msk_rect],\
        allfr_PMs['pmDE'][msk_rect], allfr_PMs['RA'][msk_rect],\
        allfr_PMs['DE'][msk_rect]
    PM_total = np.sqrt(pmRA**2 + pmDE**2)

    msk = np.argsort(PM_total)
    xRA, yDE, PM_total = xRA[msk], yDE[msk], PM_total[msk]

    PMmed, PMstd = np.median(PM_total), np.std(PM_total)
    vmin, vmax = PMmed - 3. * PMstd, PMmed + 3. * PMstd
    vmin, vmax = max(vmin, min(PM_total)), min(vmax, max(PM_total))

    cbar = plt.scatter(
        xRA, yDE, s=15, c=PM_total, vmin=vmin, vmax=vmax, alpha=.75, zorder=5)
    cbar = plt.colorbar(pad=.01, fraction=.02, aspect=40)
    cbar.ax.minorticks_off()
    cbar.set_label('Total PM')

    circle = plt.Circle(
        (kde_cent[0], kde_cent[1]), clust_rad, color='red', fill=False,
        zorder=10)
    fig.gca().add_artist(circle)

    plt.xlim(x_max, x_min)
    plt.ylim(y_min, y_max)

    ax.set_aspect(aspect=asp_ratio)


def pms_VPD_zoom(
    gs, plot_style, xlabel, coord, y_ax, clreg_PMs, fregs_PMs, raPMrng,
        dePMrng, flag_no_fl_regs):
    """
    """
    ax = plt.subplot(gs[2:4, 0:2])
    if plot_style == 'asteca':
        ax.grid()
    plt.xlabel(xlabel)
    plt.ylabel(r"$\mu_{{\delta}} \, \mathrm{[mas/yr]}$")

    pmRA, pmDE, mmag = clreg_PMs['pmRA'], clreg_PMs['pmDE'], clreg_PMs['mmag']
    msk = np.argsort(mmag)[::-1]
    pmRA, pmDE, mmag = pmRA[msk], pmDE[msk], mmag[msk]
    cmap = plt.cm.get_cmap('viridis')
    plt.scatter(
        pmRA, pmDE, marker='o', c=mmag, s=25, cmap=cmap, edgecolor='k',
        lw=.5, zorder=4, label="Cluster reg")
    cbar = plt.colorbar(pad=.01, fraction=.02, aspect=50)
    cbar.ax.minorticks_off()
    cbar.ax.invert_yaxis()
    cbar.set_label(y_ax)

    if not flag_no_fl_regs:
        pmRA, pmDE = fregs_PMs['pmRA'], fregs_PMs['pmDE']
        plt.scatter(
            pmRA, pmDE, marker='^', s=20, c='grey', linewidth=0.,
            alpha=.85, zorder=1, label="Field regs")

    plt.xlim(raPMrng)
    plt.ylim(dePMrng)
    plt.legend()


def pms_VPD_zoom_KDE(
    gs, plot_style, xlabel, coord, cr_KDE_PMs, fr_KDE_PMs, raPMrng, dePMrng,
        PMs_cent, PMs_width, PMs_height, PMs_theta, Nsigma):
    """
    KDE of PMs for cluster and field regions.
    """
    ax = plt.subplot(gs[2:4, 2:4])
    ax.set_title(
        r'Error-weighted KDE + {:.0f}$\sigma$ ellipse'.format(Nsigma))
    if plot_style == 'asteca':
        ax.grid()

    plt.xlabel(xlabel)
    plt.ylabel(r"$\mu_{{\delta}} \, \mathrm{[mas/yr]}$")

    # Cluster region contour
    CS = plt.contour(
        cr_KDE_PMs['x'], cr_KDE_PMs['y'], cr_KDE_PMs['z'], 5, colors='g',
        linewidths=1., extent=(raPMrng[0], raPMrng[1], dePMrng[0], dePMrng[1]),
        zorder=3)
    # xmax, ymax = cr_KDE_PMs['zmax_x'], cr_KDE_PMs['zmax_y']
    # zmax_x, zmax_y = cr_KDE_PMs['x'][xmax][ymax],\
    #     cr_KDE_PMs['y'][xmax][ymax]
    # CS.collections[0].set_label(
    #     r"$KDE_{{max}},\, clust: ({:.3f}, {:.3f})$".format(zmax_x, zmax_y))

    # Filed regions contours
    if bool(fr_KDE_PMs):
        CS = plt.contour(
            fr_KDE_PMs['x'], fr_KDE_PMs['y'], fr_KDE_PMs['z'], 10, colors='k',
            linewidths=.5,
            extent=(raPMrng[0], raPMrng[1], dePMrng[0], dePMrng[1]), zorder=2)
        # xmax, ymax = fr_KDE_PMs['zmax_x'], fr_KDE_PMs['zmax_y']
        # zmax_x, zmax_y = fr_KDE_PMs['x'][xmax][ymax],\
        #     fr_KDE_PMs['y'][xmax][ymax]
        # CS.collections[0].set_label(
        #     r"$KDE_{{max}},\, field: ({:.3f}, {:.3f})$".format(
        #         zmax_x, zmax_y))

    ellipse = Ellipse(
        xy=(PMs_cent[0], PMs_cent[1]), width=PMs_width,
        height=PMs_height, edgecolor='r', fc='None', lw=1.5,
        ls='--', zorder=4, angle=PMs_theta)
    ax.add_patch(ellipse)
    t1 = r"$cent=[{:.3f},{:.3f}]$".format(PMs_cent[0], PMs_cent[1])
    t2 = r"$(w,h)=[{:.3f},{:.3f}]$".format(PMs_width, PMs_height)
    lab = t1 + '\n' + t2
    plt.scatter(
        PMs_cent[0], PMs_cent[1], c='r', marker='x', s=25, zorder=4,
        label=lab)

    plt.xlim(*raPMrng)
    plt.ylim(*dePMrng)
    plt.legend()


def pms_VPD_zoom_MP(
    gs, plot_style, xlabel, coord, clreg_PMs, fregs_PMs, raPMrng,
        dePMrng):
    """
    """
    ax = plt.subplot(gs[2:4, 4:6])
    if plot_style == 'asteca':
        ax.grid()

    plt.xlabel(xlabel)
    plt.ylabel(r"$\mu_{{\delta}} \, \mathrm{[mas/yr]}$")

    cmap = plt.cm.get_cmap('RdYlBu_r')
    pmMP = clreg_PMs['MP']
    norm = Normalize(vmin=pmMP.min(), vmax=pmMP.max())

    ax.errorbar(
        clreg_PMs['pmRA'], clreg_PMs['pmDE'], yerr=clreg_PMs['epmDE'],
        xerr=clreg_PMs['epmRA'], fmt='none',
        elinewidth=.65, ecolor=cmap(norm(pmMP)), zorder=4)
    plt.scatter(clreg_PMs['pmRA'], clreg_PMs['pmDE'], c=pmMP, cmap=cmap, s=0.)
    cbar = plt.colorbar(pad=.01, fraction=.02, aspect=50)
    cbar.ax.minorticks_off()
    cbar.set_label("MPs")

    ax.errorbar(
        fregs_PMs['pmRA'], fregs_PMs['pmDE'], yerr=fregs_PMs['epmDE'],
        xerr=fregs_PMs['epmRA'], fmt='none',
        elinewidth=.35, ecolor='grey', zorder=1)

    plt.xlim(raPMrng)
    plt.ylim(dePMrng)


# DEPRECATED 04/2021
# def pms_vs_mag(
#     gs, plot_style, xlabel, coord, y_ax, clreg_PMs, fregs_PMs, raPMrng,
#         dePMrng):
#     """
#     """
#     pmMP = clreg_PMs['MP']
#     ymin, ymax = max(clreg_PMs['mmag']) + .5, min(clreg_PMs['mmag']) + .5

#     def axplot(ax, x_range, xlabel, cl_x, cl_y, fr_x, fr_y):
#         if plot_style == 'asteca':
#             ax.grid()

#         ax.set_xlabel(xlabel)
#         ax.set_ylabel(y_ax)
#         cmap = plt.cm.get_cmap('RdYlBu_r')

#         ax.scatter(
#             cl_x, cl_y, marker='o', c=pmMP, s=30, edgecolors='black',
#             cmap=cmap, lw=0.35, zorder=5)
#         ax.scatter(fr_x, fr_y, c='grey', s=15, marker='^', zorder=2)
#         ax.set_xlim(x_range)
#         ax.set_ylim(ymin, ymax)

#     ax1 = plt.subplot(gs[4:6, 0:2])
#     ax2 = plt.subplot(gs[4:6, 2:4])

#     axplot(
#         ax1, raPMrng, xlabel, clreg_PMs['pmRA'], clreg_PMs['mmag'],
#         fregs_PMs['pmRA'], fregs_PMs['mmag'])
#     xlabel = r"$\mu_{{\delta}} \, \mathrm{[mas/yr]}$"
#     axplot(ax2, dePMrng, xlabel, clreg_PMs['pmDE'], clreg_PMs['mmag'],
#            fregs_PMs['pmDE'], fregs_PMs['mmag'])
#
# def pms_dist(gs, plot_style, y_ax, clreg_PMs, pm_dist_max):
#     """
#     """
#     ax = plt.subplot(gs[4:6, 4:6])
#     if plot_style == 'asteca':
#         ax.grid()

#     ax.set_title("Distance to KDE's cluster max")
#     plt.xlabel("PM dist [mas/yr]")
#     plt.ylabel(y_ax)

#     cmap = plt.cm.get_cmap('RdYlBu_r')
#     plt.scatter(
#         pm_dist_max, clreg_PMs['mmag'], marker='o', c=clreg_PMs['MP'],
#         s=30, edgecolors='black', cmap=cmap, lw=0.35, zorder=4)

#     plx_mean, plx_median, plx_std = sigma_clipped_stats(pm_dist_max)
#     plt.xlim(-.05, min(max(pm_dist_max), plx_median + 4. * plx_std))
#     plt.gca().invert_yaxis()


def plot(N, *args):
    """
    Handle each plot separately.
    """

    plt_map = {
        0: [pms_VPD_all, 'VPD for all stars'],
        1: [pms_VPD_KDE_all, 'VPD KDE for all stars'],
        2: [pms_coords_all, 'x/RA vs y/DEC PMs colored'],
        #
        3: [pms_VPD_zoom, 'VPD cluster region'],
        4: [pms_VPD_zoom_KDE, 'PMs KDE diagram'],
        5: [pms_VPD_zoom_MP, 'PMs VPD MP colored']
        # 6: [pms_vs_mag, 'PMs mag MP colored'],
        # 7: [pms_dist, 'PMs distance to center']
    }

    fxn = plt_map.get(N, None)[0]
    if fxn is None:
        raise ValueError("  ERROR: there is no plot {}.".format(N))

    try:
        fxn(*args)
    except Exception:
        import traceback
        print(traceback.format_exc())
        print("  WARNING: error when plotting {}".format(plt_map.get(N)[1]))


# DEPRECATED 18/01/20
# def pms_NN_all(
#     gs, coord, y_ax, pmRA_all, pmDE_all, pmMag_all, PMs_d_median, nnmax,
#         nnperc):
#     """
#     N-N of PMs for all the frame.
#     """
#     ax = plt.subplot(gs[0:2, 2:4])
#     ax.set_title(
#         'All stars in frame, d<{:.1f}% prcntl'.format(nnperc), fontsize=9)
#     ax.grid(b=True, which='major', color='gray', linestyle='--', lw=.5,
#             zorder=-1)
#     if coord == 'deg':
#         plt.xlabel(
#             r"$\mu_{{\alpha}} \, cos \delta \, \mathrm{[mas/yr]}$",
#             fontsize=12)
#     else:
#         plt.xlabel(r"$\mu_{{\alpha}} \, \mathrm{[mas/yr]}$", fontsize=12)
#     plt.ylabel(r"$\mu_{{\delta}} \, \mathrm{[mas/yr]}$", fontsize=12)

#     dmax = np.percentile(PMs_d_median, nnperc)
#     msk = PMs_d_median < dmax
#     pmRA, pmDE, mag = pmRA_all[msk], pmDE_all[msk], pmMag_all[msk]

#     cmap = plt.cm.get_cmap('viridis')
#     msk = np.argsort(mag)[::-1]
#     pmRA, pmDE, mag = pmRA[msk], pmDE[msk], mag[msk]
#     plt.scatter(
#         pmRA, pmDE, s=25, c=mag, cmap=cmap, edgecolor='k', lw=.5,
#         label="N-N={}\nN={}".format(nnmax, pmRA.size), zorder=2)
#     cbar = plt.colorbar(pad=.01, fraction=.02, aspect=50)
#     cbar.ax.invert_yaxis()
#     cbar.set_label(y_ax, size=8)

#     plt.legend(fontsize='small')


# DEPRECATED 18/01/20
# def pms_coords_all(
#     fig, gs, cld_i, y_ax, x_name, y_name, coord, PMs_d_median, nnperc,
#         xRA_all, yDE_all, pmMag_all, kde_cent, clust_rad):
#     """
#     """
#     ax = plt.subplot(gs[0:2, 4:6])
#     ax.set_title('All stars in frame, filtered by N-N', fontsize=9)
#     ax.grid(b=True, which='major', color='gray', linestyle='--', lw=.5,
#             zorder=-1)
#     plt.xlabel('{} ({})'.format(x_name, coord), fontsize=12)
#     plt.ylabel('{} ({})'.format(y_name, coord), fontsize=12)

#     dmax = np.percentile(PMs_d_median, nnperc)
#     msk = PMs_d_median < dmax
#     xRA, yDE = xRA_all[msk], yDE_all[msk]

#     circle = plt.Circle(
#         (kde_cent[0], kde_cent[1]), clust_rad, color='green', fill=False)
#     fig.gca().add_artist(circle)

#     x_min, x_max, y_min, y_max = prep_plots.frame_max_min(
#         cld_i['x'], cld_i['y'])
#     asp_ratio = prep_plots.aspect_ratio(x_min, x_max, y_min, y_max)

#     N_in_r = sum(
#         np.sqrt((x - kde_cent[0])**2 + (y - kde_cent[1])**2) < clust_rad
#         for x, y in np.array([xRA, yDE]).T)

#     st_sizes_arr = prep_plots.star_size(pmMag_all)
#     # zmin, zmax = st_sizes_arr.min(), st_sizes_arr.max()
#     # st_sizes_arr = prep_plots.star_size(pmMag_all, zmin=zmin, zmax=zmax)

#     plt.scatter(
#         xRA, yDE, c='k', s=st_sizes_arr[msk],
#         label=r"$N_{{r<r_{{cl}}}}$={}".format(N_in_r), zorder=2)
#     # mag, cmap=cmap, edgecolor='k', lw=.5)
#     # cbar = plt.colorbar(pad=.01, fraction=.02, aspect=50)
#     # cbar.ax.invert_yaxis()
#     # cbar.set_label(y_ax, size=8)
#     plt.legend(fontsize='small')

#     ax.set_aspect(aspect=asp_ratio)
#     plt.xlim(x_min, x_max)
#     plt.ylim(y_min, y_max)

#     if coord == 'deg':
#         ax.invert_xaxis()

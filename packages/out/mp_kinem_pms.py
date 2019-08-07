
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Ellipse
import numpy as np
from astropy.stats import sigma_clipped_stats
from . import prep_plots


def pms_KDE_all(
    gs, coord, y_ax, pmRA_all, pmDE_all, pmMag_all, PM_kde_all,
        PM_KDE_std):
    """
    KDE of PMs for all the frame.
    """
    ax = plt.subplot(gs[0:2, 0:2])
    ax.set_title(
        r'All stars in frame + KDE ($\mu\pm\sigma={:.2f}$)'.format(PM_KDE_std),
        fontsize=9)
    ax.minorticks_on()
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=.5,
            zorder=-1)
    if coord == 'deg':
        plt.xlabel(
            r"$\mu_{{\alpha}} \, cos \delta \, \mathrm{[mas/yr]}$",
            fontsize=12)
    else:
        plt.xlabel(r"$\mu_{{\alpha}} \, \mathrm{[mas/yr]}$", fontsize=12)
    plt.ylabel(r"$\mu_{{\delta}} \, \mathrm{[mas/yr]}$", fontsize=12)

    PMx, PMy, PMz, extent = PM_kde_all

    # plt.imshow(np.rot90(PMz), cmap='Blues', zorder=1, extent=extent)

    cmap = plt.cm.get_cmap('viridis')
    msk = np.argsort(pmMag_all)[::-1]
    pmRA_all, pmDE_all, mag = pmRA_all[msk], pmDE_all[msk], pmMag_all[msk]
    plt.scatter(
        pmRA_all, pmDE_all, s=25, c=mag, cmap=cmap,
        edgecolor='w', lw=.5, zorder=3)
    cbar = plt.colorbar(pad=.01, fraction=.02, aspect=50)
    cbar.ax.tick_params(labelsize=7)
    cbar.ax.invert_yaxis()
    cbar.set_label(y_ax, size=8)

    # Cluster region contour
    plt.contour(
        PMx, PMy, PMz, 10, colors='k', linewidths=.75, zorder=5,
        extent=extent)

    plt.xlim(extent[0], extent[1])
    plt.ylim(extent[2], extent[3])


def pms_NN_all(
    gs, coord, y_ax, pmRA_all, pmDE_all, pmMag_all, PMs_d_median, nnmax,
        nnperc):
    """
    N-N of PMs for all the frame.
    """
    ax = plt.subplot(gs[0:2, 2:4])
    ax.set_title(
        'All stars in frame, d<{:.1f}% prcntl'.format(nnperc), fontsize=9)
    ax.minorticks_on()
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=.5,
            zorder=-1)
    if coord == 'deg':
        plt.xlabel(
            r"$\mu_{{\alpha}} \, cos \delta \, \mathrm{[mas/yr]}$",
            fontsize=12)
    else:
        plt.xlabel(r"$\mu_{{\alpha}} \, \mathrm{[mas/yr]}$", fontsize=12)
    plt.ylabel(r"$\mu_{{\delta}} \, \mathrm{[mas/yr]}$", fontsize=12)

    dmax = np.percentile(PMs_d_median, nnperc)
    msk = PMs_d_median < dmax
    pmRA, pmDE, mag = pmRA_all[msk], pmDE_all[msk], pmMag_all[msk]

    cmap = plt.cm.get_cmap('viridis')
    msk = np.argsort(mag)[::-1]
    pmRA, pmDE, mag = pmRA[msk], pmDE[msk], mag[msk]
    plt.scatter(
        pmRA, pmDE, s=25, c=mag, cmap=cmap, edgecolor='k', lw=.5,
        label="N-N={}\nN={}".format(nnmax, pmRA.size), zorder=2)
    cbar = plt.colorbar(pad=.01, fraction=.02, aspect=50)
    cbar.ax.tick_params(labelsize=7)
    cbar.ax.invert_yaxis()
    cbar.set_label(y_ax, size=8)

    plt.legend(fontsize='small')


def pms_coords_all(
    fig, gs, cld_i, y_ax, x_name, y_name, coord, PMs_d_median, nnperc,
        xRA_all, yDE_all, pmMag_all, kde_cent, clust_rad):
    """
    """
    ax = plt.subplot(gs[0:2, 4:6])
    ax.set_title('All stars in frame, filtered by N-N', fontsize=9)
    ax.minorticks_on()
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=.5,
            zorder=-1)
    plt.xlabel('{} ({})'.format(x_name, coord), fontsize=12)
    plt.ylabel('{} ({})'.format(y_name, coord), fontsize=12)

    dmax = np.percentile(PMs_d_median, nnperc)
    msk = PMs_d_median < dmax
    xRA, yDE = xRA_all[msk], yDE_all[msk]

    circle = plt.Circle(
        (kde_cent[0], kde_cent[1]), clust_rad, color='green', fill=False)
    fig.gca().add_artist(circle)

    x_min, x_max, y_min, y_max = prep_plots.frame_max_min(
        cld_i['x'], cld_i['y'])
    asp_ratio = prep_plots.aspect_ratio(x_min, x_max, y_min, y_max)

    N_in_r = sum(
        np.sqrt((x - kde_cent[0])**2 + (y - kde_cent[1])**2) < clust_rad
        for x, y in np.array([xRA, yDE]).T)

    st_sizes_arr = prep_plots.star_size(pmMag_all)
    # zmin, zmax = st_sizes_arr.min(), st_sizes_arr.max()
    # st_sizes_arr = prep_plots.star_size(pmMag_all, zmin=zmin, zmax=zmax)

    plt.scatter(
        xRA, yDE, c='k', s=st_sizes_arr[msk],
        label=r"$N_{{r<r_{{cl}}}}$={}".format(N_in_r), zorder=2)
    # mag, cmap=cmap, edgecolor='k', lw=.5)
    # cbar = plt.colorbar(pad=.01, fraction=.02, aspect=50)
    # cbar.ax.tick_params(labelsize=7)
    # cbar.ax.invert_yaxis()
    # cbar.set_label(y_ax, size=8)
    plt.legend(fontsize='small')

    ax.set_aspect(aspect=asp_ratio)
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)

    if coord == 'deg':
        ax.invert_xaxis()


def pms_vpd_mag(
    gs, coord, y_ax, pmRA_DE, pmDE, mmag_pm, pmRA_fl_DE, pmDE_fl,
        pm_mag_fl, raPMrng, dePMrng, flag_no_fl_regs_i):
    '''
    '''
    ax = plt.subplot(gs[2:4, 0:2])
    ax.minorticks_on()
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=.5,
            zorder=1)

    if coord == 'deg':
        plt.xlabel(
            r"$\mu_{{\alpha}} \, cos \delta \, \mathrm{[mas/yr]}$",
            fontsize=12)
    else:
        plt.xlabel(r"$\mu_{{\alpha}} \, \mathrm{[mas/yr]}$", fontsize=12)
    plt.ylabel(r"$\mu_{{\delta}} \, \mathrm{[mas/yr]}$", fontsize=12)

    msk = np.argsort(mmag_pm)[::-1]
    pmRA_DE, pmDE, mmag_pm = pmRA_DE[msk], pmDE[msk], mmag_pm[msk]
    cmap = plt.cm.get_cmap('viridis')
    plt.scatter(
        pmRA_DE, pmDE, marker='o', c=mmag_pm, s=25, cmap=cmap, edgecolor='k',
        lw=.5, zorder=4, label="Cluster reg")
    cbar = plt.colorbar(pad=.01, fraction=.02, aspect=50)
    cbar.ax.tick_params(labelsize=7)
    cbar.ax.invert_yaxis()
    cbar.set_label(y_ax, size=8)

    if not flag_no_fl_regs_i:
        msk = np.argsort(pm_mag_fl)[::-1]
        pmRA_fl_DE, pmDE_fl, pm_mag_fl = pmRA_fl_DE[msk], pmDE_fl[msk],\
            pm_mag_fl[msk]
        plt.scatter(
            pmRA_fl_DE, pmDE_fl, marker='^', s=20, c='grey', linewidth=0.,
            cmap=cmap, alpha=.85, zorder=1, label="Field regs")

    plt.xlim(raPMrng)
    plt.ylim(dePMrng)
    plt.legend(fontsize='small')


def pms_KDE_diag(
    gs, coord, PM_cl_x, PM_cl_y, PM_cl_z, PM_fl_x, PM_fl_y, PM_fl_z,
    PMs_cl_cx, PMs_cl_cy, PMs_fl_cx, PMs_fl_cy, raPMrng, dePMrng, PMs_cent,
        PMs_width, PMs_height, PMs_theta, CI_prob):
    """
    KDE of PMs for cluster and field regions.
    """
    ax = plt.subplot(gs[2:4, 2:4])
    ax.set_title(
        'Error-weighted KDE + {:.0f}% CI ellipse'.format(CI_prob * 100.),
        fontsize=9)
    ax.minorticks_on()
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=.5,
            zorder=1)

    if coord == 'deg':
        plt.xlabel(
            r"$\mu_{{\alpha}} \, cos \delta \, \mathrm{[mas/yr]}$",
            fontsize=12)
    else:
        plt.xlabel(r"$\mu_{{\alpha}} \, \mathrm{[mas/yr]}$", fontsize=12)
    plt.ylabel(r"$\mu_{{\delta}} \, \mathrm{[mas/yr]}$", fontsize=12)

    # Cluster region contour
    CS = plt.contour(
        PM_cl_x, PM_cl_y, PM_cl_z, 5, colors='g', linewidths=1.,
        extent=(raPMrng[0], raPMrng[1], dePMrng[0], dePMrng[1]),
        zorder=3)
    CS.collections[0].set_label(
        r"$KDE_{{max}},\, clust: ({:.3f}, {:.3f})$".format(
            PM_cl_x[PMs_cl_cx][PMs_cl_cy], PM_cl_y[PMs_cl_cx][PMs_cl_cy]))

    # Filed regions contours
    if PM_fl_z.any():
        CS = plt.contour(
            PM_fl_x, PM_fl_y, PM_fl_z, 10, colors='k', linewidths=.5,
            extent=(raPMrng[0], raPMrng[1], dePMrng[0], dePMrng[1]),
            zorder=2)
        CS.collections[0].set_label(
            r"$KDE_{{max}},\, field: ({:.3f}, {:.3f})$".format(
                PM_fl_x[PMs_fl_cx][PMs_fl_cy], PM_fl_y[PMs_fl_cx][PMs_fl_cy]))

    ellipse = Ellipse(
        xy=(PMs_cent[0], PMs_cent[1]), width=PMs_width,
        height=PMs_height, edgecolor='r', fc='None', lw=1.5,
        ls='--', zorder=4, angle=PMs_theta)
    ax.add_patch(ellipse)
    plt.scatter(
        PMs_cent[0], PMs_cent[1], c='r', marker='x', s=25, zorder=4,
        label=(r"$cent=[{:.3f},{:.3f}]$" + '\n' +
               r"$(w,h)=[{:.3f},{:.3f}]$").format(
                   PMs_cent[0], PMs_cent[1], PMs_width, PMs_height))

    plt.xlim(*raPMrng)
    plt.ylim(*dePMrng)
    plt.legend(fontsize='small')


def pms_vpd_mp(
    gs, coord, pmMP, pmRA_DE, e_pmRA_DE, pmDE, e_pmDE, pmRA_fl_DE,
        e_pmRA_fl_DE, pmDE_fl, e_pmDE_fl, raPMrng, dePMrng):
    '''
    '''
    ax = plt.subplot(gs[2:4, 4:6])
    ax.minorticks_on()
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=.5,
            zorder=1)

    if coord == 'deg':
        plt.xlabel(
            r"$\mu_{{\alpha}} \, cos \delta \, \mathrm{[mas/yr]}$",
            fontsize=12)
    else:
        plt.xlabel(r"$\mu_{{\alpha}} \, \mathrm{[mas/yr]}$", fontsize=12)
    plt.ylabel(r"$\mu_{{\delta}} \, \mathrm{[mas/yr]}$", fontsize=12)

    cmap = plt.cm.get_cmap('RdYlBu_r')
    norm = Normalize(vmin=pmMP.min(), vmax=pmMP.max())

    ax.errorbar(
        pmRA_DE, pmDE, yerr=e_pmDE, xerr=e_pmRA_DE, fmt='none',
        elinewidth=.65, ecolor=cmap(norm(pmMP)), zorder=4)
    plt.scatter(pmRA_DE, pmDE, c=pmMP, cmap=cmap, s=0.)
    cbar = plt.colorbar(pad=.01, fraction=.02, aspect=50)
    cbar.ax.tick_params(labelsize=7)
    cbar.set_label("MPs", size=8)

    ax.errorbar(
        pmRA_fl_DE, pmDE_fl, yerr=e_pmDE_fl, xerr=e_pmRA_fl_DE, fmt='none',
        elinewidth=.35, ecolor='grey', zorder=1)

    plt.xlim(raPMrng)
    plt.ylim(dePMrng)


def pms_vs_mag(
    gs, coord, y_ax, pmMP, pmRA_DE, pmDE, mmag_pm, pmRA_fl_DE,
        pmDE_fl, pm_mag_fl, raPMrng, dePMrng):
    '''
    '''
    def axplot(ax, x_range, xlabel, cl_x, cl_y, fr_x, fr_y, ymax, ymin):
        ax.minorticks_on()
        ax.grid(b=True, which='major', color='gray', linestyle='--', lw=.5,
                zorder=1)

        ax.set_xlabel(xlabel, fontsize=12)
        ax.set_ylabel(y_ax, fontsize=12)
        cmap = plt.cm.get_cmap('RdYlBu_r')

        ax.scatter(
            cl_x, cl_y, marker='o', c=pmMP, s=30, edgecolors='black',
            cmap=cmap, lw=0.35, zorder=5)
        ax.scatter(fr_x, fr_y, c='grey', s=15, marker='^', zorder=2)
        ax.set_xlim(x_range)
        ax.set_ylim(ax.get_ylim()[::-1])
        ax.set_ylim(ymin, ymax)

    ax1 = plt.subplot(gs[4:6, 0:2])
    ax2 = plt.subplot(gs[4:6, 2:4])

    if coord == 'deg':
        xlabel = r"$\mu_{{\alpha}} \, cos \delta \, \mathrm{[mas/yr]}$"
    else:
        xlabel = r"$\mu_{{\alpha}} \, \mathrm{[mas/yr]}$"
    ymax, ymin = min(mmag_pm) - .5, max(mmag_pm) + .5
    axplot(
        ax1, raPMrng, xlabel, pmRA_DE, mmag_pm, pmRA_fl_DE, pm_mag_fl, ymax,
        ymin)
    xlabel = r"$\mu_{{\delta}} \, \mathrm{[mas/yr]}$"
    axplot(ax2, dePMrng, xlabel, pmDE, mmag_pm, pmDE_fl, pm_mag_fl, ymax, ymin)


def pms_dist(gs, y_ax, pmMP, pm_dist_max, mmag_pm):
    '''
    '''
    ax = plt.subplot(gs[4:6, 4:6])
    ax.minorticks_on()
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=.5,
            zorder=1)

    ax.set_title("Distance to KDE's cluster max", fontsize=9)
    plt.xlabel("PM dist [mas/yr]", fontsize=12)
    plt.ylabel(y_ax, fontsize=12)

    cmap = plt.cm.get_cmap('RdYlBu_r')
    plt.scatter(
        pm_dist_max, mmag_pm, marker='o', c=pmMP, s=30, edgecolors='black',
        cmap=cmap, lw=0.35, zorder=4)

    plx_mean, plx_median, plx_std = sigma_clipped_stats(pm_dist_max)
    plt.xlim(-.05, plx_median + 2. * plx_std)
    plt.gca().invert_yaxis()


def plot(N, *args):
    '''
    Handle each plot separately.
    '''

    plt_map = {
        0: [pms_KDE_all, 'PMs VPD KDE all'],
        1: [pms_NN_all, 'PMs VPD NN all'],
        2: [pms_coords_all, 'x/RA vs y/DEC NN all'],
        3: [pms_vpd_mag, 'PMs VPD mag colored'],
        4: [pms_KDE_diag, 'PMs KDE diagram'],
        5: [pms_vpd_mp, 'PMs VPD MP colored'],
        6: [pms_vs_mag, 'PMs mag MP colored'],
        7: [pms_dist, 'PMs distance to center']
    }

    fxn = plt_map.get(N, None)[0]
    if fxn is None:
        raise ValueError("  ERROR: there is no plot {}.".format(N))

    try:
        fxn(*args)
    except Exception:
        import traceback
        print(traceback.format_exc())
        print("  WARNING: error when plotting {}.".format(plt_map.get(N)[1]))

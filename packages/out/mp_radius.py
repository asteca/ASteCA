
import numpy as np
from scipy.stats import median_abs_deviation as MAD
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib import patches as mpatches
from itertools import cycle
from ..structure.king_profile import KingProf as kpf
from ..structure.king_profile import KP_memb_x


def pl_rad_find(
    gs, plot_style, coord, clust_rad, membs_ratio, membs_ratio_smooth, CI_vals,
        rad_radii):
    """
    Radius estimation plot.
    """
    # Convert from deg to arcmin
    rad_radii, clust_rad = rad_radii * 60., clust_rad * 60.
    coord2 = 'arcmin'

    # # Cut when (if) N_nembs = 0
    # try:
    #     icut = np.where(rad_N_membs == 0.)[0][0]
    #     if icut == 0:
    #         icut = len(rad_N_membs)
    #     else:
    #         # Use 30 points since the number of radii values in
    #         # radius.rdpAreasDists() is hardcoded to 300.
    #         icut = icut + 30
    # except IndexError:
    #     icut = len(rad_N_membs)
    # rad_rads, rad_N_membs, rad_N_field, rad_CI = rad_rads[:icut],\
    #     rad_N_membs[:icut], rad_N_field[:icut], rad_CI[:icut]

    ax = plt.subplot(gs[0:2, 0:2])
    plt.xlabel(r'radius $[{}]$'.format(coord2))
    # plt.ylabel(r"$N_{memb}\;|\;N_{field}\;|\;CI$"))
    if plot_style == 'asteca':
        ax.grid()

    msk = membs_ratio > 0
    plt.scatter(rad_radii[msk], membs_ratio[msk], c='g',
                label=r'$N_{memb}/N_{ring}$', alpha=.3)
    if membs_ratio_smooth.any():
        if membs_ratio_smooth[msk].any():
            plt.plot(rad_radii[msk], membs_ratio_smooth[msk], c='b')

    plt.plot(rad_radii[msk], CI_vals[msk], ls=':', c='k', label='CI')
    ymin, _ = ax.get_ylim()
    ymin = max(ymin, 0)
    plt.axvline(x=clust_rad, lw=1.5, color='r', label=r"$r_{cl}$")

    # Legends.
    leg = plt.legend(fancybox=True)
    leg.get_frame().set_alpha(0.7)

    # xmin, xmax = ax.get_xlim()
    xmax = min(clust_rad + clust_rad * 3, rad_radii[msk][-1])
    plt.ylim(ymin, 1.02)
    plt.xlim(0., xmax)


def pl_mag_membs(gs, plot_style, y_ax, membvsmag):
    """
    """
    ax = plt.subplot(gs[0:2, 2:4])
    if plot_style == 'asteca':
        ax.grid()
    ax.set_title(r"$N_{{memb}}$ vs magnitude cut (phot incomp); $r_{{cl}}$")
    plt.xlabel('$' + y_ax + '$')
    plt.ylabel(r'$N_{memb}$')
    if membvsmag.any():
        plt.bar(*membvsmag, zorder=5)


def pl_cl_fl_regions(
    gs, fig, plot_style, x_name, y_name, coord, x_min, x_max, y_min, y_max,
    asp_ratio, kde_cent, clust_rad, field_regions_i, field_regions_rjct_i,
        cl_region_i, cl_region_rjct_i, flag_no_fl_regs_i):
    """
    Cluster and field regions defined.
    """
    ax = plt.subplot(gs[0:2, 4:6])

    ax.set_aspect(aspect=asp_ratio)
    ax.invert_xaxis()
    plt.xlabel('{} ({})'.format(x_name, coord))
    plt.ylabel('{} ({})'.format(y_name, coord))
    if plot_style == 'asteca':
        ax.grid(which='both')
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
        len(field_regions_i)))


def pl_rad_dens(
    gs, plot_style, coord, rdp_radii, rdp_points, rdp_stddev, rad_max,
    field_dens, e_fdens, clust_rad, kp_ndim, KP_Bys_rc, KP_Bys_rt,
        KP_plot, KP_conct_par):
    """
    Radial density plot.
    """

    KP_cent_dens, _16_84_rang, _84_kp, _16_kp = 0., 0., 0., 0.
    if kp_ndim in (2, 4):
        KP_cent_dens, _16_84_rang, _84_kp, _16_kp = KP_plot['KP_cent_dens'],\
            KP_plot['_16_84_rang'], KP_plot['_84_kp'], KP_plot['_16_kp']

    # Convert from deg to arcmin
    rdp_radii = np.array(rdp_radii) * 60.
    clust_rad, rad_max = clust_rad * 60., rad_max * 60.
    KP_Bys_rc, KP_Bys_rt = KP_Bys_rc * 60., KP_Bys_rt * 60.
    field_dens, e_fdens, KP_cent_dens = field_dens / 3600.,\
        e_fdens / 3600., KP_cent_dens / 3600.
    _16_84_rang, _84_kp, _16_kp = _16_84_rang * 60.,\
        _84_kp / 3600., _16_kp / 3600.
    rdp_points = np.array(rdp_points) / 3600.
    rdp_stddev = np.array(rdp_stddev) / 3600.
    coord2 = 'arcmin'

    ax = plt.subplot(gs[2:4, 0:4])

    plt.xlabel(r'radius $[{}]$'.format(coord2))
    plt.ylabel(r"$\rho$ $[st/{}^{{2}}]$".format(coord2))
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    if plot_style == 'asteca':
        ax.grid()

    # Legend texts
    r_frmt = '{:.0f}' if coord2 == 'px' else '{:.2f}'
    t_rad = r"$r_{{{}}}=" + r_frmt + r"\,[{}]$"
    kp_rad = r"$r_{{{}}}=" + r_frmt + r"_{{" + r_frmt + r"}}^{{" +\
        r_frmt + r"}}\,[{}]$"

    # Plot density profile
    ax.plot(rdp_radii, rdp_points, marker='o', ms=5, lw=0., zorder=3,
            label="RDP")
    # Plot error bars
    plt.errorbar(
        rdp_radii, rdp_points, yerr=rdp_stddev, fmt='none', ecolor='grey',
        lw=1., zorder=1)

    # Plot background level.
    ax.hlines(y=field_dens, xmin=0, xmax=max(rdp_radii), color='k', ls='--',
              zorder=5)
    if not np.isnan(e_fdens):
        ax.hlines(
            y=field_dens - e_fdens, xmin=0, xmax=max(rdp_radii),
            color='k', ls=':', zorder=5)
        ax.hlines(
            y=field_dens + e_fdens, xmin=0, xmax=max(rdp_radii),
            color='k', ls=':', zorder=5)

    # Set plot limits
    delta_backg = 0.15 * (max(rdp_points) - field_dens)
    y_min = field_dens - delta_backg
    y_min = 0. if y_min <= 0. else y_min
    N_half = int(len(rdp_points) * .5)
    y_max = min(ax.get_ylim()[1], max(rdp_points[:N_half]) + delta_backg)
    plt.xlim(0, rad_max)
    plt.ylim(y_min, y_max)
    y_mid_point = np.mean(ax.get_ylim())

    # Plot radius.
    ax.vlines(x=clust_rad, ymin=field_dens, ymax=y_mid_point, lw=1.5,
              color='r', label=t_rad.format("cl", clust_rad, coord2), zorder=5)

    # Plot King profile. Use median values
    if kp_ndim in (2, 4):
        txts = [
            'King prof ({:.2f})'.format(KP_conct_par),
            kp_rad.format(
                "c", KP_Bys_rc[1], KP_Bys_rc[0], KP_Bys_rc[2], coord2),
            kp_rad.format(
                "t", KP_Bys_rt[1], KP_Bys_rt[0], KP_Bys_rt[2], coord2)
        ]
        # Plot curve. Values outside of rt contribute 'fd'.
        kpf_xvals = np.linspace(rdp_radii[0], KP_Bys_rt[1], 100)
        kpf_yvals = KP_cent_dens * kpf(
            kpf_xvals, KP_Bys_rc[1], KP_Bys_rt[1]) + field_dens
        ax.plot(kpf_xvals, kpf_yvals, 'g--', label=txts[0], lw=2., zorder=3)
        # 16-84th range
        idx = (np.abs(_16_84_rang - kpf_xvals[-1])).argmin()
        # 16-84 region
        ax.fill_between(
            _16_84_rang[:idx], _84_kp[:idx], _16_kp[:idx], facecolor='green',
            alpha=0.1)

        # Core radius
        rc_ymax = KP_cent_dens * kpf(
            KP_Bys_rc[1], KP_Bys_rc[1], KP_Bys_rt[1]) + field_dens
        ax.vlines(
            x=KP_Bys_rc[1], ymin=field_dens, ymax=rc_ymax, label=txts[1],
            color='g', linestyles=':', lw=2., zorder=5)
        # Tidal radius
        ax.vlines(x=KP_Bys_rt[1], ymin=field_dens, ymax=y_mid_point,
                  label=txts[2], color='g')

    # get handles
    handles, labels = ax.get_legend_handles_labels()
    # use them in the legend
    ax.legend(handles, labels, loc='upper center', numpoints=2)

    #
    # Log-log plot
    axins = inset_axes(ax, width=2.5, height=2.5)
    if plot_style == 'asteca':
        axins.grid(which='both')

    # plt.xlabel(r'log(radius) $[{}]$'.format(coord2), fontsize=8)
    # plt.ylabel(r"log($\rho$) $[st/{}^{{2}}]$".format(coord2), fontsize=8)

    axins.scatter(rdp_radii, rdp_points, s=5, zorder=5)
    axins.hlines(y=field_dens, xmin=min(rdp_radii), xmax=max(rdp_radii),
                 color='k', ls='--', zorder=5)

    # Plot King profile.
    if kp_ndim in (2, 4):
        axins.plot(kpf_xvals, kpf_yvals, 'g--', lw=1., zorder=3)
        axins.vlines(
            x=KP_Bys_rc[3], ymin=field_dens, ymax=rc_ymax, color='g',
            linestyles=':', lw=1.)
        axins.vlines(
            x=KP_Bys_rt[3], ymin=field_dens, ymax=y_mid_point, color='g', lw=1)
    else:
        axins.vlines(
            x=clust_rad, ymin=field_dens, ymax=y_mid_point, lw=1.5, color='r',
            zorder=10)

    axins.set_xscale('log')
    axins.set_yscale('log')
    axins.minorticks_off()
    # y0, y1 = axins.set_ylim()
    # axins.set_ylim(y0, y_max)
    if y_min > 1:
        axins.set_ylim(y_min, y_max)


def pl_zoom_frame(
    gs, fig, x_name, y_name, coord, x_zmin, x_zmax, y_zmin, y_zmax, cont_index,
    x_data, y_data, st_sizes_arr, kde_cent, clust_rad, KP_Bys_rc,
        KP_Bys_rt, KP_Bys_ecc, KP_Bys_theta, frac_cl_area, kp_ndim):
    """
    Zoom on x,y finding chart.
    """
    ax = plt.subplot(gs[2:4, 4:6])

    # Force square plot.
    ax.set_aspect('equal')
    # Set plot limits.
    plt.xlim(x_zmin, x_zmax)
    plt.ylim(y_zmin, y_zmax)
    ax.invert_xaxis()
    # Set axis labels
    plt.xlabel('{} ({})'.format(x_name, coord))
    plt.ylabel('{} ({})'.format(y_name, coord))

    # Plot radius.
    circle = plt.Circle((kde_cent[0], kde_cent[1]), clust_rad, color='r',
                        fill=False, lw=1.5, zorder=5)
    # Stars inside the cluster
    N_in = 0
    for x, y in zip(*[x_data, y_data]):
        if np.sqrt((x - kde_cent[0])**2 + (y - kde_cent[1])**2) <= clust_rad:
            N_in += 1

    ax.set_title(
        (r"$A_{{fr}}$={:.0f}%, $N_{{r<r_{{cl}}}}$={}, "
            "[$CI = {:.2f}$] (phot incomp)").format(
            100. * frac_cl_area, N_in, cont_index))
    fig.gca().add_artist(circle)

    # Core and tidal radii
    if kp_ndim in (2, 4):
        # Plot tidal radius.
        rt_mode, ecc_mode, theta_mode = KP_Bys_rt[4], KP_Bys_ecc[4],\
            KP_Bys_theta[4]
        b = rt_mode * np.sqrt(1. - ecc_mode**2)
        ellipse = mpatches.Ellipse(
            xy=kde_cent, width=2. * rt_mode, height=2. * b,
            angle=np.rad2deg(theta_mode), facecolor='None', edgecolor='g',
            linewidth=1.5, transform=ax.transData)
        # circle = plt.Circle(
        #     (kde_cent[0], kde_cent[1]), KP_Bys_rt[3], color='g', fill=False,
        #     lw=1.5)
        fig.gca().add_artist(ellipse)
        # Plot core radius.
        rc_mode = KP_Bys_rc[4]
        # circle = plt.Circle(
        #     (kde_cent[0], kde_cent[1]), KP_Bys_rc[3], color='g',
        #     fill=False, ls='dashed', lw=1.)
        b = rc_mode * np.sqrt(1. - ecc_mode**2)
        ellipse = mpatches.Ellipse(
            xy=kde_cent, width=2. * rc_mode, height=2. * b,
            angle=np.rad2deg(theta_mode), facecolor='None', edgecolor='g',
            ls='dashed', linewidth=1.5, transform=ax.transData)
        fig.gca().add_artist(ellipse)

    # # Plot contour levels if it was obtained.
    # if kde_plot:
    #     ext_range, x, y, k_pos = kde_plot
    #     # Number of contour lines depends on how large the cluster area is
    #     # compared with the area where the positional KDE was obtained.
    #     frac_xy = (ext_range[1] - ext_range[0]) / (2. * clust_rad)
    #     if frac_xy <= 1.5:
    #         c_lines = 5
    #     elif 1.5 < frac_xy <= 2.:
    #         c_lines = 10
    #     elif 2. < frac_xy <= 2.5:
    #         c_lines = 15
    #     elif 2.5 < frac_xy <= 3.:
    #         c_lines = 20
    #     elif 3 < frac_xy <= 3.5:
    #         c_lines = 25
    #     elif 3.5 < frac_xy:
    #         c_lines = 30
    #     kde = np.reshape(k_pos.T, x.shape)
    #     plt.contour(x, y, kde, c_lines, colors='b', linewidths=0.6)

    # Plot stars.
    plt.scatter(x_data, y_data, marker='o', c='black', s=st_sizes_arr,
                zorder=1)
    # Plot center.
    plt.scatter(kde_cent[0], kde_cent[1], color='r', s=40, lw=1.5,
                marker='x', zorder=5)


def pl_memb_vs_rad(
    gs, plot_style, coord, clust_rad, rad_radii, N_membs, N_membs_16,
        N_membs_84, rc, rt, kp_ndim, KP_plot):
    """
    """
    cent_dens, _16_84_rang, cd_rc_rt_sampled = 0., 0., (0., 0., 0.)
    if kp_ndim in (2, 4):
        cent_dens, _16_84_rang, cd_rc_rt_sampled = KP_plot['KP_cent_dens'],\
            KP_plot['_16_84_rang'], KP_plot['cd_rc_rt_sampled']

    # Convert from deg to arcmin
    cent_dens, cd_sampled = cent_dens / 3600., cd_rc_rt_sampled[0] / 3600.
    clust_rad, rc, rt, _16_84_rang, rc_sampled, rt_sampled,\
        rad_radii = clust_rad * 60., rc * 60.,\
        rt * 60., _16_84_rang * 60., cd_rc_rt_sampled[1] * 60.,\
        cd_rc_rt_sampled[2] * 60., rad_radii * 60.
    coord2 = 'arcmin'

    ax = plt.subplot(gs[4:6, 0:4])
    plt.xlabel(r'radius $[{}]$'.format(coord2))
    plt.ylabel("Number of members")
    if plot_style == 'asteca':
        ax.grid()

    # Minimum arcmin
    r_min = 1
    rads = np.linspace(r_min, rt, 50)
    # Plot King profile.
    if kp_ndim in (2, 4):
        N_KP_best, N_KP_16_84 = [], []
        for r in rads:
            N_KP_best.append(KP_memb_x(cent_dens, rc, rt, r))

            N_KP = []
            # Only use 100 samples to reduce the impact on performance
            for i, cd in enumerate(cd_sampled[:100]):
                rc_s, rt_s = rc_sampled[i], rt_sampled[i]
                N_KP.append(KP_memb_x(cd, rc_s, rt_s, r))
            N_KP_16_84.append(np.nanpercentile(N_KP, (16, 84)))
        N_KP_16_84 = np.array(N_KP_16_84).T
        # 16-84 region
        ax.fill_between(rads, N_KP_16_84[0], N_KP_16_84[1], facecolor='green',
                        alpha=0.4, label=r"$16th-84th$ (KP)")
        plt.axvline(rt, color='g', ls='-', label=r"$r_{t}$" + " (N={})".format(
            int(np.mean([N_KP_16_84[0][-1], N_KP_16_84[1][-1]]))))

        N_KP_best = np.array(N_KP_best)
        idx = np.argmin(abs(N_KP_best - .90 * N_KP_best.max()))
        plt.axvline(
            rads[idx], ls=':', label="90% " + r"$r_{t}$", c='orange', lw=2)
        idx = np.argmin(abs(N_KP_best - .95 * N_KP_best.max()))
        plt.axvline(
            rads[idx], ls=':', label="95% " + r"$r_{t}$", c='cyan', lw=2)
        idx = np.argmin(abs(N_KP_best - .99 * N_KP_best.max()))
        plt.axvline(
            rads[idx], ls=':', label="99% " + r"$r_{t}$", c='blue', lw=2)

    plt.axvline(clust_rad, color='r', ls='-', label=r"$r_{cl}$")
    plt.plot(rad_radii, N_membs, c='r', ls="--", label=r"$N_{{memb}}$")

    # 16-84 region
    if N_membs_16.any():
        ax.fill_between(
            rad_radii, N_membs_16, N_membs_84, facecolor='orange', alpha=0.2,
            label=r"$16th-84th$")
        N_membs_max = N_membs_84
    else:
        N_membs_max = N_membs

    if kp_ndim in (2, 4):
        rad_max = rt + rt * .2
    else:
        rad_max = clust_rad + clust_rad * .5
    rad_max = min(2 * rad_max, rad_radii[-1])

    idx = np.argmin(abs(rad_radii - rad_max))
    N_memb_ymax = min(N_membs_max[:idx].max(), N_membs.max() * 4)
    if kp_ndim in (2, 4):
        ymax = max(N_memb_ymax, N_KP_16_84[1][-1])
    else:
        ymax = N_memb_ymax

    plt.legend()
    plt.ylim(0, ymax + .1 * ymax)
    plt.xlim(rad_radii[0], rad_max)


def pl_membs_dist(gs, fig, members_dist, n_memb_i):
    """
    Coordinates 2D KDE zoom.
    """
    if members_dist.any():

        ax = plt.subplot(gs[4:6, 4:6])

        median, mad_v = np.median(members_dist), MAD(members_dist)
        msk = (members_dist > 0.) & (members_dist < median + 5. * mad_v)
        if msk.sum() > 10:
            members_dist = members_dist[msk]

        plt.hist(members_dist, 25, density=True, alpha=.5)
        xmin, xmax = ax.get_xlim()

        _16, _84 = np.percentile(members_dist, (16, 84))
        plt.axvline(n_memb_i, c='r', ls='-', lw=2,
                    label=r"$N_{{memb}}$={}".format(n_memb_i))
        plt.axvline(_16, c='orange', ls=':', lw=3)
        plt.axvline(_84, c='orange', ls=':', lw=3, label=r"$16th-84th$")
        plt.legend()

        ax.set_yticks([])
        plt.xlim(max(0., median - 5. * mad_v), xmax)
        plt.xlabel(r"Number of members (for $r_{{cl}}$)")


def plot(N, *args):
    """
    Handle each plot separately.
    """

    plt_map = {
        0: [pl_rad_find, 'Radius estimation'],
        1: [pl_mag_membs, 'estimated members vs magnitude cut'],
        2: [pl_cl_fl_regions, 'cluster and field regions defined'],
        3: [pl_rad_dens, 'radial density function'],
        4: [pl_zoom_frame, 'zoomed frame'],
        5: [pl_memb_vs_rad, 'estimated members vs radius'],
        6: [pl_membs_dist, 'distribution of the number of members']
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

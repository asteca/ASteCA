
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from itertools import cycle
from ..structure.king_profile import KingProf as kpf
from . import BayesPlots
from .. import aux_funcs


def pl_rad_find(
    gs, coord, clust_rad, e_rad, rad_rads, rad_N_membs, rad_N_field,
        rad_CI):
    """
    Radius estimation plot.
    """
    # Convert from deg to arcmin if (ra,dec) were used.
    if coord == 'deg':
        rad_rads = np.array(rad_rads) * 60.
        clust_rad, e_rad = clust_rad * 60., e_rad * 60.
        coord2 = 'arcmin'
    else:
        coord2 = 'px'

    # Cut when N_nembs = 0
    try:
        icut = np.where(rad_N_membs == 0.)[0][0]
        if icut == 0:
            icut = len(rad_N_membs)
        else:
            # Use 30 points since the number of radii values in
            # radius.rdpAreasDists() is hardcoded to 300.
            icut = icut + 30
    except IndexError:
        icut = len(rad_N_membs)
    rad_rads, rad_N_membs, rad_N_field, rad_CI = rad_rads[:icut],\
        rad_N_membs[:icut], rad_N_field[:icut], rad_CI[:icut]

    ax = plt.subplot(gs[0:2, 0:2])
    ax.minorticks_on()
    plt.xlabel(r'radius $[{}]$'.format(coord2), fontsize=10)
    plt.ylabel(r"$N_{memb}\;|\;N_{field}\;|\;CI$", fontsize=10)
    ax.tick_params(axis='both', which='major', labelsize=8)

    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=.5)

    plt.plot(rad_rads, rad_N_membs, c='g', label=r'$N_{memb}$')
    plt.plot(rad_rads, rad_N_field, c='b', ls='--', label=r'$N_{field}$')
    # ymin, ymax = ax.get_ylim()

    plt.plot(rad_rads, rad_CI, ls=':', c='k', label='CI')
    plt.axvline(x=clust_rad, lw=1.5, color='r', label=r"$r_{cl}$")
    if not np.isnan(e_rad[0]):
        plt.axvspan((e_rad[0]), (e_rad[1]), facecolor='grey', alpha=0.25)
    # Legends.
    leg = plt.legend(fancybox=True, fontsize=10)
    leg.get_frame().set_alpha(0.7)

    xmin, xmax = ax.get_xlim()
    plt.xlim(max(xmin, 0.), xmax)
    plt.ylim(0., 1.02)
    # plt.ylim(ymin, ymax)


def pl_rad_dens(
    gs, coord, rdp_radii, rdp_points, rdp_stddev, field_dens, e_fdens,
    clust_rad, e_rad, kp_flag, KP_Bys_rc, KP_Bys_rt, KP_cent_dens,
        KP_conct_par):
    """
    Radial density plot.
    """

    # Convert from deg to arcmin if (ra,dec) were used.
    if coord == 'deg':
        rdp_radii = np.array(rdp_radii) * 60.
        clust_rad, e_rad = clust_rad * 60., e_rad * 60.
        KP_Bys_rc, KP_Bys_rt = KP_Bys_rc * 60., KP_Bys_rt * 60.
        field_dens, e_fdens, KP_cent_dens = field_dens / 3600.,\
            e_fdens / 3600., KP_cent_dens / 3600.
        rdp_points = np.array(rdp_points) / 3600.
        rdp_stddev = np.array(rdp_stddev) / 3600.
        coord2 = 'arcmin'
    else:
        coord2 = 'px'

    ax = plt.subplot(gs[0:2, 2:6])
    # Get max and min values in x,y
    x_min = max(min(rdp_radii) - (max(rdp_radii) / 20.), 0)
    x_max = min(max(rdp_radii) + (max(rdp_radii) / 20.), 3. * clust_rad)

    N_h = int(.5 * len(rdp_points))
    delta_total = (max(rdp_points) - field_dens)
    delta_backg = 0.1 * delta_total
    y_min = max((field_dens - delta_backg) -
                (max(rdp_points[:N_h]) - min(rdp_points[:N_h])) / 10., 0)
    y_max = max(rdp_points[:N_h]) + (
        max(rdp_points[:N_h]) - min(rdp_points[:N_h])) / 10.
    # Set plot limits
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    ax.minorticks_on()
    # Set axes labels
    plt.xlabel(r'radius $[{}]$'.format(coord2), fontsize=10)
    plt.ylabel(r"$\rho$ $[st/{}^{{2}}]$".format(coord2), fontsize=10)
    ax.tick_params(axis='both', which='major', labelsize=8)
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=.5)

    # Legend texts
    r_frmt = '{:.0f}' if coord2 == 'px' else '{:.2f}'
    t_rad = r"$r_{{{}}}=" + r_frmt + r"_{{" + r_frmt + r"}}^{{" +\
        r_frmt + r"}}\,[{}]$"

    # Plot density profile
    ax.plot(rdp_radii, rdp_points, marker='o', ms=5, lw=1., zorder=3,
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
    # Plot radius.
    y_mid_point = (y_max + y_min) * .5
    ax.vlines(x=clust_rad, ymin=field_dens, ymax=y_mid_point, lw=1.5,
              color='r', label=t_rad.format(
                  "cl", clust_rad, e_rad[0], e_rad[1], coord2), zorder=5)
    # Plot radius error zone.
    if not np.isnan(e_rad[0]):
        plt.axvspan(e_rad[0], e_rad[1], facecolor='grey', alpha=0.25)

    # Plot King profile.
    if kp_flag:
        txts = [
            'King prof ({:.2f})'.format(KP_conct_par),
            t_rad.format(
                "c", KP_Bys_rc[1], KP_Bys_rc[0], KP_Bys_rc[2], coord2),
            t_rad.format("t", KP_Bys_rt[1], KP_Bys_rt[0], KP_Bys_rt[2], coord2)
        ]
        # Plot curve. Values outside of rt contribute 'fd'.
        kpf_xvals = np.linspace(rdp_radii[0], KP_Bys_rt[1], 100)
        kpf_yvals = KP_cent_dens * kpf(
            kpf_xvals, KP_Bys_rc[1], KP_Bys_rt[1]) + field_dens
        ax.plot(kpf_xvals, kpf_yvals, 'g--', label=txts[0], lw=2., zorder=3)
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
    ax.legend(handles, labels, loc='upper center', numpoints=2, fontsize=10)

    #
    # Log-log plot
    # ax = plt.subplot(gs[2:4, 4:6])
    axins = inset_axes(ax, width=2.5, height=2.5)
    axins.minorticks_on()
    axins.grid(b=True, which='both', color='gray', linestyle='--', lw=.25)
    # plt.xlabel(r'log(radius) $[{}]$'.format(coord2), fontsize=8)
    # plt.ylabel(r"log($\rho$) $[st/{}^{{2}}]$".format(coord2), fontsize=8)
    axins.tick_params(axis='both', which='both', labelsize=6)

    axins.scatter(rdp_radii, rdp_points, s=5, zorder=5)
    axins.hlines(y=field_dens, xmin=min(rdp_radii), xmax=max(rdp_radii),
                 color='k', ls='--', zorder=5)
    if not np.isnan(e_fdens):
        axins.hlines(
            y=field_dens - e_fdens, xmin=min(rdp_radii),
            xmax=max(rdp_radii), color='k', ls=':', zorder=5)
        axins.hlines(
            y=field_dens + e_fdens, xmin=min(rdp_radii),
            xmax=max(rdp_radii), color='k', ls=':', zorder=5)
    axins.vlines(
        x=clust_rad, ymin=field_dens, ymax=y_mid_point, lw=1.5, color='r',
        zorder=10)
    # Plot King profile.
    if kp_flag:
        axins.plot(kpf_xvals, kpf_yvals, 'g--', lw=1., zorder=3)
        axins.vlines(
            x=KP_Bys_rc[1], ymin=field_dens, ymax=rc_ymax, color='g',
            linestyles=':', lw=1.)
        axins.vlines(
            x=KP_Bys_rt[1], ymin=field_dens, ymax=y_mid_point, color='g', lw=1)

    axins.set_xscale('log')
    axins.set_yscale('log')


def pl_KP_Bys(
    gs, coord, kp_flag, kp_nburn, KP_steps, KP_mean_afs, KP_tau_autocorr,
        KP_ESS, KP_samples, KP_Bys_rc, KP_Bys_rt):
    """
    """
    if kp_flag:

        # Convert from deg to arcmin if (ra,dec) were used.
        if coord == 'deg':
            KP_samples = KP_samples * 60.
            KP_Bys_rc, KP_Bys_rt = KP_Bys_rc * 60., KP_Bys_rt * 60.
            coord2 = 'arcmin'
        else:
            coord2 = 'px'

        plt.style.use('seaborn-darkgrid')

        gsy, gsx = (2, 3), (0, 2)
        BayesPlots.autocorr(
            gs, gsx, gsy, KP_steps, KP_tau_autocorr, KP_ESS)

        gsy, gsx = (3, 4), (0, 2)
        BayesPlots.meanAF(gs, gsx, gsy, KP_steps, KP_mean_afs)

        gsy, gsx = (2, 3), (2, 6)
        xylabel = r"$r_{{c}}$ [{}]".format(coord2)
        BayesPlots.traceplot(
            gs, gsx, gsy, KP_samples[:, :, 0], KP_Bys_rc, kp_nburn, xylabel,
            False)

        gsy, gsx = (3, 4), (2, 6)
        xylabel = r"$r_{{t}}$ [{}]".format(coord2)
        BayesPlots.traceplot(
            gs, gsx, gsy, KP_samples[:, :, 1], KP_Bys_rt, kp_nburn, xylabel)

        gsy, gsx = (4, 6), (0, 2)
        xylabel = r"$r_{{c}}$ [{}]".format(coord2)
        mu_kde_x, mu_kde = aux_funcs.kde1D(KP_samples[:, :, 0].flatten())
        BayesPlots.histogram(
            gs, gsx, gsy, KP_samples[:, :, 0], mu_kde_x, mu_kde, KP_Bys_rc,
            xylabel)

        gsy, gsx = (4, 6), (2, 4)
        xylabel = r"$r_{{t}}$ [{}]".format(coord2)
        mu_kde_x, mu_kde = aux_funcs.kde1D(KP_samples[:, :, 1].flatten())
        BayesPlots.histogram(
            gs, gsx, gsy, KP_samples[:, :, 1], mu_kde_x, mu_kde, KP_Bys_rt,
            xylabel)

        gsy, gsx = (4, 6), (4, 6)
        xylabel = (
            r"$r_{{c}}$ [{}]".format(coord2), r"$r_{{t}}$ [{}]".format(coord2))
        BayesPlots.twoParDens(
            gs, gsx, gsy, KP_samples, KP_Bys_rc, KP_Bys_rt, xylabel)

        plt.style.use('default')


def pl_zoom_frame(
    gs, fig, x_name, y_name, coord, x_zmin, x_zmax, y_zmin, y_zmax, cont_index,
    x_data, y_data, st_sizes_arr, kde_cent, clust_rad, KP_Bys_rc,
        KP_Bys_rt, frac_cl_area, kp_flag):
    """
    Zoom on x,y finding chart.
    """
    if kp_flag:
        ax = plt.subplot(gs[6:8, 0:2])
    else:
        ax = plt.subplot(gs[2:4, 0:2])

    # Force square plot.
    ax.set_aspect('equal')
    # Set plot limits.
    plt.xlim(x_zmin, x_zmax)
    plt.ylim(y_zmin, y_zmax)
    if coord == 'deg':
        # If RA is used, invert axis.
        ax.invert_xaxis()
    # Set axis labels
    plt.xlabel('{} ({})'.format(x_name, coord), fontsize=10)
    plt.ylabel('{} ({})'.format(y_name, coord), fontsize=10)
    ax.tick_params(axis='both', which='major', labelsize=8)

    # Set minor ticks
    ax.minorticks_on()
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
            100. * frac_cl_area, N_in, cont_index), fontsize=9)
    fig.gca().add_artist(circle)

    # Core and tidal radii
    if kp_flag:
        # Plot tidal radius.
        circle = plt.Circle(
            (kde_cent[0], kde_cent[1]), KP_Bys_rt[1], color='g', fill=False,
            lw=1.5)
        fig.gca().add_artist(circle)
        # Plot core radius.
        if KP_Bys_rc[1] > 0:
            circle = plt.Circle(
                (kde_cent[0], kde_cent[1]), KP_Bys_rc[1], color='g',
                fill=False, ls='dashed', lw=1.)
            fig.gca().add_artist(circle)

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


def pl_mag_membs(gs, y_ax, kp_flag, membvsmag):
    """
    """
    if kp_flag:
        ax = plt.subplot(gs[6:8, 2:4])
    else:
        ax = plt.subplot(gs[2:4, 2:4])
    ax.minorticks_on()
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=.5,
            zorder=0)
    ax.set_title(r"$N_{{memb}}$ vs magnitude cut (phot incomp)", fontsize=9)
    plt.xlabel('$' + y_ax + '$', fontsize=10)
    plt.ylabel(r'$N_{memb}$', fontsize=10)
    ax.tick_params(axis='both', which='major', labelsize=8)
    plt.bar(*membvsmag, zorder=5)


def pl_cl_fl_regions(
    gs, fig, x_name, y_name, coord, x_min, x_max, y_min, y_max, asp_ratio,
        kde_cent, clust_rad, field_regions_i, field_regions_rjct_i,
        cl_region_i, cl_region_rjct_i, flag_no_fl_regs_i, kp_flag):
    """
    Cluster and field regions defined.
    """
    if kp_flag:
        ax = plt.subplot(gs[6:8, 4:6])
    else:
        ax = plt.subplot(gs[2:4, 4:6])

    ax.set_aspect(aspect=asp_ratio)
    # Set plot limits
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    # If RA is used, invert axis.
    if coord == 'deg':
        ax.invert_xaxis()
    # Set axis labels
    plt.xlabel('{} ({})'.format(x_name, coord), fontsize=10)
    plt.ylabel('{} ({})'.format(y_name, coord), fontsize=10)
    ax.tick_params(axis='both', which='major', labelsize=8)

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
    """
    Handle each plot separately.
    """

    plt_map = {
        0: [pl_rad_find, 'Radius estimation'],
        1: [pl_rad_dens, 'radial density function'],
        2: [pl_KP_Bys, 'King profile Bayes plots'],
        3: [pl_zoom_frame, 'zoomed frame'],
        4: [pl_mag_membs, 'estimated members vs magnitude cut'],
        5: [pl_cl_fl_regions, 'cluster and field regions defined']
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

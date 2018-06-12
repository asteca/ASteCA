
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib.offsetbox as offsetbox
import numpy as np


def pl_mp_histo(
        gs, n_memb_da, memb_prob_avrg_sort, flag_decont_skip, cl_reg_fit,
        min_prob, mode_fld_clean, local_bin):
    '''
    Histogram for the distribution of membership probabilities from the
    decontamination algorithm.
    '''
    # Only attempt to plot if the DA was applied.
    if flag_decont_skip is False:
        # Reduced membership.
        ax = plt.subplot(gs[0:2, 0:2])
        plt.xlim(0., 1.)
        plt.xlabel('MP (membership probability)', fontsize=12)
        plt.ylabel('N (normalized)', fontsize=12)
        ax.minorticks_on()
        ax.grid(b=True, which='major', color='gray', linestyle='--', lw=1,
                zorder=1)
        prob_data = [star[9] for star in memb_prob_avrg_sort]
        # Histogram of the data.
        n_bins = int((max(prob_data) - min(prob_data)) / 0.025)
        if n_bins > 0:
            # Normalized histogram.
            n, bins, patches = plt.hist(prob_data, n_bins, density=True)
            # Get bin centers.
            bin_centers = 0.5 * (bins[:-1] + bins[1:])
            # scale values to interval [0,1]
            col = bin_centers - min(bin_centers)
            col /= max(col)
            cm = plt.cm.get_cmap('RdYlBu_r')
            # Plot histo colored according to colormap.
            for c, p in zip(col, patches):
                plt.setp(p, 'facecolor', cm(c), zorder=3)
                plt.setp(p, 'edgecolor', 'k')
        else:
            print("  WARNING: all MPs are equal valued. "
                  "Can not plot MPs histogram.")
        # Plot minimum probability line.
        plt.axvline(x=min_prob, linestyle='--', color='green', lw=2.5,
                    zorder=3)

        # Add text box.
        str_pm = ['MP', '\geq', 'prob']
        if mode_fld_clean == 'local':
            str_pm.append(mode_fld_clean + ';\,' + local_bin)
        else:
            str_pm.append(mode_fld_clean.replace('_', '\_'))
        text0 = r'$N_{{total}}={}$'.format(len(prob_data))
        text1 = r'$n_{{memb-DA}}={}\,(MP \geq 0.5)$'.format(n_memb_da)
        text2 = r'${}_{{min}}={:.2f}\,({})$'.format(str_pm[2], min_prob,
                                                    str_pm[3])
        text3 = r'$N_{{fit}}={} \, ({} {} {}_{{min}})$'.format(
            len(cl_reg_fit), str_pm[0], str_pm[1], str_pm[2])
        text = text0 + '\n' + text1 + '\n' + text2 + '\n' + text3
        ob = offsetbox.AnchoredText(text, loc=6, prop=dict(size=10))
        ob.patch.set(boxstyle='square,pad=0.05', alpha=0.85)
        ax.add_artist(ob)
        # Avoid showing the value 0.0 in the y axis.
        plt.ylim(0.001, plt.ylim()[1])


def pl_chart_mps(gs, fig, x_name, y_name, coord, x_zmin, x_zmax, y_zmin,
                 y_zmax, kde_cent, clust_rad, flag_decont_skip, v_min_mp,
                 v_max_mp, chart_fit_inv, chart_no_fit_inv, out_clust_rad,
                 mode_fld_clean, local_bin):
    '''
    Finding chart of cluster region with decontamination algorithm
    applied and colors assigned according to the probabilities obtained.
    '''
    ax = plt.subplot(gs[0:2, 2:4])
    # Set plot limits, Use 'zoom' x,y ranges.
    plt.xlim(x_zmin, x_zmax)
    plt.ylim(y_zmin, y_zmax)
    ax.set_title('Cluster region'.format(), fontsize=9)
    # If RA is used, invert axis.
    if coord == 'deg':
        ax.invert_xaxis()
    # Set axis labels
    plt.xlabel('{} ({})'.format(x_name, coord), fontsize=12)
    plt.ylabel('{} ({})'.format(y_name, coord), fontsize=12)
    # Set minor ticks
    ax.minorticks_on()
    # Radius
    circle = plt.Circle((kde_cent[0], kde_cent[1]), clust_rad, color='red',
                        fill=False)
    fig.gca().add_artist(circle)
    # If DA was skipped, print info on 'local' method here.
    if flag_decont_skip and mode_fld_clean == 'local':
        text = r'$({})$'.format(mode_fld_clean + ';\,' + local_bin)
        ob = offsetbox.AnchoredText(text, pad=0.2, loc=2, prop=dict(size=12))
        ob.patch.set(alpha=0.85)
        ax.add_artist(ob)
    # Color map, higher prob stars look redder.
    cm = plt.cm.get_cmap('RdYlBu_r')
    # If stars have a range of colors, use list of colors. Else use a single
    # color.
    if v_min_mp != v_max_mp:
        col_select_fit, col_select_no_fit = chart_fit_inv[2], \
            chart_no_fit_inv[2]
    else:
        col_select_fit, col_select_no_fit = '#4682b4', '#4682b4'
    # Plot stars *not* used in the best fit process.
    plt.scatter(
        chart_no_fit_inv[0], chart_no_fit_inv[1], marker='o',
        c=col_select_no_fit, s=30, edgecolors='black', cmap=cm,
        alpha=0.5, lw=0.5, vmin=v_min_mp, vmax=v_max_mp)
    # Add line to stars not used in the best bit process.
    plt.scatter(chart_no_fit_inv[0], chart_no_fit_inv[1], marker='_',
                c='k', lw=0.5, alpha=0.5)
    # Plot stars selected to be used in the best bit process.
    plt.scatter(chart_fit_inv[0], chart_fit_inv[1], marker='o',
                c=col_select_fit, s=30, edgecolors='black', cmap=cm,
                lw=0.5, vmin=v_min_mp, vmax=v_max_mp, zorder=4)
    # Plot stars outside the cluster's radius as empty circles.
    plt.scatter(out_clust_rad[0], out_clust_rad[1], marker='o',
                s=30, edgecolors='black', facecolors='none', lw=0.5)


def pl_mps_phot_diag(gs, fig, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd,
                     x_ax, y_ax, v_min_mp, v_max_mp, diag_fit_inv,
                     diag_no_fit_inv, err_bar, mode_fld_clean, bin_edges):
    '''
    Star's membership probabilities on cluster's photometric diagram.
    '''
    x_val, mag_y, xy_err = err_bar
    ax = plt.subplot(gs[0:2, 4:6])
    # Set plot limits
    plt.xlim(x_min_cmd, x_max_cmd)
    plt.ylim(y_min_cmd, y_max_cmd)
    # Set axis labels
    plt.xlabel('$' + x_ax + '$', fontsize=18)
    plt.ylabel('$' + y_ax + '$', fontsize=18)
    # Add text box.
    text = '$N_{{fit}}={}$'.format(len(diag_fit_inv[2]))
    ob = offsetbox.AnchoredText(text, loc=2, prop=dict(size=14))
    ob.patch.set(boxstyle='square,pad=-0.2', alpha=0.85)
    ax.add_artist(ob)
    # Set minor ticks
    ax.minorticks_on()
    ax.xaxis.set_major_locator(MultipleLocator(1.0))
    # Plot grid. If bin_edges == 0., it means the 'local' method was not used.
    if mode_fld_clean == 'local' and bin_edges != 0.:
        # TODO using first magnitude and color. Generalize to N-dimensions.
        for x_ed in bin_edges[1]:
            # vertical lines
            ax.axvline(x_ed, linestyle=':', lw=.8, color='k', zorder=1)
        for y_ed in bin_edges[0]:
            # horizontal lines
            ax.axhline(y_ed, linestyle=':', lw=.8, color='k', zorder=1)
    else:
        ax.grid(b=True, which='major', color='gray', linestyle='--', lw=1,
                zorder=1)
    # This reversed colormap means higher prob stars will look redder.
    cm = plt.cm.get_cmap('RdYlBu_r')
    # If stars have a range of colors, use list of colors. Else use a single
    # color.
    if v_min_mp != v_max_mp:
        col_select_fit, col_select_no_fit = diag_fit_inv[2], \
            diag_no_fit_inv[2]
    else:
        col_select_fit, col_select_no_fit = '#4682b4', '#4682b4'
    # Plot stars *not* used in the best fit process.
    plt.scatter(diag_no_fit_inv[1][0], diag_no_fit_inv[0][0], marker='o',
                c=col_select_no_fit, s=25, cmap=cm, lw=0.5, edgecolor='k',
                alpha=0.5, vmin=v_min_mp, vmax=v_max_mp, zorder=2)
    # Draw horizontal line over stars discarded.
    plt.scatter(diag_no_fit_inv[1][0], diag_no_fit_inv[0][0],
                marker='_', c='k', lw=0.5, alpha=0.5, zorder=3)
    # Plot stars used in the best fit process.
    sca = plt.scatter(diag_fit_inv[1][0], diag_fit_inv[0][0], marker='o',
                      c=col_select_fit, s=30, cmap=cm, lw=0.5, edgecolor='k',
                      vmin=v_min_mp, vmax=v_max_mp, zorder=4)
    # If list is not empty, plot error bars at several values.
    if x_val:
        plt.errorbar(
            x_val, mag_y, yerr=xy_err[0], xerr=xy_err[1], fmt='k.', lw=0.8,
            ms=0., zorder=4)
    # For plotting the colorbar (see bottom of make_plots file).
    trans = ax.transAxes + fig.transFigure.inverted()

    return sca, trans


def pl_plx_histo(gs, plx_flag, plx_clrg, plx_xmin, plx_xmax, plx_x_kde, kde_pl,
                 plx_flrg, flag_no_fl_regs_i):
    '''
    Histogram for the distribution of parallaxes within the cluster region.
    '''
    if plx_flag:
        ax = plt.subplot(gs[2:4, 0:2])
        plt.xlim(plx_xmin, plx_xmax)
        plt.xlabel('Plx [mas]', fontsize=12)
        plt.ylabel('N', fontsize=12)
        ax.minorticks_on()
        ax.grid(b=True, which='major', color='gray', linestyle='--', lw=1,
                zorder=1)
        # Normalized histogram for cluster region.
        plt.hist(
            plx_clrg, 120, density=True, zorder=4, color='#9aafd1',
            label="Cluster region (fit)")
        # Plot histogram for the parallaxes of the field regions.
        if not flag_no_fl_regs_i:
            plt.hist(
                plx_flrg, 120, density=True, zorder=4, color='#ef703e',
                label="Field regions", alpha=0.5)

        plt.plot(plx_x_kde, kde_pl / max(kde_pl), color='b', lw=1., zorder=4)
        # Maximum KDE value.
        p_max_mas = plx_x_kde[np.argmax(kde_pl)]
        plt.axvline(x=p_max_mas, linestyle='--', color='b', lw=.7, zorder=5)
        d_max_pc = 1000. / p_max_mas
        plx_lt_zero = 100. * plx_clrg[plx_clrg < 0.].size / plx_clrg.size
        ob = offsetbox.AnchoredText(
            r"$Plx_{{max}}$={:.3f} [mas]".format(p_max_mas) +
            "\n({:.0f} [pc])\n".format(d_max_pc) +
            r"$Plx<0 \rightarrow$ {:.1f}%".format(plx_lt_zero),
            pad=0.2, loc=1, prop=dict(size=9))
        ob.patch.set(alpha=0.85)
        ax.add_artist(ob)
        # Avoid showing the value 0.0 in the y axis.
        plt.ylim(0.01, plt.ylim()[1])
        ax.legend(fontsize='small', loc=7)


def pl_plx_chart(gs, plx_flag, x_name, y_name, coord, cl_reg_fit, plx_x_kde,
                 kde_pl):
    '''
    Finding chart of cluster region with colors assigned according to the
    probabilities obtained and sizes according to parallaxes.
    '''
    if plx_flag:
        ax = plt.subplot(gs[2:4, 2:4])
        ax.grid(b=True, which='major', color='gray', linestyle='--', lw=1,
                zorder=1)
        ax.set_title('Cluster region (fit)'.format(), fontsize=9)

        # If RA is used, invert axis.
        if coord == 'deg':
            ax.invert_xaxis()
        # Set axis labels
        plt.xlabel('{} ({})'.format(x_name, coord), fontsize=12)
        plt.ylabel('{} ({})'.format(y_name, coord), fontsize=12)
        # Set minor ticks
        ax.minorticks_on()

        # Prepare data.
        x = np.array(zip(*cl_reg_fit)[1])
        y = np.array(zip(*cl_reg_fit)[2])
        mp = np.array(zip(*cl_reg_fit)[9])
        plx = np.array(zip(*zip(*cl_reg_fit)[7])[0])
        msk = (~np.isnan(x)) & (~np.isnan(y)) & (~np.isnan(mp)) &\
            (~np.isnan(plx))
        x, y, mp, plx = x[msk], y[msk], mp[msk], plx[msk]
        # Clip parallax values.
        np.clip(plx, a_min=0., a_max=10., out=plx)

        # Maximum KDE parallax value.
        p_max_mas = plx_x_kde[np.argmax(kde_pl)]
        # Distance to max value. Stars closer to the max value are larger.
        plx_d = 2. + 1. / (abs(plx - p_max_mas) + .1) ** 2.3

        # Re-arrange so stars closer to the max Plx value are on top
        plx_i = plx_d.argsort()
        x, y, mp, plx_d = x[plx_i], y[plx_i], mp[plx_i], plx_d[plx_i]

        # Color map, higher prob stars look redder.
        cm = plt.cm.get_cmap('viridis')  # RdYlBu_r
        # Plot stars selected to be used in the best fit process.
        plt.scatter(
            x, y, marker='o', c=mp, s=plx_d, edgecolors='black',
            cmap=cm, lw=0.35, zorder=4)

        ob = offsetbox.AnchoredText(
            "MP: ({:.2f}, {:.2f})".format(np.min(mp), np.max(mp)),
            pad=0.2, loc=1, prop=dict(size=9))
        ob.patch.set(alpha=0.85)
        ax.add_artist(ob)


def pl_plx_vs_MP(gs, plx_flag, cl_reg_fit, plx_x_kde, kde_pl):
    '''
    Finding chart of cluster region with colors assigned according to the
    probabilities obtained and sizes according to parallaxes.
    '''
    if plx_flag:
        ax = plt.subplot(gs[2:4, 4:6])
        ax.grid(b=True, which='major', color='gray', linestyle='--', lw=1,
                zorder=1)

        # Set axis labels
        plt.xlabel('Plx [mas]', fontsize=12)
        plt.ylabel('MP', fontsize=12)
        # Set minor ticks
        ax.minorticks_on()

        # Prepare data.
        mp = np.array(zip(*cl_reg_fit)[9])
        plx = np.array(zip(*zip(*cl_reg_fit)[7])[0])
        e_plx = np.array(zip(*zip(*cl_reg_fit)[8])[0])

        # Re-arrange so stars closer to the max Plx value are on top
        plx_i = plx.argsort()
        mp, plx, e_plx = mp[plx_i], plx[plx_i], e_plx[plx_i]

        # Color map, higher prob stars look redder.
        cm = plt.cm.get_cmap('viridis')
        # Plot stars selected to be used in the best fit process.
        plt.scatter(
            plx, mp, marker='o', c=mp, s=30, edgecolors='black',
            cmap=cm, lw=0.35, zorder=4)
        ax.errorbar(
            plx, mp, xerr=e_plx, fmt='none', elinewidth=.35, ecolor='grey')

        plt.plot(plx_x_kde, kde_pl / max(kde_pl), color='k', lw=1., zorder=4)
        # Maximum KDE value.
        p_max_mas = plx_x_kde[np.argmax(kde_pl)]
        plt.axvline(x=p_max_mas, linestyle='--', color='r', lw=.85, zorder=5)

        plt.xlim(0., np.median(plx) + 2. * np.std(plx))
        plt.ylim(np.min(mp), 1.01)


def pl_pms_plot(gs, coord, plx_flag, cl_reg_fit):
    '''
    Finding chart of cluster region with colors assigned according to the
    probabilities obtained and sizes according to parallaxes.
    '''

    pmRA = np.array(zip(*zip(*cl_reg_fit)[7])[1])
    if pmRA.any():
        if plx_flag:
            ax = plt.subplot(gs[4:6, 0:2])
        else:
            ax = plt.subplot(gs[2:4, 0:2])
        ax.grid(b=True, which='major', color='gray', linestyle='--', lw=1,
                zorder=1)

        # Set minor ticks
        ax.minorticks_on()

        # Prepare data.
        mp = np.array(zip(*cl_reg_fit)[9])
        e_pmRA, pmDE, e_pmDE = \
            np.array(zip(*zip(*cl_reg_fit)[8])[1]),\
            np.array(zip(*zip(*cl_reg_fit)[7])[2]),\
            np.array(zip(*zip(*cl_reg_fit)[8])[2])

        if coord == 'deg':
            plt.xlabel(
                r"$\mu_{{\alpha}} \, cos \delta \, [mas/yr]$", fontsize=12)
            DE = np.array(zip(*zip(*cl_reg_fit)[7])[2])
        else:
            plt.xlabel(r"$\mu_{{\alpha}} \, [mas/yr]$", fontsize=12)
            DE = np.zeros(pmRA.size)
        plt.ylabel(r"$\mu_{{\delta}} \, [mas/yr]$", fontsize=12)

        # msk = (~np.isnan(x)) & (~np.isnan(y)) & (~np.isnan(mp)) &\
        #     (~np.isnan(plx))
        # x, y, mp, plx = x[msk], y[msk], mp[msk], plx[msk]
        # # Clip parallax values.
        # np.clip(plx, a_min=0., a_max=10., out=plx)

        # Re-arrange so stars closer to the max Plx value are on top
        mp_i = mp.argsort()
        pmRA, e_pmRA, pmDE, e_pmDE, DE, mp = pmRA[mp_i], e_pmRA[mp_i],\
            pmDE[mp_i], e_pmDE[mp_i], DE[mp_i], mp[mp_i]

        # Color map, higher prob stars look redder.
        # cm = plt.cm.get_cmap('RdYlBu_r')

        import matplotlib.cm as cm
        from matplotlib.colors import Normalize
        cmap = cm.viridis
        norm = Normalize(vmin=mp.min(), vmax=mp.max())

        # Plot stars selected to be used in the best bit process.
        # plt.scatter(
        #     pmRA * np.cos(np.deg2rad(DE)), pmDE, marker='o', c=mp, s=plx_d,
        #     edgecolors='black',
        #     cmap=cm, lw=0.35, zorder=4)

        ax.errorbar(
            pmRA * np.cos(np.deg2rad(DE)), pmDE, yerr=e_pmDE, xerr=e_pmRA,
            fmt='none', elinewidth=.35, ecolor=cmap(norm(mp)))


def plot(N, *args):
    '''
    Handle each plot separately.
    '''

    plt_map = {
        0: [pl_mp_histo, 'MPs histogram'],
        1: [pl_chart_mps, 'frame with MPs coloring'],
        2: [pl_plx_histo, 'Plx histogram'],
        3: [pl_plx_chart, 'Plx chart'],
        4: [pl_plx_vs_MP, 'Plx vs MP'],
        5: [pl_pms_plot, 'PMs plot']
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


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.offsetbox as offsetbox
from . prep_plots import xylabelsize, xytickssize, titlesize, legendsize


def pl_mp_histo(
    gs, n_memb_da, memb_prob_avrg_sort, flag_decont_skip, cl_reg_fit,
        mode_fld_clean, local_bin):
    """
    Histogram for the distribution of membership probabilities from the
    decontamination algorithm.
    """
    # Only attempt to plot if the DA was applied.
    if flag_decont_skip is False:
        # Reduced membership.
        ax = plt.subplot(gs[0:2, 0:2])
        plt.xlim(0., 1.)
        plt.xlabel('MP', fontsize=xylabelsize)
        plt.ylabel('N (normalized)', fontsize=xylabelsize)
        ax.minorticks_on()
        ax.tick_params(axis='both', which='major', labelsize=xytickssize)
        ax.grid(b=True, which='major', color='gray', linestyle='--', lw=.5,
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
            for c, p in list(zip(col, patches)):
                plt.setp(p, 'facecolor', cm(c), zorder=3)
                plt.setp(p, 'edgecolor', 'k')
        else:
            print("  WARNING: all MPs are equal valued. "
                  "Can not plot MPs histogram.")
        # Plot minimum probability line.
        min_prob = cl_reg_fit[-1][-1]
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
        text2 = r'${}_{{min}}={:.2f}\,({})$'.format(
            str_pm[2], min_prob, str_pm[3])
        text3 = r'$N_{{fit}}={} \, ({} {} {}_{{min}})$'.format(
            len(cl_reg_fit), str_pm[0], str_pm[1], str_pm[2])
        text = text0 + '\n' + text1 + '\n' + text2 + '\n' + text3
        plt.plot([], label=text)
        plt.legend(fontsize=legendsize, handlelength=0)
        # Avoid showing the value 0.0 in the y axis.
        plt.ylim(0.001, plt.ylim()[1])


def pl_chart_mps(
    gs, fig, x_name, y_name, coord, x_zmin, x_zmax, y_zmin, y_zmax, kde_cent,
    clust_rad, flag_decont_skip, v_min_mp, v_max_mp, chart_fit_inv,
        chart_no_fit_inv, out_clust_rad, mode_fld_clean, local_bin):
    """
    Finding chart of cluster region with decontamination algorithm
    applied and colors assigned according to the probabilities obtained.
    """
    ax = plt.subplot(gs[0:2, 2:4])
    # Set plot limits, Use 'zoom' x,y ranges.
    plt.xlim(x_zmin, x_zmax)
    plt.ylim(y_zmin, y_zmax)
    ax.set_title('Cluster region'.format(), fontsize=titlesize)
    # If RA is used, invert axis.
    if coord == 'deg':
        ax.invert_xaxis()
    # Set axis labels
    plt.xlabel('{} ({})'.format(x_name, coord), fontsize=xylabelsize)
    plt.ylabel('{} ({})'.format(y_name, coord), fontsize=xylabelsize)
    # Set minor ticks
    ax.minorticks_on()
    ax.tick_params(axis='both', which='major', labelsize=xytickssize)
    # Radius
    circle = plt.Circle((kde_cent[0], kde_cent[1]), clust_rad, color='red',
                        fill=False)
    fig.gca().add_artist(circle)
    # If DA was skipped, print info on 'local' method here.
    if flag_decont_skip and mode_fld_clean == 'local':
        text = r'$({})$'.format(mode_fld_clean + ';\,' + local_bin)
        ob = offsetbox.AnchoredText(
            text, pad=0.2, loc=2, prop=dict(size=legendsize))
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


def pl_mps_phot_diag(
    gs, fig, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd, x_ax, y_ax, v_min_mp,
    v_max_mp, diag_fit_inv, diag_no_fit_inv, err_bar, mode_fld_clean,
        local_rm_edges):
    """
    Star's membership probabilities on cluster's photometric diagram.
    """
    x_val, mag_y, xy_err = err_bar
    ax = plt.subplot(gs[0:2, 4:6])
    # Set plot limits
    plt.xlim(x_min_cmd, x_max_cmd)
    plt.ylim(y_min_cmd, y_max_cmd)
    # Set axis labels
    plt.xlabel('$' + x_ax + '$', fontsize=xylabelsize)
    plt.ylabel('$' + y_ax + '$', fontsize=xylabelsize)
    # Add text box.
    text = '$N_{{fit}}={}$'.format(len(diag_fit_inv[2]))
    ob = offsetbox.AnchoredText(text, loc=2, prop=dict(size=legendsize))
    ob.patch.set(boxstyle='square,pad=-0.2', alpha=0.85)
    ax.add_artist(ob)
    # Set minor ticks
    ax.minorticks_on()
    ax.tick_params(axis='both', which='major', labelsize=xytickssize)
    # Plot grid.
    if mode_fld_clean == 'local' and local_rm_edges is not None:
        # Use first magnitude and color.
        for x_ed in local_rm_edges[1]:
            # vertical lines
            ax.axvline(x_ed, linestyle=':', lw=.8, color='k', zorder=1)
        for y_ed in local_rm_edges[0]:
            # horizontal lines
            ax.axhline(y_ed, linestyle=':', lw=.8, color='k', zorder=1)
    else:
        ax.grid(b=True, which='major', color='gray', linestyle='--', lw=.5,
                zorder=1)
    # This reversed colormap means higher prob stars will look redder.
    cm = plt.cm.get_cmap('RdYlBu_r')
    # If stars have a range of colors, use list of colors. Else use a single
    # color.
    if v_min_mp != v_max_mp:
        col_select_fit, col_select_no_fit = diag_fit_inv[2], \
            diag_no_fit_inv[2]
        # Decide if colorbar should be plotted.
        plot_colorbar = True
    else:
        col_select_fit, col_select_no_fit = '#4682b4', '#4682b4'
        plot_colorbar = False
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
    # For plotting the colorbar
    trans = ax.transAxes + fig.transFigure.inverted()

    return plot_colorbar, sca, trans


def pl_mps_incomp_diags(
    gs, fig, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd, x_ax, y_ax,
    xdata_c, ydata_c, coldata_c, xdata_i, ydata_i, coldata_i,
        mode_fld_clean, local_rm_edges, c1, c2, gspos_idx):
    """
    Star's membership probabilities on cluster's photometric diagram.
    """

    pos = (
        2 * (gspos_idx // 3) + 2, 2 * (gspos_idx // 3) + 4,
        2 * (gspos_idx % 3), 2 * (gspos_idx % 3) + 2)

    ax = plt.subplot(gs[pos[0]:pos[1], pos[2]:pos[3]])
    # Set plot limits
    plt.xlim(x_min_cmd, x_max_cmd)
    plt.ylim(y_min_cmd, y_max_cmd)
    # Set axis labels
    plt.xlabel('$' + x_ax + '$', fontsize=xylabelsize)
    plt.ylabel('$' + y_ax + '$', fontsize=xylabelsize)
    # Add text box.
    text = r'$N_{{cmpl}}={},\, N_{{incmpl}}={}$'.format(
        sum(~np.isnan(xdata_c) & ~np.isnan(ydata_c)),
        sum(~np.isnan(xdata_i) & ~np.isnan(ydata_i)))
    ax.set_title(text, fontsize=titlesize)
    # Set minor ticks
    ax.minorticks_on()
    ax.tick_params(axis='both', which='major', labelsize=xytickssize)

    # Plot grid.
    if mode_fld_clean == 'local' and local_rm_edges is not None:
        # Use first magnitude and color.
        for x_ed in local_rm_edges[c1]:
            # vertical lines
            ax.axvline(x_ed, linestyle=':', lw=.8, color='k', zorder=1)
        for y_ed in local_rm_edges[c2]:
            # horizontal lines
            ax.axhline(y_ed, linestyle=':', lw=.8, color='k', zorder=1)
    else:
        ax.grid(b=True, which='major', color='gray', linestyle='--', lw=.5,
                zorder=1)

    # This reversed colormap means higher prob stars will look redder.
    cm = plt.cm.get_cmap('RdYlBu_r')
    if sum(~np.isnan(xdata_i) & ~np.isnan(ydata_i)) > 0:
        v_min_mp = round(min(np.nanmin(coldata_i), np.nanmin(coldata_c)), 2)
        v_max_mp = round(max(np.nanmax(coldata_i), np.nanmax(coldata_c)), 2)
    else:
        v_min_mp = round(np.nanmin(coldata_c), 2)
        v_max_mp = round(np.nanmax(coldata_c), 2)
    plot_colorbar = True if v_min_mp != v_max_mp else False
    # Plot stars used in the best fit process.
    sca = plt.scatter(
        xdata_c, ydata_c, marker='o', c=coldata_c, s=30, cmap=cm, lw=0.5,
        edgecolor='k', vmin=v_min_mp, vmax=v_max_mp, zorder=4)
    if sum(~np.isnan(xdata_i) & ~np.isnan(ydata_i)) > 0:
        plt.scatter(
            xdata_i, ydata_i, marker='^', c=coldata_i, s=30, cmap=cm, lw=0.5,
            edgecolor='k', zorder=6)
    # For plotting the colorbar
    trans = ax.transAxes + fig.transFigure.inverted()

    return plot_colorbar, sca, trans, v_min_mp, v_max_mp


def plot(N, *args):
    """
    Handle each plot separately.
    """

    plt_map = {
        0: [pl_mp_histo, 'MPs histogram'],
        1: [pl_chart_mps, 'frame with MPs coloring']
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

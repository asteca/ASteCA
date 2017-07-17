
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib.offsetbox as offsetbox
import numpy as np


def pl_mp_histo(
        gs, n_memb_da, memb_prob_avrg_sort, flag_decont_skip, cl_reg_fit,
        min_prob, mode_red_memb, local_bin):
    '''
    Histogram for the distribution of membership probabilities from the
    decontamination algorithm.
    '''
    # Only attempt to plot if the DA was applied.
    if flag_decont_skip is False:
        # Reduced membership.
        ax = plt.subplot(gs[0:2, 6:8])
        plt.xlim(0., 1.)
        plt.xlabel('MP (membership probability)', fontsize=12)
        plt.ylabel('N (normalized)', fontsize=12)
        ax.minorticks_on()
        ax.grid(b=True, which='major', color='gray', linestyle='--', lw=1,
                zorder=1)
        prob_data = [star[7] for star in memb_prob_avrg_sort]
        # Histogram of the data.
        n_bins = int((max(prob_data) - min(prob_data)) / 0.025)
        if n_bins > 0:
            # Normalized histogram.
            weights = np.ones_like(prob_data) / len(prob_data)
            n, bins, patches = plt.hist(prob_data, n_bins, weights=weights,
                                        normed=0)
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
        # Add text box.
        if mode_red_memb == 'mag':
            str_pm = ['mag', '\leq', 'mag']
        else:
            str_pm = ['MP', '\geq', 'prob']
        if mode_red_memb == 'local':
            str_pm.append(mode_red_memb + ';\,' + local_bin)
        else:
            str_pm.append(mode_red_memb.replace('_', '\_'))
        text1 = r'$n_{{memb-DA}}={}\,(MP \geq 0.5)$'.format(n_memb_da)
        text2 = r'${}_{{min}}={:.2f}\,({})$'.format(str_pm[2], min_prob,
                                                    str_pm[3])
        text3 = r'$N_{{fit}}={} \, ({} {} {}_{{min}})$'.format(
            len(cl_reg_fit), str_pm[0], str_pm[1], str_pm[2])
        text = text1 + '\n' + text2 + '\n' + text3
        # Plot minimum probability line.
        plt.axvline(x=min_prob, linestyle='--', color='green', lw=2.5,
                    zorder=3)
        ob = offsetbox.AnchoredText(text, loc=2, prop=dict(size=10))
        ob.patch.set(boxstyle='square,pad=0.05', alpha=0.85)
        ax.add_artist(ob)
        # Avoid showing the value 0.0 in the y axis.
        plt.ylim(0.001, plt.ylim()[1])


def pl_chart_mps(gs, fig, x_name, y_name, coord, x_zmin, x_zmax, y_zmin,
                 y_zmax, clust_cent, clust_rad, field_dens, flag_decont_skip,
                 v_min_mp, v_max_mp, chart_fit_inv, chart_no_fit_inv,
                 out_clust_rad, mode_red_memb, local_bin):
    '''
    Finding chart of cluster region with decontamination algorithm
    applied and colors assigned according to the probabilities obtained.
    '''
    ax = plt.subplot(gs[0:2, 8:10])
    # Set plot limits, Use 'zoom' x,y ranges.
    plt.xlim(x_zmin, x_zmax)
    plt.ylim(y_zmin, y_zmax)
    # If RA is used, invert axis.
    if coord == 'deg':
        ax.invert_xaxis()
    # Set axis labels
    plt.xlabel('{} ({})'.format(x_name, coord), fontsize=12)
    plt.ylabel('{} ({})'.format(y_name, coord), fontsize=12)
    # Set minor ticks
    ax.minorticks_on()
    # Radius
    circle = plt.Circle((clust_cent[0], clust_cent[1]), clust_rad, color='red',
                        fill=False)
    fig.gca().add_artist(circle)
    ob = offsetbox.AnchoredText('Cluster region', loc=1, prop=dict(size=12))
    ob.patch.set(boxstyle='square,pad=-0.2', alpha=0.85)
    ax.add_artist(ob)
    # If DA was skipped, print info on 'local' method here.
    if flag_decont_skip and mode_red_memb == 'local':
        text = r'$({})$'.format(mode_red_memb + ';\,' + local_bin)
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
    # Star sizes for dense and not dense regions.
    st_size = 20 if field_dens > 0.005 else 35
    # Plot stars *not* used in the best fit process.
    plt.scatter(
        chart_no_fit_inv[0], chart_no_fit_inv[1], marker='o',
        c=col_select_no_fit, s=st_size, edgecolors='black', cmap=cm,
        alpha=0.5, lw=0.5, vmin=v_min_mp, vmax=v_max_mp)
    # Add line to stars not used in the best bit process.
    plt.scatter(chart_no_fit_inv[0], chart_no_fit_inv[1], marker='_',
                c='k', lw=0.5, alpha=0.5)
    # Plot stars selected to be used in the best bit process.
    plt.scatter(chart_fit_inv[0], chart_fit_inv[1], marker='o',
                c=col_select_fit, s=st_size, edgecolors='black', cmap=cm,
                lw=0.5, vmin=v_min_mp, vmax=v_max_mp, zorder=4)
    # Plot stars outside the cluster's radius as empty circles.
    plt.scatter(out_clust_rad[0], out_clust_rad[1], marker='o',
                s=st_size, edgecolors='black', facecolors='none', lw=0.5)


def pl_mps_phot_diag(gs, fig, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd,
                     x_ax, y_ax, v_min_mp, v_max_mp, diag_fit_inv,
                     diag_no_fit_inv, err_bar, mode_red_memb, bin_edges):
    '''
    Star's membership probabilities on cluster's photometric diagram.
    '''
    x_val, mag_y, x_err, y_err = err_bar
    ax = plt.subplot(gs[0:2, 10:12])
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
    if mode_red_memb == 'local' and bin_edges != 0.:
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
                c=col_select_no_fit, s=35, cmap=cm, lw=0.5, edgecolor='k',
                alpha=0.5, vmin=v_min_mp, vmax=v_max_mp, zorder=2)
    # Draw horizontal line over stars discarded.
    plt.scatter(diag_no_fit_inv[1][0], diag_no_fit_inv[0][0],
                marker='_', c='k', lw=0.5, alpha=0.5, zorder=3)
    # Plot stars used in the best fit process.
    sca = plt.scatter(diag_fit_inv[1][0], diag_fit_inv[0][0], marker='o',
                      c=col_select_fit, s=40, cmap=cm, lw=0.5, edgecolor='k',
                      vmin=v_min_mp, vmax=v_max_mp, zorder=4)
    # If list is not empty, plot error bars at several values.
    if x_val:
        plt.errorbar(x_val, mag_y, yerr=y_err, xerr=x_err, fmt='k.', lw=0.8,
                     ms=0., zorder=4)
    # For plotting the colorbar (see bottom of make_plots file).
    trans = ax.transAxes + fig.transFigure.inverted()

    return sca, trans


def plot(N, *args):
    '''
    Handle each plot separately.
    '''

    plt_map = {
        0: [pl_mp_histo, 'MPs histogram'],
        1: [pl_chart_mps, 'frame with MPs coloring']
    }

    fxn = plt_map.get(N, None)[0]
    if fxn is None:
        raise ValueError("  ERROR: there is no plot {}.".format(N))

    try:
        fxn(*args)
    except:
        # import traceback
        # print traceback.format_exc()
        print("  WARNING: error when plotting {}.".format(plt_map.get(N)[1]))

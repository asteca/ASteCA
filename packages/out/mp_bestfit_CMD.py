
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.offsetbox as offsetbox
from matplotlib.colors import LinearSegmentedColormap, LogNorm, ListedColormap
import matplotlib.patches as mpatches


def pl_mps_phot_diag(
    gs, gs_y1, gs_y2, fig, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd,
    x_ax, y_ax, v_min_mp, v_max_mp, obs_x, obs_y, obs_MPs, err_bar,
    cl_sz_pt, hess_xedges, hess_yedges, x_isoch, y_isoch, phot_Nsigma,
        lkl_method, redVect):
    """
    Star's membership probabilities on cluster's photometric diagram.
    """
    ax = plt.subplot(gs[gs_y1:gs_y2, 0:2])
    # Set axis labels
    plt.xlabel('$' + x_ax + '$')
    plt.ylabel('$' + y_ax + '$')
    # Add text box.
    if gs_y1 == 0:
        text = '$N_{{fit}}={}$'.format(len(obs_MPs))
        ob = offsetbox.AnchoredText(text, loc=3)
        ob.patch.set(boxstyle='square,pad=-0.2', alpha=0.85)
        ax.add_artist(ob)

    if gs_y1 == 0:
        txt = r" + $\mathcal{N}(\mu,\sigma^2)$" if phot_Nsigma else ""
        ax.set_title("Observed" + txt)

    # Plot grid.
    gls, glw, gc = plt.rcParams['grid.linestyle'],\
        plt.rcParams['grid.linewidth'], plt.rcParams['grid.color']
    for x_ed in hess_xedges:
        # vertical lines
        ax.axvline(x_ed, linestyle=gls, lw=glw, color=gc, zorder=1)
    for y_ed in hess_yedges:
        # horizontal lines
        ax.axhline(y_ed, linestyle=gls, lw=glw, color=gc, zorder=1)

    # This reversed colormap means higher prob stars will look redder.
    rmap = plt.cm.get_cmap('RdYlBu_r')
    # If the 'tolstoy' method was used AND the stars have a range of colors.
    # The 'dolphin / tremmel' likelihoods do not use MPs in the fit, so it's
    # confusing to color stars as if they did.
    if (v_min_mp != v_max_mp) and lkl_method == 'tolstoy':
        col_select_fit = obs_MPs
        plot_colorbar = True
    else:
        col_select_fit = '#519ddb'
        plot_colorbar = False
    # Plot stars used in the best fit process.
    sca = plt.scatter(
        obs_x, obs_y, marker='o', c=col_select_fit, s=cl_sz_pt, cmap=rmap,
        lw=0.3, edgecolor='k', vmin=v_min_mp, vmax=v_max_mp, zorder=4)

    # Plot sigma region
    if phot_Nsigma:
        cGreys = plt.cm.get_cmap('Greys', 100)
        cmap = ListedColormap(cGreys(range(65)))
        # Extend one bin upwards and to the left
        ybin = abs(hess_yedges[1] - hess_yedges[0])
        hess_yedges = [hess_yedges[0] - ybin] + list(hess_yedges)
        xbin = abs(hess_xedges[1] - hess_xedges[0])
        hess_xedges = [hess_xedges[0] - xbin] + list(hess_xedges)
        plt.hist2d(*phot_Nsigma, bins=(
            hess_xedges, hess_yedges), cmap=cmap, norm=LogNorm(), zorder=-1)
    # Plot isochrone.
    plt.plot(x_isoch, y_isoch, c='k', lw=1., zorder=6)

    plt.xlim(x_min_cmd, x_max_cmd)
    plt.ylim(y_min_cmd, y_max_cmd)

    # Reddening vector
    x0, y0, dx, dy = redVect
    prop = dict(arrowstyle="-|>,head_width=0.4,head_length=0.8",
                color='indianred', shrinkA=0, shrinkB=0, lw=2)
    plt.annotate("", xy=(x0 + dx, y0 + dy), xytext=(x0, y0), arrowprops=prop)

    # If list is not empty, plot error bars at several values. The
    # prep_plots.error_bars() is not able to handle the color-color diagram.
    x_val, mag_y, xy_err = err_bar
    if x_val:
        xye_i = {
            '0': (mag_y, 0, 1), '2': (mag_y, 0, 2),
            '4': (np.linspace(min(obs_y), max(obs_y), len(x_val)), 1, 2)}
        plt.errorbar(
            x_val, xye_i[str(gs_y1)][0], yerr=xy_err[xye_i[str(gs_y1)][1]],
            xerr=xy_err[xye_i[str(gs_y1)][2]],
            fmt='k.', lw=0.8, ms=0., zorder=4)
    # For plotting the colorbar (see bottom of make_D_plot file).
    trans = ax.transAxes + fig.transFigure.inverted()

    return plot_colorbar, sca, trans


def pl_hess_diag(
    gs, gs_y1, gs_y2, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd, x_ax, y_ax,
    lkl_method, hess_xedges, hess_yedges, hess_x, hess_y, HD, x_isoch,
        y_isoch):
    """
    Hess diagram of observed minus best match synthetic cluster.
    """
    ax = plt.subplot(gs[gs_y1:gs_y2, 2:4])
    # Set plot limits
    plt.xlim(x_min_cmd, x_max_cmd)
    plt.ylim(y_min_cmd, y_max_cmd)
    # Set axis labels
    plt.xlabel('$' + x_ax + '$')
    if gs_y1 == 0:
        ax.set_title("Hess diagram (observed - synthetic)")

    gls, glw, gc = plt.rcParams['grid.linestyle'],\
        plt.rcParams['grid.linewidth'], plt.rcParams['grid.color']
    for x_ed in hess_xedges:
        # vertical lines
        ax.axvline(x_ed, linestyle=gls, lw=glw, color=gc, zorder=3)
    for y_ed in hess_yedges:
        # horizontal lines
        ax.axhline(y_ed, linestyle=gls, lw=glw, color=gc, zorder=3)
    if HD.any():
        # Add text box.
        if HD.min() < 0:
            plt.scatter(-100., -100., marker='s', lw=0., s=60, c='#0B02F8',
                        label='{}'.format(int(HD.min())))
        if HD.max() > 0:
            plt.scatter(-100., -100., marker='s', lw=0., s=60, c='#FB0605',
                        label='{}'.format(int(HD.max())))
        # Define custom colorbar.
        if HD.min() == 0:
            cmap = LinearSegmentedColormap.from_list(
                'mycmap', [(0, 'white'), (1, 'red')])
        else:
            # Zero point for empty bins which should be colored in white.
            zero_pt = (0. - HD.min()) / float(HD.max() - HD.min())
            N = 256.
            zero_pt0 = np.floor(zero_pt * (N - 1)) / (N - 1)
            zero_pt1 = np.ceil(zero_pt * (N - 1)) / (N - 1)
            cmap = LinearSegmentedColormap.from_list(
                'mycmap', [(0, 'blue'), (zero_pt0, 'white'),
                           (zero_pt1, 'white'), (1, 'red')], N=N)
        ax.pcolormesh(hess_x, hess_y, HD, cmap=cmap, vmin=HD.min(),
                      vmax=HD.max(), zorder=1)
        # Legend.
        handles, labels = ax.get_legend_handles_labels()
        leg = ax.legend(
            handles, labels, loc='lower left', scatterpoints=1, ncol=2,
            columnspacing=.2, handletextpad=-.3)
        leg.get_frame().set_alpha(0.7)

    plt.plot(x_isoch, y_isoch, c='k', lw=1., zorder=6)


def pl_bf_synth_cl(
    gs, gs_y1, gs_y2, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd, x_ax, y_ax,
    hess_xedges, hess_yedges, x_synth, y_synth, sy_sz_pt, binar_idx, best_sol,
    x_isoch, y_isoch, IMF_name, DR_percentage, alpha, gamma, lkl_method,
    bin_method, evol_track, D3_sol, MassT_dist_vals, binar_dist_vals,
        p_err):
    """
    Best fit synthetic cluster obtained.
    """
    ax = plt.subplot(gs[gs_y1:gs_y2, 4:6])
    # Set plot limits
    plt.xlim(x_min_cmd, x_max_cmd)
    plt.ylim(y_min_cmd, y_max_cmd)
    # Set axis labels
    plt.xlabel('$' + x_ax + '$')

    if gs_y1 == 0:
        ax.set_title("Synthetic ({} solution)".format(D3_sol))
        # Add text box
        text = r'$({};\,{})$'.format(lkl_method, bin_method)
        ob = offsetbox.AnchoredText(text, pad=.2, loc=1)
        ob.patch.set(alpha=0.85)
        ax.add_artist(ob)

    gls, glw, gc = plt.rcParams['grid.linestyle'],\
        plt.rcParams['grid.linewidth'], plt.rcParams['grid.color']
    for x_ed in hess_xedges:
        # vertical lines
        ax.axvline(x_ed, linestyle=gls, lw=glw, color=gc, zorder=1)
    for y_ed in hess_yedges:
        # horizontal lines
        ax.axhline(y_ed, linestyle=gls, lw=glw, color=gc, zorder=1)

    # Single systems
    plt.scatter(
        x_synth[~binar_idx], y_synth[~binar_idx], marker='o', s=sy_sz_pt,
        c='#519ddb', lw=0.3, edgecolor='k', zorder=2)
    # Binary systems
    plt.scatter(
        x_synth[binar_idx], y_synth[binar_idx], marker='o', s=sy_sz_pt,
        c='#F34C4C', lw=0.3, edgecolor='k', zorder=3)
    # Plot isochrone.
    plt.plot(x_isoch, y_isoch, c='k', lw=1., zorder=6)
    if gs_y1 == 0:
        # Add text box
        text1 = '$N_{{synth}} = {}$'.format(len(x_synth))
        ob = offsetbox.AnchoredText(text1, pad=.2, loc=3)
        ob.patch.set(alpha=0.85)
        ax.add_artist(ob)

        # Add text box to the right of the synthetic cluster.
        ax_t = plt.subplot(gs[gs_y1:gs_y2, 6:7])
        ax_t.axis('off')  # Remove axis from frame.
        t1 = r'Tracks: {}'.format(evol_track)
        t2 = r'IMF: {}'.format(IMF_name.replace('_', ' ').title())

        t3 = r'$M\,(M_{{\odot}}) \hspace{{1.}}=\;{:.0f}\pm {:.0f}$'.format(
            MassT_dist_vals[D3_sol + '_sol'], MassT_dist_vals['errors'][2])

        t4 = r'$z \hspace{{3.9}}=\;{:.4f}\pm {:.4f}$'.format(
            best_sol[0], p_err[0][2])
        t5 = r'$\log(age) \hspace{{0.17}}=\;{:.3f}\pm {:.3f}$'.format(
            best_sol[1], p_err[1][2])

        t6 = r"$\alpha \hspace{{3.85}}=\;{:.3f}$".format(alpha)
        t7 = r'$\beta \hspace{{3.85}}=\;{:.3f}\pm {:.3f}$'.format(
            best_sol[2], p_err[2][2])
        t8 = r"$\gamma \hspace{{3.85}}=\;{}$".format(gamma)

        # If 'beta' was not fitted, inherit the 'nan' uncertainty
        if np.isnan(p_err[2][2]):
            e_bfr = np.nan
        else:
            e_bfr = binar_dist_vals['errors'][2]
        t9 = r'$b_{{frac}} \hspace{{2.4}}=\;{:.3f}\pm {:.3f}$'.format(
            binar_dist_vals[D3_sol + '_sol'], e_bfr)

        t10 = r'$A_{{V}} \hspace{{3.2}}=\;{:.3f}\pm {:.3f}$'.format(
            best_sol[3], p_err[3][2])
        t11 = r'$DR \hspace{{2.9}} =\;{:.3f}\pm {:.3f}$'.format(  # \,({})
            best_sol[4], p_err[4][2])  # DR_percentage
        t12 = r'$R_{{V}} \hspace{{3.2}}=\;{:.3f}\pm {:.3f}$'.format(
            best_sol[5], p_err[5][2])
        t13 = r'$(m-M)_{{0}} =\;{:.2f} \pm {:.2f}$'.format(
            best_sol[6], p_err[6][2])

        text = t1 + '\n' + t2 + '\n\n' + t3 + '\n' + t4 + '\n' + t5 + '\n' +\
            t6 + '\n' + t7 + '\n' + t8 + '\n' + t9 + '\n\n' + t10 + '\n' +\
            t11 + '\n' + t12 + '\n' + t13
        ob = offsetbox.AnchoredText(text, pad=1, loc=6, borderpad=-5)
        ob.patch.set(alpha=0.85)
        ax_t.add_artist(ob)


def plot(N, *args):
    """
    Handle each plot separately.
    """
    plt_map = {
        0: [pl_hess_diag, 'Hess diagram'],
        1: [pl_bf_synth_cl, 'synthetic cluster']
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


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.offsetbox as offsetbox
from matplotlib.ticker import MultipleLocator
from matplotlib.colors import LinearSegmentedColormap


def pl_mps_phot_diag(gs, fig, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd,
                     x_ax, y_ax, v_min_mp, v_max_mp, diag_fit_inv,
                     err_bar, hess_data, shift_isoch):
    '''
    Star's membership probabilities on cluster's photometric diagram.
    '''
    x_val, mag_y, x_err, y_err = err_bar
    ax = plt.subplot(gs[0:2, 0:2])
    # Set plot limits
    plt.xlim(x_min_cmd, x_max_cmd)
    plt.ylim(y_min_cmd, y_max_cmd)
    # Set axis labels
    plt.xlabel('$' + x_ax + '$', fontsize=18)
    plt.ylabel('$' + y_ax + '$', fontsize=18)
    # Add text box.
    text = '$N_{{fit}}={}$'.format(len(diag_fit_inv[2]))
    ob = offsetbox.AnchoredText(text, loc=4, prop=dict(size=14))
    ob.patch.set(boxstyle='square,pad=-0.2', alpha=0.85)
    ax.add_artist(ob)
    # Set minor ticks
    ax.minorticks_on()
    ax.xaxis.set_major_locator(MultipleLocator(1.0))
    ax.set_title("Observed", fontsize=10)
    # TODO this assumes that the first photom dimen is the main mag, and the
    # second the color.
    # Plot grid.
    for x_ed in hess_data['hess_edges'][1]:
        # vertical lines
        ax.axvline(x_ed, linestyle=':', lw=.8, color='k', zorder=1)
    for y_ed in hess_data['hess_edges'][0]:
        # horizontal lines
            ax.axhline(y_ed, linestyle=':', lw=.8, color='k', zorder=1)
    # This reversed colormap means higher prob stars will look redder.
    cm = plt.cm.get_cmap('RdYlBu_r')
    # If the 'tolstoy method was used AND the stars have a range of colors.
    # Currently the 'dolphin' likelihood does not use MPs in the fit, so it's
    # confusing to color stars is if it did.
    if v_min_mp != v_max_mp:
        col_select_fit, isoch_col = diag_fit_inv[2], 'g'
    else:
        col_select_fit, isoch_col = '#4682b4', 'r'
    # Plot stars used in the best fit process.
    sca = plt.scatter(diag_fit_inv[1][0], diag_fit_inv[0][0], marker='o',
                      c=col_select_fit, s=40, cmap=cm, lw=0.5, edgecolor='k',
                      vmin=v_min_mp, vmax=v_max_mp, zorder=4)
    # TODO using main magnitude and first color.
    # Plot isochrone.
    if shift_isoch:
        plt.plot(shift_isoch[1], shift_isoch[0], isoch_col, lw=1.2, zorder=6)
    # If list is not empty, plot error bars at several values.
    if x_val:
        plt.errorbar(x_val, mag_y, yerr=y_err, xerr=x_err, fmt='k.', lw=0.8,
                     ms=0., zorder=4)
    # For plotting the colorbar (see bottom of make_D_plot file).
    trans = ax.transAxes + fig.transFigure.inverted()

    return sca, trans


def pl_hess_diag(
    gs, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd, x_ax, y_ax, lkl_method,
        hess_data):
    """
    Hess diagram of observed minus best match synthetic cluster.
    """
    ax = plt.subplot(gs[0:2, 2:4])
    # Set plot limits
    plt.xlim(x_min_cmd, x_max_cmd)
    plt.ylim(y_min_cmd, y_max_cmd)
    # Set axis labels
    plt.xlabel('$' + x_ax + '$', fontsize=18)
    plt.ylabel('$' + y_ax + '$', fontsize=18)
    # Set minor ticks
    ax.minorticks_on()
    ax.xaxis.set_major_locator(MultipleLocator(1.0))
    ax.set_title("Hess diagram (observed - synthetic)", fontsize=10)
    # TODO this assumes that the first photom dimen is the main mag, and the
    # second the color.
    for x_ed in hess_data['hess_edges'][1]:
        # vertical lines
        ax.axvline(x_ed, linestyle=':', lw=.8, color='k', zorder=1)
    for y_ed in hess_data['hess_edges'][0]:
        # horizontal lines
        ax.axhline(y_ed, linestyle=':', lw=.8, color='k', zorder=1)
    # Grid for pcolormesh, using first two dimensions
    Y, X = np.meshgrid(hess_data['hess_edges'][0], hess_data['hess_edges'][1])
    if hess_data['hess_diag'].size:
        HD = np.rot90(hess_data['hess_diag'])
        HD = np.flipud(HD)
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
                'mycmap', [(0, 'blue'), (zero_pt0, 'white'), (zero_pt1,
                           'white'), (1, 'red')], N=N)
        ax.pcolormesh(X, Y, HD, cmap=cmap, vmin=HD.min(), vmax=HD.max(),
                      zorder=1)
        # Legend.
        handles, labels = ax.get_legend_handles_labels()
        leg = ax.legend(
            handles, labels, loc='lower right', scatterpoints=1, ncol=2,
            columnspacing=.2, handletextpad=-.3, fontsize=10)
        leg.get_frame().set_alpha(0.7)


def pl_bf_synth_cl(
    gs, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd, x_ax, y_ax, synth_clst,
        hess_data, IMF_name, R_V, cp_r, cp_e, shift_isoch, lkl_method,
        bin_method, cmd_evol_tracks, evol_track):
    '''
    Best fit synthetic cluster obtained.
    '''
    ax = plt.subplot(gs[0:2, 4:6])
    # Set plot limits
    plt.xlim(x_min_cmd, x_max_cmd)
    plt.ylim(y_min_cmd, y_max_cmd)
    # Set axis labels
    plt.xlabel('$' + x_ax + '$', fontsize=18)
    plt.ylabel('$' + y_ax + '$', fontsize=18)
    # Set minor ticks
    ax.minorticks_on()
    ax.xaxis.set_major_locator(MultipleLocator(1.0))
    ax.set_title("Synthetic (best match)", fontsize=10)
    # Plot grid.
    # TODO this assumes that the first photom dimen is the main mag, and the
    # second the color.
    for x_ed in hess_data['hess_edges'][1]:
        # vertical lines
        ax.axvline(x_ed, linestyle=':', lw=.8, color='k', zorder=1)
    for y_ed in hess_data['hess_edges'][0]:
        # horizontal lines
        ax.axhline(y_ed, linestyle=':', lw=.8, color='k', zorder=1)
    # Add text box
    text = '$({};\,{})$'.format(lkl_method, bin_method)
    ob = offsetbox.AnchoredText(text, pad=.2, loc=1, prop=dict(size=12))
    ob.patch.set(alpha=0.85)
    ax.add_artist(ob)
    if synth_clst:
        # TODO using main magnitude and first color
        # Plot synthetic cluster.
        plt.scatter(synth_clst[0][0][1], synth_clst[0][0][0], marker='o', s=40,
                    c='#4682b4', lw=0.5, edgecolor='k', zorder=4)
        text1 = '$N_{{synth}} = {}$'.format(len(synth_clst[0][0][0]))
        # Plot isochrone.
        plt.plot(shift_isoch[1], shift_isoch[0], 'r', lw=1.2, zorder=6)
    else:
        text1 = '$N_{{synth}} = {}$'.format(0)
    # Add text box
    ob = offsetbox.AnchoredText(text1, pad=.2, loc=3, prop=dict(size=14))
    ob.patch.set(alpha=0.85)
    ax.add_artist(ob)

    # Add text box to the right of the synthetic cluster.
    ax_t = plt.subplot(gs[0:2, 6:7])
    ax_t.axis('off')  # Remove axis from frame.
    # Map isochrones set selection to proper name.
    iso_print = cmd_evol_tracks[evol_track][1]
    t1 = r'$Synthetic\;cluster\;parameters$' + '\n' + \
        r'$[Tracks:\;{}]$'.format(iso_print.replace(' ', '\;'))
    t2 = r'$IMF \hspace{{3.}}:\;{}$'.format(
        IMF_name.replace('_', '\;').title())
    t3 = r'$R_{{V}} \hspace{{3.2}}=\;{}$'.format(R_V)
    t4 = r'$z \hspace{{3.9}}=\;{}\pm {}$'.format(cp_r[0], cp_e[0])
    t5 = r'$\log(age) \hspace{{0.17}}=\;{}\pm {}$'.format(cp_r[1], cp_e[1])
    t6 = r'$E_{{(B-V)}} \hspace{{1.35}}=\;{}\pm {}$'.format(cp_r[2], cp_e[2])
    t7 = r'$(m-M)_o=\;{} \pm {}$'.format(cp_r[3], cp_e[3])
    t8 = r'$M\,(M_{{\odot}}) \hspace{{1.07}} =\;{}\pm {}$'.format(
        cp_r[4], cp_e[4])
    t9 = r'$b_{{frac}} \hspace{{2.37}}=\;{}\pm {}$'.format(cp_r[5], cp_e[5])
    text = t1 + '\n\n' + t2 + '\n' + t3 + '\n' + t4 + '\n' + t5 + '\n' + t6 +\
        '\n' + t7 + '\n' + t8 + '\n' + t9
    ob = offsetbox.AnchoredText(
        text, pad=1, loc=6, borderpad=-5, prop=dict(size=11))
    ob.patch.set(alpha=0.85)
    ax_t.add_artist(ob)


def plot(N, *args):
    '''
    Handle each plot separately.
    '''
    plt_map = {
        0: [pl_hess_diag, 'Hess diagram'],
        1: [pl_bf_synth_cl, 'synthetic cluster']
    }

    fxn = plt_map.get(N, None)[0]
    if fxn is None:
        raise ValueError("  ERROR: there is no plot {}.".format(N))

    try:
        fxn(*args)
    except:
        import traceback
        print traceback.format_exc()
        print("  WARNING: error when plotting {}.".format(plt_map.get(N)[1]))

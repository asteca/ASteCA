
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.offsetbox as offsetbox
from . import prep_plots


def starsPlot(boundary, x_data, y_data):
    """
    Plot accepted/rejected stars outside/inside the cluster region.
    """
    if boundary == 'rjct':
        if len(y_data) > 0:
            # Only attempt to plot if any star is stored in the list.
            plt.scatter(
                x_data, y_data, marker='x', c='teal', s=35, lw=.5, zorder=1)
    if boundary == 'accpt_in':
        if len(y_data) > 0:
            plt.scatter(
                x_data, y_data, marker='o', c='r', s=20, zorder=3,
                lw=0.3, edgecolor='k', label=r'$r \leq r_{cl}}$')
    if boundary == 'accpt_out':
        if len(y_data) > 0:
            plt.scatter(
                x_data, y_data, marker='o', c='b', s=15, zorder=2,
                lw=0.1, edgecolor='k', label=r'$field$')


def pl_phot_err(
    gs, colors, filters, id_kinem, mags, em_float, cl_region_c,
        cl_region_rjct_c, stars_out_c, stars_out_rjct_c, err_bar_all):
    '''
    Photometric + kinematic error rejection.
    '''

    # Main magnitude (x) data for accepted/rejected stars.
    mmag_out_acpt, mmag_out_rjct, mmag_in_acpt, mmag_in_rjct =\
        np.array([]), np.array([]), np.array([]), np.array([])
    if stars_out_c:
        mmag_out_acpt = np.array(list(zip(*list(zip(*stars_out_c))[3]))[0])
    if stars_out_rjct_c:
        mmag_out_rjct = np.array(list(zip(*list(zip(
            *stars_out_rjct_c))[3]))[0])
    if cl_region_c:
        mmag_in_acpt = np.array(list(zip(*list(zip(*cl_region_c))[3]))[0])
    if cl_region_rjct_c:
        mmag_in_rjct = np.array(list(zip(*list(zip(*cl_region_rjct_c))[3]))[0])

    # Define parameters for main magnitude error plot.
    y_ax, x_ax = prep_plots.ax_names(filters[0], filters[0], 'mag')
    err_plot = [[x_ax, y_ax, 4, 0]]
    # For all defined colors.
    for i, _ in enumerate(colors):
        y_ax, _ = prep_plots.ax_names(colors[i], filters[0], 'mag')
        err_plot.append([x_ax, y_ax, 6, i])

    pd_Plx, pd_PMRA, pd_PMDE, pd_RV = id_kinem[0], id_kinem[2], id_kinem[4],\
        id_kinem[6]
    # For the kinematic data
    if pd_Plx != 'n':
        err_plot.append([x_ax, "Plx", 8, 0])
    if pd_PMRA != 'n':
        err_plot.append([x_ax, "PMra", 8, 1])
    if pd_PMDE != 'n':
        err_plot.append([x_ax, "PMde", 8, 2])
    if pd_RV != 'n':
        err_plot.append([x_ax, "RV", 8, 3])

    # Set plot limits
    x_min, x_max = min(mags[0]) - 0.5, max(mags[0]) + 0.5
    for i, pl in enumerate(err_plot):
        x_ax, y_ax, j, k = pl

        ax = plt.subplot(
            gs[i // 2: (i // 2) + 1, 3 * (i % 2):3 * (i % 2) + 3])
        ax.grid(b=True, which='major', color='gray', linestyle='--', lw=.5,
                zorder=1)
        plt.xlim(x_min, x_max)
        # Set axis labels
        plt.xlabel(r'$' + x_ax + r'$', fontsize=14)
        plt.ylabel(r'$\sigma_{{{}}}$'.format(y_ax), fontsize=14)
        ax.set_facecolor('#EFF0F1')
        # Set minor ticks
        ax.minorticks_on()

        # Rejected stars outside the cluster region
        if any(mmag_out_rjct) and any(stars_out_rjct_c):
            starsPlot('rjct', mmag_out_rjct,
                      list(zip(*list(zip(*stars_out_rjct_c))[j]))[k])
        if any(mmag_in_rjct) and any(cl_region_rjct_c):
            # Rejected stars inside the cluster region
            starsPlot('rjct', mmag_in_rjct,
                      list(zip(*list(zip(*cl_region_rjct_c))[j]))[k])
        if any(mmag_in_acpt) and any(cl_region_c):
            # Accepted stars inside the cluster region.
            starsPlot('accpt_in', mmag_in_acpt,
                      list(zip(*list(zip(*cl_region_c))[j]))[k])
        if any(mmag_out_acpt) and any(stars_out_c):
            # Accepted stars outside the cluster region.
            starsPlot('accpt_out', mmag_out_acpt,
                      list(zip(*list(zip(*stars_out_c))[j]))[k])

        if j == 4:
            # Plot legend in the main magnitude plot.
            leg = plt.legend(fancybox=True, loc='upper left', scatterpoints=1,
                             fontsize=16, markerscale=2.5, prop={'size': 13})
            # Set the alpha value of the legend.
            leg.get_frame().set_alpha(0.7)
            # Max error cut
            max_cut_y = em_float[0]
            ax.hlines(y=max_cut_y, xmin=x_min, xmax=x_max, color='k',
                      linestyles='dashed', zorder=4)
            ob = offsetbox.AnchoredText(
                r"$max={}$ mag".format(em_float[0]), loc=1, prop=dict(size=11))
            ob.patch.set(alpha=0.7)
            ax.add_artist(ob)
            # Plot error curve
            plt.plot(err_bar_all[1], err_bar_all[2][0], color='yellow',
                     ls='--', lw=2, zorder=5)

        elif j == 6:
            max_cut_y = em_float[1 + k]
            ax.hlines(y=max_cut_y, xmin=x_min, xmax=x_max, color='k',
                      linestyles='dashed', zorder=4)
            ob = offsetbox.AnchoredText(
                r"$max={}$ mag".format(em_float[1 + k]), loc=2,
                prop=dict(size=11))
            ob.patch.set(alpha=0.7)
            ax.add_artist(ob)
            plt.plot(err_bar_all[1], err_bar_all[2][k + 1], color='yellow',
                     ls='--', lw=2, zorder=5)
        else:
            unit = '[mas]' if k == 0 else '[mas/yr]'
            # Use single PM max error value defined, for PMde
            k = 1 if k == 2 else 1
            max_cut_y = em_float[-(3 - k)]
            ax.hlines(y=max_cut_y, xmin=x_min, xmax=x_max, color='k',
                      linestyles='dashed', zorder=4)
            ob = offsetbox.AnchoredText(
                r"$max={}$ {}".format(em_float[-(3 - k)], unit), loc=2,
                prop=dict(size=11))
            ob.patch.set(alpha=0.7)
            ax.add_artist(ob)
        # Maximum error limit of 1.
        plt.ylim(-0.0024, min(plt.ylim()[1], 2. * max_cut_y, 1.))


def plot(N, *args):
    '''
    Handle each plot separately.
    '''

    plt_map = {
        0: [pl_phot_err, 'error rejection function']
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

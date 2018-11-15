
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
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
                x_data, y_data, marker='x', c='teal', s=15, lw=.5, zorder=1)
    if boundary == 'accpt_in':
        if len(y_data) > 0:
            plt.scatter(
                x_data, y_data, marker='o', c='r', s=10, zorder=3,
                lw=0.3, edgecolor='k', label='$r \leq r_{cl}}$')
    if boundary == 'accpt_out':
        if len(y_data) > 0:
            plt.scatter(
                x_data, y_data, marker='o', c='b', s=5, zorder=2,
                lw=0.1, edgecolor='k', label='$r > r_{cl}$')


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

    pd_Plx, pd_PMRA, pd_RV = id_kinem[0], id_kinem[2], id_kinem[6]
    # Define first row depending on whether kinematic data was defined.
    gs0 = 0
    if np.array([_ == 'n' for _ in (pd_Plx, pd_PMRA, pd_RV)]).all():
        gs0 = 1

    # Define parameters for main magnitude error plot.
    y_ax, x_ax = prep_plots.ax_names(filters[0], filters[0], 'mag')
    err_plot = [[gs[gs0, 0:2], x_ax, y_ax, 4, 0]]
    # For up to two defined colors.
    for i, _ in enumerate(colors[:2]):
        y_ax, _ = prep_plots.ax_names(colors[i], filters[0], 'mag')
        err_plot.append([gs[gs0, (2 * i + 2):(2 * i + 4)], x_ax, y_ax, 6, i])

    # For the kinematic data
    if pd_Plx != 'n':
        err_plot.append([gs[1, 0:2], x_ax, "Plx", 8, 0])
    if pd_PMRA != 'n':
        err_plot.append([gs[1, 2:4], x_ax, "PMxy", 8, 1])
    if pd_RV != 'n':
        err_plot.append([gs[1, 4:6], x_ax, "RV", 8, 3])

    # Set plot limits
    x_min, x_max = min(mags[0]) - 0.5, max(mags[0]) + 0.5
    for pl in err_plot:
        gs_pos, x_ax, y_ax, j, k = pl

        ax = plt.subplot(gs_pos)
        plt.xlim(x_min, x_max)
        # Set axis labels
        plt.xlabel('$' + x_ax + '$', fontsize=12)
        plt.ylabel('$\sigma_{' + y_ax + '}$', fontsize=12)
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
        # For the PM data, add y coordinates to the same plot.
        if j == 8 and k == 1:
            if any(mmag_out_rjct) and any(stars_out_rjct_c):
                starsPlot('rjct', mmag_out_rjct,
                          list(zip(*list(zip(*stars_out_rjct_c))[j]))[k + 1])
            if any(mmag_in_rjct) and any(cl_region_rjct_c):
                starsPlot('rjct', mmag_in_rjct,
                          list(zip(*list(zip(*cl_region_rjct_c))[j]))[k + 1])
            if any(mmag_in_acpt) and any(cl_region_c):
                starsPlot('accpt_in', mmag_in_acpt,
                          list(zip(*list(zip(*cl_region_c))[j]))[k + 1])
            if any(mmag_out_acpt) and any(stars_out_c):
                starsPlot('accpt_out', mmag_out_acpt,
                          list(zip(*list(zip(*stars_out_c))[j]))[k + 1])

        # Plot legend in the main magnitude plot.
        if j == 4:
            # Legends.
            leg = plt.legend(fancybox=True, loc='upper left', scatterpoints=1,
                             fontsize=16, markerscale=2.5, prop={'size': 13})
            # Set the alpha value of the legend.
            leg.get_frame().set_alpha(0.7)

        if j in [4, 6]:
            ax.hlines(y=em_float[0], xmin=x_min, xmax=x_max, color='k',
                      linestyles='dashed', zorder=4)
        if j == 4:
            # Plot error curve
            plt.plot(err_bar_all[1], err_bar_all[2][0], color='#ffff00',
                     ls='--', zorder=5)
        elif j == 6:
            plt.plot(err_bar_all[1], err_bar_all[2][k + 1], color='#ffff00',
                     ls='--', zorder=5)
            if k == 0:
                ax.set_title(
                    r"$N_{{accpt}}={}$ , $N_{{rjct}}={}$ (compl frame)".format(
                        len(cl_region_c) + len(stars_out_c),
                        len(stars_out_rjct_c) + len(cl_region_rjct_c)),
                    fontsize=9)
        else:
            idx = k + 1 if k in (0, 1) else k
            ax.hlines(y=em_float[idx], xmin=x_min, xmax=x_max, color='k',
                      linestyles='dashed', zorder=4)
        # Maximum error limit of 1.
        plt.ylim(-0.005, min(plt.ylim()[1], 1.))


def pl_cl_fl_regions(
    gs, x_name, y_name, coord, x_min, x_max, y_min, y_max, asp_ratio,
        field_regions_rjct_c, cl_region_rjct_c, flag_no_fl_regs_c):
    '''
    Cluster and field regions defined.
    '''
    ax = plt.subplot(gs[2:4, 0:2])
    ax.set_aspect(aspect=asp_ratio)
    # Set plot limits
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    # If RA is used, invert axis.
    if coord == 'deg':
        ax.invert_xaxis()
    # Set axis labels
    plt.xlabel('{} ({})'.format(x_name, coord), fontsize=12)
    plt.ylabel('{} ({})'.format(y_name, coord), fontsize=12)
    # Set minor ticks
    ax.minorticks_on()
    ax.grid(b=True, which='both', color='gray', linestyle='--', lw=.5)

    # Plot cluster region.
    if len(cl_region_rjct_c) > 0:
        plt.scatter(
            list(zip(*cl_region_rjct_c))[1], list(zip(*cl_region_rjct_c))[2],
            marker='x', c='teal', s=15, lw=.5, edgecolors='none')

    N_flrg = 0
    if not flag_no_fl_regs_c:
        # Stars inside the field regions with rejected errors.
        for i, reg in enumerate(field_regions_rjct_c):
            fl_reg = list(zip(*reg))
            N_flrg += len(fl_reg[0])
            plt.scatter(fl_reg[1], fl_reg[2], marker='x',
                        c='teal', s=15, lw=.5, edgecolors='none')

    ax.set_title(r"$N_{{rjct}}$={} (phot compl)".format(
        len(cl_region_rjct_c) + N_flrg), fontsize=9)


def pl_fl_diag(
    gs, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd, x_ax, y_ax,
        field_regions_c, stars_f_rjct, stars_f_acpt, f_sz_pt, err_bar):
    '''
    Field stars CMD/CCD diagram.
    '''
    ax = plt.subplot(gs[2:4, 2:4])
    ax.set_title(r"$N_{{accpt}}={}$ , $N_{{rjct}}={}$ (fields compl)".format(
        len(stars_f_acpt[0]), len(stars_f_rjct[0])), fontsize=9)
    # Set plot limits
    plt.xlim(x_min_cmd, x_max_cmd)
    plt.ylim(y_min_cmd, y_max_cmd)
    # Set axis labels
    plt.xlabel('$' + x_ax + '$', fontsize=12)
    plt.ylabel('$' + y_ax + '$', fontsize=12)
    # Set minor ticks
    ax.minorticks_on()
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=.5,
            zorder=1)
    # Plot accepted/rejected stars within the field regions defined.
    if stars_f_rjct[0]:
        plt.scatter(stars_f_rjct[0], stars_f_rjct[1], marker='x',
                    c='teal', s=15, lw=.5, zorder=2)
    if stars_f_acpt[0]:
        plt.scatter(stars_f_acpt[0], stars_f_acpt[1], marker='o', c='b',
                    s=f_sz_pt, lw=0.3, edgecolor='k', zorder=3)
        n_field = int(len(stars_f_acpt[0]) / float(len(field_regions_c)))
        # Add text box.
        text = r'$n_{{field}} \approx {}$'.format(n_field)
        ob = offsetbox.AnchoredText(text, pad=0.2, loc=1, prop=dict(size=12))
        ob.patch.set(alpha=0.7)
        ax.add_artist(ob)
    # If list is not empty, plot error bars at several values.
    x_val, mag_y, xy_err = err_bar
    if x_val:
        plt.errorbar(
            x_val, mag_y, yerr=xy_err[0], xerr=xy_err[1], fmt='k.', lw=0.8,
            ms=0., zorder=4)


def pl_cl_diag(
    gs, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd, x_ax, y_ax,
        cl_region_rjct_c, cl_region_c, n_memb, cl_sz_pt, err_bar):
    '''
    Cluster's stars diagram (stars inside cluster's radius)
    '''
    ax = plt.subplot(gs[2:4, 4:6])
    ax.set_title(
        r"$N_{{accpt}}={}$ , $N_{{rjct}}={}$"
        r" ($r \leq r_{{cl}}$ compl)".format(
            len(cl_region_c), len(cl_region_rjct_c)), fontsize=9)
    # Set plot limits
    plt.xlim(x_min_cmd, x_max_cmd)
    plt.ylim(y_min_cmd, y_max_cmd)
    # Set axis labels
    plt.xlabel('$' + x_ax + '$', fontsize=12)
    plt.ylabel('$' + y_ax + '$', fontsize=12)
    # Set minor ticks
    ax.minorticks_on()
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=.5,
            zorder=1)
    # Add text box.
    text = r'$n_{{memb}} \approx {}$'.format(n_memb)
    ob = offsetbox.AnchoredText(text, pad=0.2, loc=1, prop=dict(size=12))
    ob.patch.set(alpha=0.7)
    ax.add_artist(ob)
    # Plot stars in CMD.
    if len(cl_region_rjct_c) > 0:
        # Only attempt to plot if any star is stored in the list.
        plt.scatter(
            list(zip(*list(zip(*cl_region_rjct_c))[5]))[0],
            list(zip(*list(zip(*cl_region_rjct_c))[3]))[0],
            marker='x', c='teal', s=12, lw=.5, zorder=2)
    plt.scatter(
        list(zip(*list(zip(*cl_region_c))[5]))[0],
        list(zip(*list(zip(*cl_region_c))[3]))[0],
        marker='o', c='r', s=cl_sz_pt, lw=0.3, edgecolor='k', zorder=3)
    # If list is not empty, plot error bars at several values.
    x_val, mag_y, xy_err = err_bar
    if x_val:
        plt.errorbar(
            x_val, mag_y, yerr=xy_err[0], xerr=xy_err[1], fmt='k.', lw=0.8,
            ms=0., zorder=4)


def pl_lum_func(gs, y_ax, flag_no_fl_regs, lum_func):
    '''
    LF of stars in cluster region and outside.
    '''
    x_cl, y_cl, x_fl, y_fl, x_all, y_all = lum_func
    ax = plt.subplot(gs[4:6, 0:2])
    ax.set_title("LF after error removal (compl)", fontsize=9)
    ax.minorticks_on()
    # Only draw units on axis (ie: 1, 2, 3)
    ax.xaxis.set_major_locator(MultipleLocator(2.0))
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=.5,
            zorder=1)
    # Set axis labels
    plt.xlabel('$' + y_ax + '$', fontsize=12)
    plt.ylabel('$N^{\star}/A_{cl}$', fontsize=12)

    # All frame.
    plt.step(x_all, y_all, where='post', color='k', lw=2.5, linestyle=':',
             label='Frame (compl)', zorder=6)
    # Cluster region LF (contaminated).
    plt.step(x_cl, y_cl, where='post', color='r', lw=1.,
             label='$LF_{cl+fl} \,(r \leq r_{cl})$', zorder=2)
    # Check if field regions were defined.
    if flag_no_fl_regs is not True:
        # Average field regions LF.
        plt.step(x_fl, y_fl, where='post', color='b', lw=1.,
                 label='$LF_{fl} \,(\star_{field})$', zorder=3)
        # Cluster region LF - average field regions LF.
        plt.step(x_cl, y_cl - y_fl, where='post', color='g', lw=1.7,
                 label='$LF_{cl}$', zorder=4)
        max_y = max(max(y_cl), max(y_fl), max(y_all))
    else:
        max_y = max(y_cl)
    # Set plot limits
    x_min, x_max = x_cl[-1] + .3, x_cl[1] - .3
    plt.xlim(x_min, x_max)
    plt.ylim(0., max_y + 0.05 * max_y)

    # Legends.
    leg = plt.legend(fancybox=True, loc='upper right', numpoints=1,
                     fontsize=11)
    # Set the alpha value of the legend.
    leg.get_frame().set_alpha(0.7)


def pl_data_rm_perc(
    gs, y_ax, phot_analy_compl, phot_data_compl, err_rm_data,
        combined_compl):
    """
    """
    ax = plt.subplot(gs[4:6, 2:4])
    ax.set_title("Percentage of stars kept after each process", fontsize=9)
    ax.minorticks_on()
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=.5,
            zorder=1)
    # Set axis labels
    plt.xlabel('$' + y_ax + '$', fontsize=12)
    plt.ylabel('perc', fontsize=12)

    edges, perc_vals = phot_analy_compl
    perc_vals_min = [min(perc_vals)]
    txt = "Photometric analysis\ncompleteness"
    plt.step(edges[:-1], perc_vals, where='post', lw=2., linestyle='--',
             label=txt)

    edges, perc_vals, perc_rmvd = phot_data_compl
    perc_vals_min.append(min(perc_vals))
    txt = "Photometric data\ncompleteness\n" +\
        "({:.1f}% rmvd)".format(perc_rmvd)
    plt.step(
        edges[:-1], perc_vals, where='post', lw=2., linestyle='--', label=txt)

    edges, perc_vals, perc_rmvd = err_rm_data
    perc_vals_min.append(min(perc_vals))
    txt = "Error removal\n({:.1f}% rmvd)".format(perc_rmvd)
    plt.step(
        edges[:-1], perc_vals, where='post', lw=2., color='teal',
        linestyle='--', label=txt)

    edges, perc_vals, perc_rmvd = combined_compl
    # Reverse.
    perc_vals = 1. - perc_vals
    perc_vals_min.append(min(perc_vals))
    txt = "Combined function\n({:.1f}% rmvd)".format(perc_rmvd)
    plt.step(
        edges[:-1], perc_vals, where='post', lw=2., color='r', linestyle='--',
        label=txt)

    # Legends.
    leg = plt.legend(
        fancybox=True, numpoints=1, loc='lower center', fontsize=9)
    # Set the alpha value of the legend.
    leg.get_frame().set_alpha(0.7)

    plt.gca().invert_xaxis()
    plt.ylim(min(.9, min(perc_vals_min)) - .05, 1.05)


# DEPRECATED 31/10/18
# def pl_p_vals(gs, flag_pval_test, pval_test_params):
#     '''
#     Distribution of KDE p_values.
#     '''
#     if flag_pval_test:
#         # Extract parameters from list.
#         prob_cl_kde, kde_cl_1d, kde_f_1d, x_kde, y_over = pval_test_params
#         ax = plt.subplot(gs[4:6, 2:4])
#         plt.xlim(-0.15, 1.15)
#         plt.ylim(0, 1.02)
#         plt.xlabel('p-values', fontsize=12)
#         plt.ylabel('Density (normalized)', fontsize=12)
#         ax.minorticks_on()
#         ax.grid(b=True, which='major', color='gray', linestyle='--', lw=1,
#                 zorder=1)
#         # Grid to background.
#         ax.set_axisbelow(True)
#         # Plot field vs field KDE.
#         if kde_f_1d.any():
#             max_kde = max(max(kde_f_1d), max(kde_cl_1d))
#             plt.plot(x_kde, kde_f_1d / max_kde, color='b', ls='-', lw=1.,
#                      label='$KDE_{fl}$', zorder=2)
#         else:
#             max_kde = max(kde_cl_1d)
#         # Plot cluster vs field KDE.
#         plt.plot(x_kde, kde_cl_1d / max_kde, color='r', ls='-', lw=1.,
#                  label='$KDE_{cl}$', zorder=2)
#         # Fill overlap.
#         if y_over:
#             plt.fill_between(x_kde, np.asarray(y_over) / max_kde, 0,
#                              color='grey', alpha='0.5')
#         text = '$P_{cl}^{KDE} = %0.2f$' % round(prob_cl_kde, 2)
#         plt.text(0.05, 0.92, text, transform=ax.transAxes,
#                  bbox=dict(facecolor='white', alpha=0.6), fontsize=12)
#         # Legend.
#         handles, labels = ax.get_legend_handles_labels()
#         leg = ax.legend(handles, labels, loc='upper right', numpoints=1,
#                         fontsize=12)
#         leg.get_frame().set_alpha(0.6)


def plot(N, *args):
    '''
    Handle each plot separately.
    '''

    plt_map = {
        0: [pl_phot_err, 'error rejection function'],
        1: [pl_cl_fl_regions, 'cluster + field regions rejected stars'],
        2: [pl_fl_diag, 'field regions photometric diagram'],
        3: [pl_cl_diag, 'cluster region photometric diagram'],
        4: [pl_lum_func, 'luminosity function'],
        5: [pl_data_rm_perc, 'error removal percentage']
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

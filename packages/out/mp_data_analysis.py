
import numpy as np
# from scipy.ndimage.filters import gaussian_filter
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
# from matplotlib.colors import LinearSegmentedColormap
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
    k_flag = np.array([_ == 'n' for _ in (pd_Plx, pd_PMRA, pd_RV)]).all()
    gs0 = 1 if k_flag else 0

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
        ax.grid(b=True, which='major', color='gray', linestyle='--', lw=.5,
                zorder=1)
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

        if j == 4:
            # Plot legend in the main magnitude plot.
            leg = plt.legend(fancybox=True, loc='upper left', scatterpoints=1,
                             fontsize=16, markerscale=2.5, prop={'size': 13})
            # Set the alpha value of the legend.
            leg.get_frame().set_alpha(0.7)
            # Max error cut
            ax.hlines(y=em_float[0], xmin=x_min, xmax=x_max, color='k',
                      linestyles='dashed', zorder=4)
            ob = offsetbox.AnchoredText(
                r"$max={}$ mag".format(em_float[0]), loc=1, prop=dict(size=9))
            ob.patch.set(alpha=0.7)
            ax.add_artist(ob)
            # Plot error curve
            plt.plot(err_bar_all[1], err_bar_all[2][0], color='#ffff00',
                     ls='--', zorder=5)

        elif j == 6:
            ax.hlines(y=em_float[1 + k], xmin=x_min, xmax=x_max, color='k',
                      linestyles='dashed', zorder=4)
            ob = offsetbox.AnchoredText(
                r"$max={}$ mag".format(em_float[1 + k]), loc=2,
                prop=dict(size=9))
            ob.patch.set(alpha=0.7)
            ax.add_artist(ob)
            plt.plot(err_bar_all[1], err_bar_all[2][k + 1], color='#ffff00',
                     ls='--', zorder=5)
            if k == 0:
                ax.set_title(
                    r"$N_{{accpt}}={}$ , $N_{{rjct}}={}$ (compl frame)".format(
                        len(cl_region_c) + len(stars_out_c),
                        len(stars_out_rjct_c) + len(cl_region_rjct_c)),
                    fontsize=9)
        else:
            unit = {0: '[mas]', 1: '[mas/yr]'}
            ax.hlines(y=em_float[-(3 - k)], xmin=x_min, xmax=x_max, color='k',
                      linestyles='dashed', zorder=4)
            ob = offsetbox.AnchoredText(
                r"$max={}$ {}".format(em_float[-(3 - k)], unit[k]), loc=2,
                prop=dict(size=9))
            ob.patch.set(alpha=0.7)
            ax.add_artist(ob)
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
            if reg:
                fl_reg = list(zip(*reg))
                N_flrg += len(fl_reg[0])
                plt.scatter(fl_reg[1], fl_reg[2], marker='x',
                            c='teal', s=15, lw=.5, edgecolors='none')

    ax.set_title(r"$N_{{rjct}}$={} (phot compl)".format(
        len(cl_region_rjct_c) + N_flrg), fontsize=9)


def pl_lum_func(gs, y_ax, flag_no_fl_regs, lum_func):
    '''
    LF of stars in cluster region and outside.
    '''
    x_cl, y_cl, x_fl, y_fl, x_all, y_all = lum_func
    ax = plt.subplot(gs[2:4, 2:4])
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
        max_y = max(max(y_cl), max(y_all))
    # Set plot limits
    x_min, x_max = x_cl[-1] + .3, x_cl[1] - .3
    plt.xlim(x_min, x_max)
    plt.ylim(0., max_y + 0.05 * max_y)

    # Legends.
    leg = plt.legend(fancybox=True, loc='upper right', numpoints=1,
                     fontsize='small')
    # Set the alpha value of the legend.
    leg.get_frame().set_alpha(0.7)


def pl_data_rm_perc(
    gs, y_ax, phot_analy_compl, phot_data_compl, err_rm_data,
        combined_compl):
    """
    """
    ax = plt.subplot(gs[2:4, 4:6])
    ax.set_title("Percentage of stars kept after each process", fontsize=9)
    ax.minorticks_on()
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=.5,
            zorder=1)
    # Set axis labels
    plt.xlabel('$' + y_ax + '$', fontsize=12)
    plt.ylabel('perc', fontsize=12)

    edges, perc_vals = phot_analy_compl
    perc_vals_min = [min(perc_vals)]
    txt = "Photometric analysis completeness"
    plt.step(edges[:-1], perc_vals, where='post', lw=2., linestyle='--',
             label=txt)

    edges, perc_vals, perc_rmvd = phot_data_compl
    perc_vals_min.append(min(perc_vals))
    txt = "Photometric data completeness ({:.1f}% rm)".format(perc_rmvd)
    plt.step(
        edges[:-1], perc_vals, where='post', lw=2., linestyle='--', label=txt)

    edges, perc_vals, perc_rmvd = err_rm_data
    perc_vals_min.append(min(perc_vals))
    txt = "Error removal ({:.1f}% rm)".format(perc_rmvd)
    plt.step(
        edges[:-1], perc_vals, where='post', lw=2., color='teal',
        linestyle='--', label=txt)

    edges, perc_vals, perc_rmvd = combined_compl
    # Reverse.
    perc_vals = 1. - perc_vals
    perc_vals_min.append(min(perc_vals))
    txt = "Combined function ({:.1f}% rm)".format(perc_rmvd)
    plt.step(
        edges[:-1], perc_vals, where='post', lw=2., color='r', linestyle='--',
        label=txt)

    # Legends.
    leg = plt.legend(
        fancybox=True, numpoints=1, loc='lower right', fontsize='small')
    # Set the alpha value of the legend.
    leg.get_frame().set_alpha(0.7)

    plt.gca().invert_xaxis()
    plt.ylim(min(.9, min(perc_vals_min)) - .05, 1.05)


def flCMD(
    ax, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd, x_ax, y_ax,
        N_fr, x_fr_rject, y_fr_rject, x_fr_accpt, y_fr_accpt, f_sz_pt,
        err_bar):
    '''
    Field stars CMD/CCD diagram.
    '''
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
    if x_fr_rject:
        plt.scatter(x_fr_rject, y_fr_rject, marker='x',
                    c='teal', s=15, lw=.5, zorder=2)
    if x_fr_accpt:
        plt.scatter(x_fr_accpt, y_fr_accpt, marker='o', c='b',
                    s=f_sz_pt, lw=0.3, edgecolor='k', zorder=3)
        n_field = int(len(x_fr_accpt) / float(N_fr))
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


def pl_fl_diag(
    gs, x_ax0, y_ax, x_min_cmd0, x_max_cmd0, y_min_cmd0, y_max_cmd0, x_ax1,
        x_min_cmd1, x_max_cmd1, y_min_cmd1, y_max_cmd1, field_regions_c,
        stars_f_rjct, stars_f_acpt, f_sz_pt, err_bar_fl0, err_bar_fl1):
    '''
    Field stars CMD/CCD diagram.
    '''
    ax = plt.subplot(gs[4:6, 4:6])

    N_fr = len(field_regions_c)
    x_fr_rject, y_fr_rject = stars_f_rjct[1], stars_f_rjct[0]
    x_fr_accpt, y_fr_accpt = stars_f_acpt[1], stars_f_acpt[0]

    ax.set_title(r"$N_{{accpt}}={}$ , $N_{{rjct}}={}$ (fields compl)".format(
        len(x_fr_accpt), len(x_fr_rject)), fontsize=9)

    flCMD(
        ax, x_min_cmd0, x_max_cmd0, y_min_cmd0, y_max_cmd0, x_ax0, y_ax,
        N_fr, x_fr_rject, y_fr_rject, x_fr_accpt, y_fr_accpt, f_sz_pt,
        err_bar_fl0)

    if x_ax1 != '':
        ax = plt.subplot(gs[6:8, 4:6])
        x_fr_rject, x_fr_accpt = stars_f_rjct[2], stars_f_acpt[2]
        flCMD(
            ax, x_min_cmd1, x_max_cmd1, y_min_cmd1, y_max_cmd1, x_ax1, y_ax,
            N_fr, x_fr_rject, y_fr_rject, x_fr_accpt, y_fr_accpt, f_sz_pt,
            err_bar_fl1)


def clCMD(
    ax, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd, x_ax, y_ax,
        xr, yr, xa, ya, n_memb, cl_sz_pt, err_bar):
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
    if yr:
        # Only attempt to plot if any star is stored in the list.
        plt.scatter(xr, yr, marker='x', c='teal', s=12, lw=.5, zorder=2)
    plt.scatter(
        xa, ya, marker='o', c='r', s=cl_sz_pt, lw=0.3, edgecolor='k', zorder=3)
    # If list is not empty, plot error bars at several values.
    x_val, mag_y, xy_err = err_bar
    if x_val:
        plt.errorbar(
            x_val, mag_y, yerr=xy_err[0], xerr=xy_err[1], fmt='k.', lw=0.8,
            ms=0., zorder=4)


def pl_cl_diag(
    gs, x_ax0, y_ax, x_min_cmd0, x_max_cmd0, y_min_cmd0, y_max_cmd0, x_ax1,
        x_min_cmd1, x_max_cmd1, y_min_cmd1, y_max_cmd1, err_bar_cl0,
        err_bar_cl1, cl_region_rjct_c, cl_region_c, n_memb, cl_sz_pt):
    '''
    Cluster's stars diagram (stars inside cluster's radius)
    '''
    ax = plt.subplot(gs[4:6, 0:2])
    ax.set_title(
        r"$N_{{accpt}}={}$ , $N_{{rjct}}={}$"
        r" ($r \leq r_{{cl}}$ compl)".format(
            len(cl_region_c), len(cl_region_rjct_c)), fontsize=9)
    xr, yr = [], []
    if len(cl_region_rjct_c) > 0:
        xr = list(zip(*list(zip(*cl_region_rjct_c))[5]))[0]
        yr = list(zip(*list(zip(*cl_region_rjct_c))[3]))[0]
    xa = list(zip(*list(zip(*cl_region_c))[5]))[0]
    ya = list(zip(*list(zip(*cl_region_c))[3]))[0]
    clCMD(
        ax, x_min_cmd0, x_max_cmd0, y_min_cmd0, y_max_cmd0, x_ax0, y_ax,
        xr, yr, xa, ya, n_memb, cl_sz_pt, err_bar_cl0)

    if x_ax1 != '':
        ax = plt.subplot(gs[6:8, 0:2])
        xr = []
        if len(cl_region_rjct_c) > 0:
            xr = list(zip(*list(zip(*cl_region_rjct_c))[5]))[1]
        xa = list(zip(*list(zip(*cl_region_c))[5]))[1]
        clCMD(
            ax, x_min_cmd1, x_max_cmd1, y_min_cmd1, y_max_cmd1, x_ax1, y_ax,
            xr, yr, xa, ya, n_memb, cl_sz_pt, err_bar_cl1)


def hessKDE(
    ax, x_ax, y_ax, x_max_cmd, x_min_cmd, y_min_cmd, y_max_cmd,
        cl_col, cl_mag, fr_col, fr_mag):

    # This bandwidth seems to produce nice results.
    bw, Nb = .2, 100

    ax.set_title("Cluster - Field (normalized)", fontsize=9)
    plt.xlabel('$' + x_ax + '$', fontsize=12)
    plt.ylabel('$' + y_ax + '$', fontsize=12)

    xx, yy = np.mgrid[x_min_cmd:x_max_cmd:complex(Nb),
                      y_max_cmd:y_min_cmd:complex(Nb)]
    positions = np.vstack([xx.ravel(), yy.ravel()])

    # # Test the impact of the bw value in the final estimate.
    # scott_f = len(cl_col)**(-1. / (2 + 4))
    # for bws in (scott_f / 2., scott_f, scott_f * 2.):
    #     kernel1 = gaussian_kde(np.vstack([cl_col, cl_mag]), bw_method=bws)
    #     kernel2 = gaussian_kde(np.vstack([fr_col, fr_mag]), bw_method=bws)
    #     f1 = np.reshape(kernel1(positions).T, xx.shape)
    #     f2 = np.reshape(kernel2(positions).T, xx.shape)
    #     diff = f1 - f2
    #     diff = np.clip(diff, 1e-9, np.inf)
    #     cell = ((x_max_cmd - x_min_cmd) * (y_min_cmd - y_max_cmd)) / Nb**2
    #     print(np.sum(diff * cell))

    # Cluster data
    values1 = np.vstack([cl_col, cl_mag])
    kernel1 = gaussian_kde(values1, bw_method=bw)
    f1 = np.reshape(kernel1(positions).T, xx.shape)

    # Field regions data
    values2 = np.vstack([fr_col, fr_mag])
    kernel2 = gaussian_kde(values2, bw_method=bw)
    f2 = np.reshape(kernel2(positions).T, xx.shape)

    # Cluster - field regions
    diff = f1 - f2
    # Clip negative values.
    diff = np.clip(diff, 1e-9, np.inf)

    # Area of the 2D cell.
    cell = ((x_max_cmd - x_min_cmd) * (y_min_cmd - y_max_cmd)) / Nb**2
    # Integral of the cluster-field KDE. This value strongly depends on
    # the selected 'bw' value, so it is not really a stable indicator of
    # field contamination in the cluster region.
    integ = np.sum(diff * cell)
    # Add text box.
    text = r'$\int \Delta KDE_{{[cl-fr]}} \approx {:.2f}$'.format(integ)
    ob = offsetbox.AnchoredText(text, pad=0.2, loc=1, prop=dict(size=9))
    ob.patch.set(alpha=0.7)
    ax.add_artist(ob)

    ax.contourf(xx, yy, diff, cmap='Blues')
    ax.contour(xx, yy, diff, colors='k', linewidths=.5)
    # ax.clabel(CS, inline=1, fontsize=10)

    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=.5,
            zorder=3)
    plt.gca().invert_yaxis()


def pl_hess_cmd(
    gs, x_ax0, x_ax1, y_ax, x_max_cmd0, x_min_cmd0, y_min_cmd0, y_max_cmd0,
        x_max_cmd1, x_min_cmd1, y_min_cmd1, y_max_cmd1, stars_f_acpt,
        cl_region_c):
    '''
    Hess diagram for CMD of field vs cluster region.
    '''
    if stars_f_acpt[0]:
        ax = plt.subplot(gs[4:6, 2:4])
        cl_col = list(zip(*list(zip(*cl_region_c))[5]))[0]
        cl_mag = list(zip(*list(zip(*cl_region_c))[3]))[0]
        fr_col, fr_mag = stars_f_acpt[1], stars_f_acpt[0]

        hessKDE(
            ax, x_ax0, y_ax, x_max_cmd0, x_min_cmd0, y_min_cmd0,
            y_max_cmd0, cl_col, cl_mag, fr_col, fr_mag)

        if stars_f_acpt[2]:
            cl_col = list(zip(*list(zip(*cl_region_c))[5]))[1]
            fr_col = stars_f_acpt[2]
            ax = plt.subplot(gs[6:8, 2:4])
            hessKDE(
                ax, x_ax1, y_ax, x_max_cmd1, x_min_cmd1, y_min_cmd1,
                y_max_cmd1, cl_col, cl_mag, fr_col, fr_mag)


def pl_ad_test(gs, b, flag_ad_test, ad_cl, ad_fr, id_kinem, ad_k_comb):
    """
    """
    if flag_ad_test:

        def adPlot(ax, d1, d2, s):
            ax.set_title('(' + s + ')', fontsize=9)
            ax.axes.yaxis.set_ticklabels([])
            plt.xlabel("A-D test")
            plt.ylabel("N (norm)")
            ax.grid(
                b=True, which='major', color='gray', linestyle='--', lw=.5,
                zorder=1)
            ax.hist(d1, bins=25, density=True, color='r', histtype='step')
            plt.plot([0, 0], label=r'$Cl$', color='r')
            ax.hist(d2, bins=25, density=True, color='b', histtype='step')
            plt.plot([0, 0], label=r'$Fr$', color='b')
            if d2:
                xmin, xmax = min(min(d1), min(d2)), max(max(d1), max(d2))
            else:
                xmin, xmax = min(d1), max(d1)
            if xmin < 0.325:
                ax.axvline(
                    x=0.325, ls=':', lw=2.5, c='orange', label=r"$p_{v}=0.25$")
            ax.axvline(x=3.752, ls=':', lw=2.5, c='g', label=r"$p_{v}=0.01$")
            if xmax > 6.546:
                ax.axvline(
                    x=6.546, ls=':', lw=2.5, c='k', label=r"$p_{v}=0.001$")
            ax.set_xscale('log')
            ax.legend(fontsize='small')

        ax = plt.subplot(gs[6 + b:7 + b, 0:2])
        adPlot(ax, ad_cl[0], ad_fr[0], 'phot')

        # Only plot if either parallax or PMs or radial velocities are defined
        pd_Plx, pd_PMRA, pd_RV = id_kinem[0], id_kinem[2], id_kinem[6]
        # Define first row depending on whether kinematic data was defined.
        k_flag = np.array([_ != 'n' for _ in (pd_Plx, pd_PMRA, pd_RV)]).any()
        if k_flag:
            s = 'all' if ad_k_comb else 'plx+pm+rv'
            ax = plt.subplot(gs[7 + b:8 + b, 0:2])
            adPlot(ax, ad_cl[1], ad_fr[1], s)


def pl_p_vals(
    ax, Ncl, Nf, prob_cl, kde_cl, kde_fr, x_cl, x_fr, x_over, y_over,
        reg_id):
    ax.set_title(
        r'$P_{{cl}}={:.2f}\;({})$'.format(prob_cl, reg_id), fontsize=9)
    ax.axes.yaxis.set_ticklabels([])
    plt.xlabel('p-values', fontsize=12)
    plt.ylabel('Density (norm)', fontsize=12)
    ax.minorticks_on()
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=.5,
            zorder=1)
    # Grid to background.
    ax.set_axisbelow(True)
    # Plot field vs field KDE.
    if kde_fr.any():
        plt.plot(x_fr, kde_fr, color='b', ls='-', lw=1.,
                 label=r'$Fr\,({})$'.format(Nf), zorder=2)
    # Plot cluster vs field KDE.
    if not kde_cl.any():
        ax.axvline(
            x=x_cl, c='r', ls='-', lw=1., label=r'$Cl\,({})$'.format(Ncl))
    else:
        plt.plot(x_cl, kde_cl, color='r', ls='-', lw=1.,
                 label=r'$Cl\,({})$'.format(Ncl), zorder=2)
    # Fill overlap.
    if y_over.any():
        plt.fill_between(x_over, y_over, 0, color='grey', alpha='0.5')
    # Legend.
    handles, labels = ax.get_legend_handles_labels()
    leg = ax.legend(handles, labels, numpoints=1, fontsize=9)
    leg.get_frame().set_alpha(0.6)
    plt.gca().set_ylim(bottom=0)


def pl_ad_pvals_phot(gs, b, flag_ad_test, ad_cl_fr_p):
    if flag_ad_test:
        Ncl, Nf, prob_cl, kde_cl, kde_fr, x_cl, x_fr, x_over, y_over =\
            ad_cl_fr_p
        ax = plt.subplot(gs[6 + b:8 + b, 2:4])
        pl_p_vals(
            ax, Ncl, Nf, prob_cl, kde_cl, kde_fr, x_cl, x_fr, x_over, y_over,
            'phot')


def pl_ad_pvals_pk(gs, b, flag_ad_test, ad_cl_fr_pk, id_kinem, ad_k_comb):
    # Only plot if either parallax or PMs or radial velocities are defined
    pd_Plx, pd_PMRA, pd_RV = id_kinem[0], id_kinem[2], id_kinem[6]
    # Define first row depending on whether kinematic data was defined.
    k_flag = np.array([_ != 'n' for _ in (pd_Plx, pd_PMRA, pd_RV)]).any()
    if flag_ad_test and k_flag:
        Ncl, Nf, prob_cl, kde_cl, kde_fr, x_cl, x_fr, x_over, y_over =\
            ad_cl_fr_pk
        ax = plt.subplot(gs[6 + b:8 + b, 4:6])
        s = 'all' if ad_k_comb else 'plx+pm+rv'
        pl_p_vals(
            ax, Ncl, Nf, prob_cl, kde_cl, kde_fr, x_cl, x_fr, x_over, y_over,
            s)


def plot(N, *args):
    '''
    Handle each plot separately.
    '''

    plt_map = {
        0: [pl_phot_err, 'error rejection function'],
        1: [pl_cl_fl_regions, 'cluster + field regions rejected stars'],
        2: [pl_fl_diag, 'field regions photometric diagram'],
        3: [pl_cl_diag, 'cluster region photometric diagram'],
        4: [pl_hess_cmd, 'Hess CMD'],
        5: [pl_lum_func, 'luminosity function'],
        6: [pl_data_rm_perc, 'error removal percentage'],
        7: [pl_ad_test, 'A-D test values'],
        8: [pl_ad_pvals_phot, 'photometric A-D pvalues'],
        9: [pl_ad_pvals_pk, 'photometric + kinem A-D pvalues']
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

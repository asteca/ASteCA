
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib.offsetbox as offsetbox


def pl_phot_err(gs, fig, up_low, x_ax, y_ax, mags, err_max, cl_region,
                stars_in_rjct, stars_out, stars_out_rjct, err_bar_all):
    '''
    Photometric error rejection.
    '''
    # Define parameters for upper and lower plots.
    if up_low == 'up':
        ax, ax_y, j = plt.subplot(gs[0, 0:2]), y_ax, 4
    else:
        ax, ax_y, j = plt.subplot(gs[1, 0:2]), x_ax, 6

    # Set plot limits
    x_min, x_max = min(mags[0]) - 0.5, max(mags[0]) + 0.5
    plt.xlim(x_min, x_max)
    # Set axis labels
    plt.ylabel('$\sigma_{' + ax_y + '}$', fontsize=18)
    plt.xlabel('$' + y_ax + '$', fontsize=18)
    ax.set_facecolor('#EFF0F1')
    # Set minor ticks
    ax.minorticks_on()
    # Plot rejected stars.
    if len(stars_out_rjct) > 0:
        # Only attempt to pot if any star is stored in the list.
        plt.scatter(
            zip(*zip(*stars_out_rjct)[3])[0], zip(*zip(*stars_out_rjct)[j])[0],
            marker='x', c='teal', s=15, zorder=1)
    if len(stars_in_rjct) > 0:
        plt.scatter(
            zip(*zip(*stars_in_rjct)[3])[0], zip(*zip(*stars_in_rjct)[j])[0],
            marker='x', c='teal', s=15, zorder=1)
    # Plot accepted stars.
    plt.scatter(
        zip(*zip(*stars_out)[3])[0], zip(*zip(*stars_out)[j])[0], marker='o',
        c='b', s=5, zorder=2, lw=0.1, edgecolor='k', label='$r > r_{cl}$')
    plt.scatter(
        zip(*zip(*cl_region)[3])[0], zip(*zip(*cl_region)[j])[0], marker='o',
        c='r', s=10, zorder=2, lw=0.3, edgecolor='k', label='$r \leq r_{cl}}$')
    if up_low == 'up':
        # Legends.
        leg = plt.legend(fancybox=True, loc='upper left', scatterpoints=1,
                         fontsize=16, markerscale=2.5, prop={'size': 13})
        # Set the alpha value of the legend.
        leg.get_frame().set_alpha(0.7)
    # Plot error curve
    if up_low == 'up':
        plt.plot(err_bar_all[1], err_bar_all[3], color='#ffff00', ls='--',
                 zorder=4)
    else:
        plt.plot(err_bar_all[1], err_bar_all[2], color='#ffff00', ls='--',
                 zorder=4)
    ax.hlines(y=err_max, xmin=x_min, xmax=x_max, color='k',
              linestyles='dashed', zorder=2)
    plt.ylim(-0.005, plt.ylim()[1])


def pl_fl_diag(gs, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd, x_ax, y_ax,
               stars_f_rjct, stars_f_acpt, f_sz_pt, err_bar):
    '''
    Field stars CMD/CCD diagram.
    '''
    ax = plt.subplot(gs[0:2, 2:4])
    # Set plot limits
    plt.xlim(x_min_cmd, x_max_cmd)
    plt.ylim(y_min_cmd, y_max_cmd)
    # Set axis labels
    plt.xlabel('$' + x_ax + '$', fontsize=18)
    plt.ylabel('$' + y_ax + '$', fontsize=18)
    # Set minor ticks
    ax.minorticks_on()
    # Only draw units on axis (ie: 1, 2, 3)
    ax.xaxis.set_major_locator(MultipleLocator(1.0))
    # Set grid
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=1,
            zorder=1)
    # Plot *all* rejected stars outside of the cluster region.
    plt.scatter(stars_f_rjct[0], stars_f_rjct[1], marker='x',
                c='teal', s=15, zorder=2)
    # Add text box.
    text = '$N_{{accpt}}={}\,(\star_{{field}})$'.format(len(stars_f_acpt[0]))
    ob = offsetbox.AnchoredText(text, pad=0.2, loc=1, prop=dict(size=12))
    ob.patch.set(alpha=0.7)
    ax.add_artist(ob)
    # Plot accepted stars within the field regions defined, if at least one
    # exists.
    if stars_f_acpt[0]:
        plt.scatter(stars_f_acpt[0], stars_f_acpt[1], marker='o', c='b',
                    s=f_sz_pt, lw=0.3, edgecolor='k', zorder=3)
    # If list is not empty, plot error bars at several values.
    x_val, mag_y, x_err, y_err = err_bar
    if x_val:
        plt.errorbar(x_val, mag_y, yerr=y_err, xerr=x_err, fmt='k.', lw=0.8,
                     ms=0., zorder=4)


def pl_cl_diag(gs, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd, x_ax, y_ax,
               cl_region_rjct, cl_region, n_memb, cl_sz_pt, err_bar):
    '''
    Cluster's stars diagram (stars inside cluster's radius)
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
    # Only draw units on axis (ie: 1, 2, 3)
    ax.xaxis.set_major_locator(MultipleLocator(1.0))
    # Set grid
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=1,
            zorder=1)
    # Add text box.
    text1 = '$N_{{accpt}}={}\,(r \leq r_{{cl}})$'.format(len(cl_region))
    text2 = r'$n_{{memb}} \approx {}$'.format(n_memb)
    text = text1 + '\n' + text2
    ob = offsetbox.AnchoredText(text, pad=0.2, loc=1, prop=dict(size=12))
    ob.patch.set(alpha=0.7)
    ax.add_artist(ob)
    # Plot stars in CMD.
    if len(cl_region_rjct) > 0:
        # Only attempt to plot if any star is stored in the list.
        plt.scatter(
            zip(*zip(*cl_region_rjct)[5])[0], zip(*zip(*cl_region_rjct)[3])[0],
            marker='x', c='teal', s=12, zorder=2)
    plt.scatter(
        zip(*zip(*cl_region)[5])[0], zip(*zip(*cl_region)[3])[0], marker='o',
        c='r', s=cl_sz_pt, lw=0.3, edgecolor='k', zorder=3)
    # If list is not empty, plot error bars at several values.
    x_val, mag_y, x_err, y_err = err_bar
    if x_val:
        plt.errorbar(x_val, mag_y, yerr=y_err, xerr=x_err, fmt='k.', lw=0.8,
                     ms=0., zorder=4)


def pl_lum_func(gs, y_ax, flag_no_fl_regs, lum_func, completeness):
    '''
    LF of stars in cluster region and outside.
    '''
    x_cl, y_cl, x_fl, y_fl, x_all, y_all = lum_func
    ax = plt.subplot(gs[2:4, 0:2])
    ax.minorticks_on()
    # Only draw units on axis (ie: 1, 2, 3)
    ax.xaxis.set_major_locator(MultipleLocator(2.0))
    # Set grid
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=1,
            zorder=1)
    # Set axis labels
    plt.xlabel('$' + y_ax + '$', fontsize=18)
    plt.ylabel('$N^{\star}/A_{cl}$', fontsize=18)

    # All frame.
    plt.step(x_all, y_all, where='post', color='k', lw=1.5, linestyle=':',
             label='Frame', zorder=2)
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
        max_y = max(max(y_cl), max(y_fl))
    else:
        max_y = max(y_cl)
    # Set plot limits
    x_min, x_max = x_cl[-1] + .3, x_cl[1] - .3
    plt.xlim(x_min, x_max)
    plt.ylim(0., max_y + 0.05 * max_y)
    # Completeness maximum value.
    # completeness = [bin_edges, max_indx, comp_perc]
    bin_edges, max_indx = completeness[0], completeness[1]
    mag_peak = bin_edges[max_indx]
    text = '$' + y_ax + r',_{compl}\,\approx\,%0.1f$' % mag_peak
    ax.vlines(x=mag_peak, ymin=0., ymax=plt.ylim()[1], color='k',
              lw=1.5, linestyles='dashed', label=text, zorder=1)
    # Legends.
    leg = plt.legend(fancybox=True, loc='upper right', numpoints=1,
                     fontsize=12)
    # Set the alpha value of the legend.
    leg.get_frame().set_alpha(0.7)


def pl_integ_mag(gs, cl_reg_imag, fl_reg_imag, integ_mag, y_ax,
                 flag_no_fl_regs):
    '''
    Integrated magnitudes.
    '''
    if cl_reg_imag:
        # Make plot
        ax = plt.subplot(gs[2:4, 2:4])
        # If field lists are not empty.
        if fl_reg_imag:
            x_min = min(min(cl_reg_imag[0][0]), min(fl_reg_imag[0][0])) - 0.2
            x_max = max(max(cl_reg_imag[0][0]), max(fl_reg_imag[0][0])) + 0.2
            y_min = max(max(cl_reg_imag[0][1]), max(fl_reg_imag[0][1])) + 0.2
            y_max = min(min(cl_reg_imag[0][1]), min(fl_reg_imag[0][1])) - 0.2
        else:
            x_min, x_max = min(cl_reg_imag[0][0]) - 0.2,\
                max(cl_reg_imag[0][0]) + 0.2
            y_min, y_max = max(cl_reg_imag[0][1]) + 0.2,\
                min(cl_reg_imag[0][1]) - 0.2
        plt.xlim(x_min, x_max)
        plt.ylim(y_min, y_max)
        ax.set_xlabel('$mag$', fontsize=18)
        ax.set_ylabel('$mag^*$', fontsize=18)
        ax.minorticks_on()
        ax.grid(b=True, which='major', color='gray', linestyle='--', lw=1,
                zorder=1)
        text1 = '$' + y_ax + '^{*}_{cl+fl}$'
        # Cluster + field integrated magnitude curve.
        plt.plot(cl_reg_imag[0][0], cl_reg_imag[0][1], 'r-', lw=1.,
                 label=text1, zorder=2)
        # Check if field regions were defined.
        if not flag_no_fl_regs:
            text3 = '$' + y_ax + '^{*}_{fl}$'
            # Field average integrated magnitude curve.
            plt.plot(fl_reg_imag[0][0], fl_reg_imag[0][1], 'b-', lw=1.,
                     label=text3, zorder=2)
        text = r'$' + y_ax + '^{{*}}_{{cl}} = {:.2f}$'.format(integ_mag[0])
        ob = offsetbox.AnchoredText(text, pad=0.2, loc=7, prop=dict(size=12))
        ob.patch.set(alpha=0.75)
        ax.add_artist(ob)

        lines, labels = ax.get_legend_handles_labels()
        leg = ax.legend(lines, labels, loc='lower right', numpoints=1,
                        fontsize=12)
        leg.get_frame().set_alpha(0.7)


def pl_p_vals(gs, flag_pval_test, pval_test_params):
    '''
    Distribution of KDE p_values.
    '''
    if flag_pval_test:
        # Extract parameters from list.
        prob_cl_kde, kde_cl_1d, kde_f_1d, x_kde, y_over = pval_test_params
        ax = plt.subplot(gs[2:4, 4:6])
        plt.xlim(-0.15, 1.15)
        plt.ylim(0, 1.02)
        plt.xlabel('p-values', fontsize=12)
        plt.ylabel('Density (normalized)', fontsize=12)
        ax.minorticks_on()
        ax.grid(b=True, which='major', color='gray', linestyle='--', lw=1,
                zorder=1)
        # Grid to background.
        ax.set_axisbelow(True)
        # Plot field vs field KDE.
        if kde_f_1d.any():
            max_kde = max(max(kde_f_1d), max(kde_cl_1d))
            plt.plot(x_kde, kde_f_1d / max_kde, color='b', ls='-', lw=1.,
                     label='$KDE_{fl}$', zorder=2)
        else:
            max_kde = max(kde_cl_1d)
        # Plot cluster vs field KDE.
        plt.plot(x_kde, kde_cl_1d / max_kde, color='r', ls='-', lw=1.,
                 label='$KDE_{cl}$', zorder=2)
        # Fill overlap.
        if y_over:
            plt.fill_between(x_kde, np.asarray(y_over) / max_kde, 0,
                             color='grey', alpha='0.5')
        text = '$P_{cl}^{KDE} = %0.2f$' % round(prob_cl_kde, 2)
        plt.text(0.05, 0.92, text, transform=ax.transAxes,
                 bbox=dict(facecolor='white', alpha=0.6), fontsize=12)
        # Legend.
        handles, labels = ax.get_legend_handles_labels()
        leg = ax.legend(handles, labels, loc='upper right', numpoints=1,
                        fontsize=12)
        leg.get_frame().set_alpha(0.6)


def plot(N, *args):
    '''
    Handle each plot separately.
    '''

    plt_map = {
        0: [pl_phot_err, 'upper error rejection function'],
        1: [pl_phot_err, 'lower error rejection function'],
        2: [pl_fl_diag, 'field regions photometric diagram'],
        3: [pl_cl_diag, 'cluster region photometric diagram'],
        4: [pl_lum_func, 'luminosity function'],
        5: [pl_integ_mag, 'integrated magnitudes'],
        6: [pl_p_vals, 'KDE p-values distribution']
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


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.offsetbox as offsetbox
from matplotlib.patches import Ellipse
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy.interpolate


def pl_GA_lkl(gs, l_min_max, lkl_old, model_done, new_bs_indx, N_pop, N_gen,
              fit_diff, cross_prob, cross_sel, mut_prob, N_el, N_ei, N_es,
              N_b):
    '''
    Likelihood evolution for the GA.
    '''
    ax = plt.subplot(gs[0:2, 0:4])
    plt.xlim(-0.5, len(lkl_old[0]) + int(0.01 * len(lkl_old[0])))
    plt.ylim(l_min_max[0], l_min_max[1])
    # Set minor ticks
    ax.minorticks_on()
    ax.tick_params(axis='y', which='major', labelsize=9)
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=0.6)
    plt.xlabel('Generation', fontsize=12)
    plt.ylabel('Likelihood', fontsize=12)

    # Plot likelihood minimum and mean lines.
    ax.plot(range(len(lkl_old[0])), lkl_old[0], lw=1., c='green',
            label='$L_{{min}}={:.1f}$'.format(min(lkl_old[0])))
    ax.plot(range(len(lkl_old[0])), lkl_old[1], lw=1., c='k',
            label='$L_{mean}$')
    # Plot line marking a new best solution found.
    for lin in new_bs_indx:
        lw_lin = 1.5 if lin > 0.05 * len(lkl_old[0]) else 0.5
        plt.axvline(x=lin, linestyle='--', lw=lw_lin, color='blue')
    # Add text box.
    text1 = '$N_{{total}}={:.1E}\,;\,N_{{btst}}={}$'.format(
        len(model_done[0]), N_b)
    text2 = '$n_{{gen}}={}\,;\,n_{{pop}}={}$'.format(N_gen, N_pop)
    text3 = '$f_{{dif}}={:.2f}\,;\,cr_{{sel}}={}$'.format(
        fit_diff, cross_sel)
    text4 = '$p_{{cross}}={:.2f}\,;\,p_{{mut}}={:.2f}$'.format(
        cross_prob, mut_prob)
    text5 = '$n_{{el}}={}\,;\,n_{{ei}}={}\,;\,n_{{es}}={}$'.format(
        N_el, N_ei, N_es)
    text = text1 + '\n' + text2 + '\n' + text3 + '\n' + text4 + '\n' +\
        text5
    ob = offsetbox.AnchoredText(text, loc=1, prop=dict(size=12))
    ob.patch.set(boxstyle='square,pad=0.', alpha=0.85)
    ax.add_artist(ob)
    # Legend.
    handles, labels = ax.get_legend_handles_labels()
    leg = ax.legend(handles, labels, loc='upper left', numpoints=1,
                    fontsize=12)
    leg.get_frame().set_alpha(0.7)


def selectMinLkl(x, y, z):
    """
    Select the xy pairs with the smallest z values.
    Source: https://stackoverflow.com/a/45887552/1391441
    """
    # Get grouped lex-sort indices
    sidx = np.lexsort([x, y])
    # Lex-sort x, y, z
    x_sorted = x[sidx]
    y_sorted = y[sidx]
    z_sorted = z[sidx]

    # Get equality mask between each sorted X and Y elem against previous
    # ones. The non-zero indices of its inverted mask gives us the indices
    # where the new groupings start. We are calling those as cut_idx.
    seq_eq_mask = (x_sorted[1:] == x_sorted[:-1]) &\
        (y_sorted[1:] == y_sorted[:-1])
    cut_idx = np.flatnonzero(np.concatenate(([True], ~seq_eq_mask)))

    # Use those cut_idx to get intervalled minimum values
    minZ = np.minimum.reduceat(z_sorted, cut_idx)

    # Make tuples of the groupings of x,y and the corresponding min Z
    # values.
    return x_sorted[cut_idx], y_sorted[cut_idx], minZ.tolist()


def pl_2_param_dens(gs, _2_params, min_max_p, cp_r, cp_e, model_done):
    '''
    Parameter vs parameters solutions density map.
    '''
    # Define parameters for upper and lower plots.
    if _2_params == 'age-metal':
        gs_x1, gs_x2, cp, d_map, mx, my = 4, 6, 'r', 'GnBu_r', 0, 1
        x_label, y_label = '$z$', '$log(age)$'
    elif _2_params == 'dist-ext':
        gs_x1, gs_x2, cp, d_map, mx, my = 6, 8, 'b', 'YlOrBr_r', 2, 3
        x_label, y_label = '$E_{(B-V)}$', '$(m-M)_o$'
    elif _2_params == 'ext-age':
        gs_x1, gs_x2, cp, d_map, mx, my = 8, 10, 'r', 'GnBu_r', 2, 1
        x_label, y_label = '$E_{(B-V)}$', '$log(age)$'
    elif _2_params == 'mass-binar':
        gs_x1, gs_x2, cp, d_map, mx, my = 10, 12, 'b', 'YlOrBr_r', 4, 5
        x_label, y_label = '$M\,(M_{{\odot}})$', '$b_{frac}$'

    ax = plt.subplot(gs[0:2, gs_x1:gs_x2])
    # Parameter values and errors.
    xp, e_xp = map(float, [cp_r[mx], cp_e[mx]])
    yp, e_yp = map(float, [cp_r[my], cp_e[my]])
    # Axis limits.
    xp_min, xp_max = min_max_p[mx]
    yp_min, yp_max = min_max_p[my]
    plt.xlim(xp_min, xp_max)
    plt.ylim(yp_min, yp_max)
    # To specify the number of ticks on both or any single axes
    ax.locator_params(nbins=5)
    plt.xlabel(x_label, fontsize=16)
    plt.ylabel(y_label, fontsize=16)
    plt.minorticks_on()
    # Check if errors in both dimensions are defined.
    if all([~np.isnan(_) for _ in [e_xp, e_yp]]):
        # Plot ellipse error.
        plt.gca()
        ellipse = Ellipse(xy=(xp, yp), width=2 * e_xp, height=2 * e_yp,
                          edgecolor=cp, fc='None', lw=1., zorder=4)
        ax.add_patch(ellipse)
    # Else plot an error bar in the corresponding dimension.
    elif np.isnan(e_xp) and ~np.isnan(e_yp):
        plt.errorbar(xp, yp, yerr=e_yp, color=cp, zorder=4)
    elif np.isnan(e_yp) and ~np.isnan(e_xp):
        plt.errorbar(xp, yp, xerr=e_xp, color=cp, zorder=4)
    # Plot best fit point.
    plt.scatter(xp, yp, marker='x', c=cp, s=30, linewidth=2, zorder=4)

    # Minimum likelihood for each (x,y) pair in the density plot.
    z_lkl = np.log(np.asarray(model_done[1]) + 1.)
    x, y, z = selectMinLkl(
        np.array(zip(*model_done[0])[mx]), np.array(zip(*model_done[0])[my]),
        z_lkl)

    # Plot density map.
    xmin, xmax, ymin, ymax = min(x), max(x), min(y), max(y)
    # If at least one of the parameters was not fixed.
    if xmin != xmax or ymin != ymax:
        if xmin == xmax:
            xmin, xmax = xp_min, xp_max
        if ymin == ymax:
            ymin, ymax = yp_min, yp_max

        # Set up a regular grid of interpolation points
        xi, yi = np.linspace(xmin, xmax, 200), np.linspace(ymin, ymax, 200)
        xi, yi = np.meshgrid(xi, yi)
        # Normalize data and grid.
        # Source: https://stackoverflow.com/a/3867302/1391441
        x_new, x_grid = (x - xmin) / (xmax - xmin), (xi - xmin) / (xmax - xmin)
        y_new, y_grid = (y - ymin) / (ymax - ymin), (yi - ymin) / (ymax - ymin)

        if xmin != xmax and ymin != ymax and len(x) > 2500:
            # Use 'griddata' if no parameter was fixed, and the number of
            # unique solutions is large.
            zi = scipy.interpolate.griddata(
                (x_new, y_new), z, (x_grid, y_grid), method='linear')
            plt.imshow(zi, vmin=min(z), vmax=max(z), origin='lower',
                       extent=[xmin, xmax, ymin, ymax],
                       cmap=plt.get_cmap(d_map), zorder=2)
            ax.set_aspect('auto')
        else:
            # Use 'Rbf' if one parameter was fixed, or if the number of
            # solutions is small.
            # Source: https://stackoverflow.com/a/9008576/1391441
            rbf = scipy.interpolate.Rbf(x_new, y_new, z, function='linear')
            zi = rbf(x_grid, y_grid)
            plt.pcolormesh(xi, yi, zi, cmap=plt.get_cmap(d_map), zorder=2)

        plt.contour(xi, yi, zi, 5, colors='#551a8b', linewidths=0.5, zorder=3)


def pl_lkl_scatt(gs, ld_p, min_max_p, cp_r, cp_e, model_done):
    '''
    Parameter likelihood scatter plot.
    '''
    # Define parameters for upper and lower plots.
    if ld_p == '$z$':
        ax, cp, ci, zlab = plt.subplot(gs[2:4, 0:2]), 0, 1, '$log(age)$'
        plt.ylabel('Likelihood', fontsize=12)
    elif ld_p == '$log(age)$':
        ax, cp, ci, zlab = plt.subplot(gs[2:4, 2:4]), 1, 2, '$E_{{(B-V)}}$'
    elif ld_p == '$E_{{(B-V)}}$':
        ax, cp, ci, zlab = plt.subplot(gs[2:4, 4:6]), 2, 3, '$(m-M)_o$'
    elif ld_p == '$(m-M)_o$':
        ax, cp, ci, zlab = plt.subplot(gs[2:4, 6:8]), 3, 1, '$log(age)$'
    elif ld_p == '$M\,(M_{{\odot}})$':
        ax, cp, ci, zlab = plt.subplot(gs[2:4, 8:10]), 4, 5, '$b_{{frac}}$'
    elif ld_p == '$b_{{frac}}$':
        ax, cp, ci, zlab = plt.subplot(gs[2:4, 10:12]), 5, 4,\
            '$M\,(M_{{\odot}})$'

    # Parameter values and errors.
    xp, e_xp = map(float, [cp_r[cp], cp_e[cp]])
    # Set x axis limit.
    xp_min, xp_max = min_max_p[cp]
    plt.xlim(xp_min, xp_max)
    ax.locator_params(nbins=5)
    # Set minor ticks
    ax.minorticks_on()
    ax.tick_params(axis='y', which='major', labelsize=9)
    plt.xlabel(ld_p, fontsize=16)
    # Add textbox.
    text = (ld_p + '$ = {} \pm {}$').format(xp, e_xp)
    ob = offsetbox.AnchoredText(text, pad=0.1, loc=2, prop=dict(size=10))
    ob.patch.set(alpha=0.8)
    ax.add_artist(ob)
    # Plot scatter points over likelihood density map.
    cm = plt.cm.get_cmap('viridis')
    col_arr = [float(_) for _ in zip(*model_done[0])[ci]]
    SC = plt.scatter(zip(*model_done[0])[cp], model_done[1], marker='o',
                     c=col_arr, s=25, edgecolors='k',
                     lw=0.2, edgecolor='w', cmap=cm, zorder=2)
    # Set y axis limit.
    min_lik, med_lik = min(model_done[1]), np.median(model_done[1])
    min_y, max_y = min_lik - min_lik * 0.1, min(2.5 * min_lik, 1.2 * med_lik)
    plt.ylim(min_y, max_y)
    ax.locator_params(nbins=5)
    # Position colorbar.
    the_divider = make_axes_locatable(ax)
    color_axis = the_divider.append_axes("right", size="2%", pad=0.1)
    # Colorbar.
    cbar = plt.colorbar(SC, cax=color_axis)
    cbar.set_label(zlab, fontsize=12, labelpad=4)
    cbar.ax.tick_params(labelsize=8)
    # Best fit and errors.
    plt.axvline(x=xp, linestyle='--', color='red', zorder=4)
    if e_xp > 0.:
        # Plot error bars only if errors where assigned.
        plt.axvline(x=xp + e_xp, linestyle='--', color='blue', zorder=4)
        plt.axvline(x=xp - e_xp, linestyle='--', color='blue', zorder=4)


def plot(N, *args):
    '''
    Handle each plot separately.
    '''
    plt_map = {
        0: [pl_GA_lkl, 'GA likelihood evolution'],
        1: [pl_2_param_dens, 'age vs metallicity density map'],
        2: [pl_2_param_dens, 'distance vs extinction density map'],
        3: [pl_2_param_dens, 'age vs extinction density map'],
        4: [pl_2_param_dens, 'mass vs binarity density map'],
        5: [pl_lkl_scatt, 'z likelihood scatter'],
        6: [pl_lkl_scatt, 'age likelihood scatter'],
        7: [pl_lkl_scatt, 'extinction likelihood scatter'],
        8: [pl_lkl_scatt, 'distance likelihood scatter'],
        9: [pl_lkl_scatt, 'mass likelihood scatter'],
        10: [pl_lkl_scatt, 'binarity likelihood scatter']
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

import numpy as np
import matplotlib.pyplot as plt
from . prep_plots import xylabelsize, xytickssize, titlesize, legendsize,\
    grid_col, grid_ls, grid_lw


def pl_ad_test(gs, flag_ad_test, ad_cl, ad_fr, id_kinem):
    """
    """
    if flag_ad_test:

        def adPlot(ax, d1, d2, s):
            ax.tick_params(axis='both', which='major', labelsize=xytickssize)
            ax.set_title('(' + s + ')', fontsize=titlesize)
            ax.axes.yaxis.set_ticklabels([])
            plt.xlabel("A-D test (log)", fontsize=xylabelsize)
            plt.ylabel("N (log)", fontsize=xylabelsize)
            ax.grid(b=True, which='major', color=grid_col, linestyle=grid_ls,
                    lw=grid_lw, zorder=1)
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
            ax.set_yscale('log')
            ax.legend(fontsize=legendsize)

        ax = plt.subplot(gs[0:2, 0:2])
        adPlot(ax, ad_cl[0], ad_fr[0], 'phot')

        # Only plot if either parallax or PMs or radial velocities are defined
        pd_Plx, pd_PMRA, pd_RV = id_kinem[0], id_kinem[2], id_kinem[6]
        # Define first row depending on whether kinematic data was defined.
        k_flag = np.array([_ != 'n' for _ in (pd_Plx, pd_PMRA, pd_RV)]).any()
        if k_flag:
            ax = plt.subplot(gs[2:4, 0:2])
            adPlot(ax, ad_cl[1], ad_fr[1], 'plx+pm')


def pl_p_vals(
    ax, Ncl, Nf, prob_cl, kde_cl, kde_fr, x_cl, x_fr, x_over, y_over,
        reg_id):
    ax.tick_params(axis='both', which='major', labelsize=xytickssize)
    ax.set_title(
        r'$P_{{cl}}={:.2f}\;({})$'.format(prob_cl, reg_id), fontsize=titlesize)
    ax.axes.yaxis.set_ticklabels([])
    plt.xlabel('p-values', fontsize=xylabelsize)
    plt.ylabel('Density (norm)', fontsize=xylabelsize)
    ax.minorticks_on()
    ax.grid(b=True, which='major', color=grid_col, linestyle=grid_ls,
            lw=grid_lw, zorder=1)
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
    leg = ax.legend(handles, labels, numpoints=1, fontsize=legendsize)
    leg.get_frame().set_alpha(0.6)
    plt.gca().set_ylim(bottom=0)


def pl_ad_pvals_phot(gs, flag_ad_test, ad_cl_fr_p):
    if flag_ad_test:
        Ncl, Nf, prob_cl, kde_cl, kde_fr, x_cl, x_fr, x_over, y_over =\
            ad_cl_fr_p
        ax = plt.subplot(gs[0:2, 2:4])
        pl_p_vals(
            ax, Ncl, Nf, prob_cl, kde_cl, kde_fr, x_cl, x_fr, x_over, y_over,
            'phot')


def pl_ad_pvals_pk(gs, flag_ad_test, ad_cl_fr_pk, id_kinem):
    # Only plot if either parallax or PMs or radial velocities are defined
    pd_Plx, pd_PMRA, pd_RV = id_kinem[0], id_kinem[2], id_kinem[6]
    # Define first row depending on whether kinematic data was defined.
    k_flag = np.array([_ != 'n' for _ in (pd_Plx, pd_PMRA, pd_RV)]).any()
    if flag_ad_test and k_flag:
        Ncl, Nf, prob_cl, kde_cl, kde_fr, x_cl, x_fr, x_over, y_over =\
            ad_cl_fr_pk
        ax = plt.subplot(gs[2:4, 2:4])
        pl_p_vals(
            ax, Ncl, Nf, prob_cl, kde_cl, kde_fr, x_cl, x_fr, x_over, y_over,
            'Plx+PM')


def plot(N, *args):
    """
    Handle each plot separately.
    """

    plt_map = {
        0: [pl_ad_test, 'A-D test values'],
        1: [pl_ad_pvals_phot, 'photometric A-D pvalues'],
        2: [pl_ad_pvals_pk, 'photometric + kinem A-D pvalues']
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

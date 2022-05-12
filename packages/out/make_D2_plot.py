
import logging
import numpy as np
import matplotlib.pyplot as plt
# import matplotlib.gridspec as gridspec
from os.path import join
from . import add_version_plot
# from . import cornerPlot
# from . import prep_plots
# from . prep_plots import figsize_x, figsize_y, grid_x, grid_y
from .prep_plots import SigmaEllipse
from matplotlib.patches import Ellipse


def main(npd, pd, clp, td):
    """
    Make D2 block plots (corner plot).
    """

    try:
        import corner
    except ModuleNotFoundError:
        logging.warning(
            "corner module is not installed. Could not generate 'D2' plot")
        return

    fit_pars = clp['isoch_fit_params']

    fig = plotCorner(fit_pars, clp['ndim'], clp['varIdxs'])

    xf, yf = .02, .999
    add_version_plot.main(x_fix=xf, y_fix=yf)

    # Generate output file.
    fig.tight_layout()
    plt.savefig(join(
        npd['output_subdir'], str(npd['clust_name']) + '_D2_'
        + pd['best_fit_algor'] + npd['ext']))

    plt.clf()
    plt.close("all")

    # DEPRECATED 05/22

    # labels = [r'$z$', r'$log(age)$', r'$\beta$', r'$E_{(B-V)}$', r'$DR$',
    #           r'$R_v$', r'$(m-M)_o$']

    # plot_dict = {
    #     'metal-age':  [0, 1, 1, 2, 0, 1],
    #     'metal-beta': [0, 1, 2, 3, 0, 2],
    #     'metal-ext':  [0, 1, 3, 4, 0, 3],
    #     'metal-dr':   [0, 1, 4, 5, 0, 4],
    #     'metal-rv':   [0, 1, 5, 6, 0, 5],
    #     'metal-dist': [0, 1, 6, 7, 0, 6],
    #     #
    #     'age-beta':   [1, 2, 2, 3, 1, 2],
    #     'age-ext':    [1, 2, 3, 4, 1, 3],
    #     'age-dr':     [1, 2, 4, 5, 1, 4],
    #     'age-rv':     [1, 2, 5, 6, 1, 5],
    #     'age-dist':   [1, 2, 6, 7, 1, 6],
    #     #
    #     'beta-ext':   [2, 3, 3, 4, 2, 3],
    #     'beta-dr':    [2, 3, 4, 5, 2, 4],
    #     'beta-rv':    [2, 3, 5, 6, 2, 5],
    #     'beta-dist':  [2, 3, 6, 7, 2, 6],
    #     #
    #     'ext-dr':     [3, 4, 4, 5, 3, 4],
    #     'ext-rv':     [3, 4, 5, 6, 3, 5],
    #     'ext-dist':   [3, 4, 6, 7, 3, 6],
    #     #
    #     'dr-rv':      [4, 5, 5, 6, 4, 5],
    #     'dr-dist':    [4, 5, 6, 7, 4, 6],
    #     #
    #     'rv-dist':    [5, 6, 6, 7, 5, 6]
    # }
    # for key, val in plot_dict.items():

    #     # Limits for the 2-dens plots.
    #     min_max_p = prep_plots.param_ranges(
    #         td['fundam_params'], clp['varIdxs'], fit_pars['pars_chains'])
    #     min_max_p2 = prep_plots.p2_ranges(key, min_max_p)

    #     trace = fit_pars['mcmc_trace']
    #     args = [key, gs, labels, val, min_max_p2, clp['varIdxs'], trace]
    #     cornerPlot.plot(0, *args)

    # par_list = ['metal', 'age', 'beta', 'ext', 'dr', 'rv', 'dist']

    # # pl_param_pf: Parameters probability functions.
    # for p in par_list:
    #     args = [
    #         p, gs, labels, min_max_p, clp['varIdxs'],
    #         fit_pars['mean_sol'], fit_pars['map_sol'],
    #         fit_pars['median_sol'], fit_pars['mode_sol'],
    #         fit_pars['pardist_kde'], trace]
    #     cornerPlot.plot(1, *args)

    # # Generate output file.
    # fig.tight_layout()
    # plt.savefig(join(
    #     npd['output_subdir'], str(npd['clust_name']) + '_D2_'
    #     + pd['best_fit_algor'] + npd['ext']))

    # plt.clf()
    # plt.close("all")


def plotCorner(fit_pars, ndim, varIdxs):
    """
    """
    import corner

    labels = np.array([
        r'$z$', r'$log(age)$', r'$\beta$', r'$A_{V}$', r'$DR$',
        r'$R_V$', r'$(m-M)_o$'])

    # samples.shape = N_dim, N_pts
    samples = fit_pars['mcmc_trace']

    figure = corner.corner(samples.T, labels=labels[varIdxs])

    mean = np.array(fit_pars['mean_sol'])[varIdxs]
    median = np.array(fit_pars['median_sol'])[varIdxs]
    mode = np.array(fit_pars['mode_sol'])[varIdxs]

    # Extract the axes
    axes = np.array(figure.axes).reshape((ndim, ndim))

    # Loop over the diagonal
    for i in range(ndim):
        ax = axes[i, i]

        # Plot KDE.
        if fit_pars['pardist_kde'][i]:
            hist, _ = np.histogram(samples[i], bins=20)
            x_kde, par_kde = fit_pars['pardist_kde'][i]
            y_kde = (par_kde / max(par_kde)) * hist.max()
            ax.plot(x_kde, y_kde, ls='--', lw=1, c="grey")

        ax.axvline(mean[i], label='Mean', color="blue")
        ax.axvline(median[i], label='Median', color="red")
        ax.axvline(mode[i], label='Mode', color="cyan")
        pl, ph = np.percentile(samples[i], (16, 84))
        ax.axvline(x=ph, linestyle=':', color='orange', label="16-84th perc")
        ax.axvline(x=pl, linestyle=':', color='orange')
        if i == 0:
            handles, labels = ax.get_legend_handles_labels()
    ax_l = axes[0, 1]
    ax_l.axis('off')
    ax_l.legend(handles, labels, loc="center")

    # Loop over the histograms
    for yi in range(ndim):
        for xi in range(yi):
            ax = axes[yi, xi]
            mean_pos, width, height, theta = SigmaEllipse(np.array([
                samples[xi], samples[yi]]).T)
            # Plot 95% confidence ellipse.
            plt.scatter(
                mean_pos[0], mean_pos[1], marker='x', c='b', s=30, linewidth=2,
                zorder=4)
            ellipse = Ellipse(
                xy=mean_pos, width=width, height=height, angle=theta,
                edgecolor='k', fc='None', lw=.7, zorder=4)
            ax.add_patch(ellipse)

    return figure

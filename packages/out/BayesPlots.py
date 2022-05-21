
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from ..aux_funcs import reject_outliers
# from . cornerPlot import hist2d


def autocorr(gs, gsx, gsy, Nsteps, tau_autocorr, ESS):
    """
    Mean autocorrelation time.
    """
    plt.subplot(gs[gsy[0]:gsy[1], gsx[0]:gsx[1]])
    plt.plot(
        Nsteps * np.arange(tau_autocorr.size), tau_autocorr,
        label=r"$N_{{ESS}}\approx${:.0f}".format(ESS))
    plt.xlabel("Steps")
    plt.ylabel(r"$\hat{\tau}$")
    plt.legend()


def meanAF(gs, gsx, gsy, Nsteps, mean_afs):
    """
    Mean acceptance fraction.
    """
    plt.subplot(gs[gsy[0]:gsy[1], gsx[0]:gsx[1]])
    plt.plot(Nsteps * np.arange(mean_afs.size), mean_afs,
             label=r"$MAF\approx${:.2f}".format(mean_afs[-1]))
    plt.xlabel("Steps")
    plt.ylabel(r"$MAF$")
    plt.legend()


def traceplot(
    gs, gsx, gsy, mcmc_samples, stats_dict, Brn_prcnt, xylabel,
        xticks=True, tp_p=.5):
    """
    Chains traceplot.

    tp_p : percentage of the final portion of the trace that will be used to
    set the limits in the y axis.
    """
    ax = plt.subplot(gs[gsy[0]:gsy[1], gsx[0]:gsx[1]])
    N_tot = mcmc_samples.shape[0]
    plt.plot(mcmc_samples, c='k', lw=.8, ls='-', alpha=0.5)
    Nburn = Brn_prcnt * N_tot
    plt.axvline(x=Nburn, linestyle=':', color='r', zorder=4)
    plt.axhline(y=stats_dict['16th'], linestyle=':', color='orange', zorder=4)
    plt.axhline(y=stats_dict['median'], linestyle=':', color='blue', zorder=4)
    plt.axhline(y=stats_dict['84th'], linestyle=':', color='orange', zorder=4)
    if xticks is True:
        plt.xlabel("Steps")
    else:
        plt.xticks([])
    plt.ylabel(xylabel)
    # Use the last X% of the chains.
    N = int(mcmc_samples.shape[0] * tp_p)
    std = np.std(mcmc_samples[-N:])
    pmin, pmax = np.min(mcmc_samples[-N:]), np.max(mcmc_samples[-N:])
    ymin, ymax = pmin - std, pmax + std
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(0, N_tot)


def histogram(
    gs, gsx, gsy, mcmc_samples, stats_dict, mu_x_kde,
        xylabel, dec_places, xlims=None):
    """
    Parameter's distribution
    """
    # Histogram
    ax = plt.subplot(gs[gsy[0]:gsy[1], gsx[0]:gsx[1]])
    plt.xlabel(xylabel, fontsize=11)
    no_outlr = reject_outliers(mcmc_samples.flatten())
    # Obtain the bin edges.
    hist, bin_edges = np.histogram(no_outlr, bins=25)
    # Plot bars with the proper positioning, height, and width.
    plt.bar(
        (bin_edges[1:] + bin_edges[:-1]) * .5, hist / float(hist.max()),
        width=(bin_edges[1] - bin_edges[0]), color='grey', alpha=0.3)
    # Plot KDE.
    plt.plot(mu_x_kde[0], mu_x_kde[1] / max(mu_x_kde[1]), color='k', lw=1.5)

    # Mean
    plt.axvline(
        x=stats_dict['mean'], linestyle='--', color='blue', zorder=4,
        label=("Mean (" + dec_places + ")").format(stats_dict['mean']))
    # Median
    plt.axvline(
        x=stats_dict['median'], linestyle='--', color='green', zorder=4,
        label=("Median (" + dec_places + ")").format(stats_dict['median']))
    # Mode
    plt.axvline(
        x=stats_dict['mode'], linestyle='--', color='cyan', zorder=4,
        label=("Mode (" + dec_places + ")").format(stats_dict['mode']))

    # 16th and 84th percentiles.
    txt = "16-84th perc\n" +\
        (r"$(" + dec_places + ", " + dec_places + ")$").format(
            stats_dict['16th'], stats_dict['84th'])
    plt.axvline(
        x=stats_dict['16th'], linestyle=':', color='orange', zorder=4,
        label=txt)
    plt.axvline(x=stats_dict['84th'], linestyle=':', color='orange', zorder=4)

    if xlims is None:
        # MAD is robust to outliers
        mad = stats.median_abs_deviation(no_outlr)
        if not np.isnan(stats_dict['median']):
            cx = stats_dict['median']
        else:
            cx = stats_dict['mean']
        plt.xlim(cx - 4. * mad, cx + 4. * mad)
    else:
        plt.xlim(xlims[0], xlims[1])
    cur_ylim = ax.get_ylim()
    ax.set_ylim([0, cur_ylim[1]])
    plt.legend()


def twoParDens(
        gs, gsx, gsy, x_samples, y_samples, KP_Bys_x, KP_Bys_y, xylabel):
    """
    """
    ax = plt.subplot(gs[gsy[0]:gsy[1], gsx[0]:gsx[1]])

    x_flat, y_flat = x_samples.flatten(), y_samples.flatten()
    x_no_outlr = reject_outliers(x_flat)
    y_no_outlr = reject_outliers(y_flat)
    bin_edges = (np.linspace(x_no_outlr.min(), x_no_outlr.max(), 25),
                 np.linspace(y_no_outlr.min(), y_no_outlr.max(), 25))
    # hist2d(ax, x_samples, y_samples, bin_edges)
    H, xe, ye, _ = plt.hist2d(x_flat, y_flat, bin_edges, cmap='Greys')
    xbins = xe[:-1] + (xe[1] - xe[0]) / 2
    ybins = ye[:-1] + (ye[1] - ye[0]) / 2
    plt.contour(xbins, ybins, H.T, 6, colors='orange')

    plt.scatter(
        KP_Bys_x['median'], KP_Bys_y['median'], marker='x', c='green', s=50,
        zorder=5)
    plt.scatter(
        KP_Bys_x['mean'], KP_Bys_y['mean'], marker='x', c='blue', s=50,
        zorder=5)
    plt.scatter(
        KP_Bys_x['mode'], KP_Bys_y['mode'], marker='x', c='cyan', s=50,
        zorder=5)

    plt.xlabel(xylabel[0])
    plt.ylabel(xylabel[1])

    xl, yl = ax.get_xlim(), ax.get_ylim()
    # MAD is robust to outliers
    xmad = stats.median_abs_deviation(x_no_outlr)
    ymad = stats.median_abs_deviation(y_no_outlr)
    xsm, ysm = min(xmad, np.std(x_no_outlr)), min(ymad, np.std(y_no_outlr))
    xmed, ymed = np.median(x_no_outlr), np.median(y_no_outlr)
    x3s_min, x3s_max = xmed - 3. * xsm, xmed + 3. * xsm
    y3s_min, y3s_max = ymed - 3. * ysm, ymed + 3. * ysm
    plt.xlim(max(xl[0], x3s_min), min(xl[1], x3s_max))
    plt.ylim(max(yl[0], y3s_min), min(yl[1], y3s_max))

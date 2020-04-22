
import numpy as np
import matplotlib.pyplot as plt
from ..aux_funcs import reject_outliers
from . cornerPlot import hist2d
from . prep_plots import xylabelsize, xytickssize, legendsize


def histogram(
    gs, gsx, gsy, mcmc_samples, mu_kde_x, mu_kde, _16_50_84,
        xylabel):
    """
    Parameter's distribution
    """
    # Histogram
    ax = plt.subplot(gs[gsy[0]:gsy[1], gsx[0]:gsx[1]])
    ax.tick_params(axis='both', which='major', labelsize=xytickssize)
    plt.xlabel(xylabel, fontsize=11)
    no_outlr = reject_outliers(mcmc_samples.flatten())
    # Obtain the bin edges.
    hist, bin_edges = np.histogram(no_outlr, bins=25)
    # Plot bars with the proper positioning, height, and width.
    plt.bar(
        (bin_edges[1:] + bin_edges[:-1]) * .5, hist / float(hist.max()),
        width=(bin_edges[1] - bin_edges[0]), color='grey', alpha=0.3)
    # Plot KDE.
    plt.plot(mu_kde_x, mu_kde / max(mu_kde), color='k', lw=1.5)
    # Mean
    plt.axvline(
        x=_16_50_84[1], linestyle='--', color='blue', zorder=4,
        label=("Mean={:.3f}").format(_16_50_84[1]))
    # 16th and 84th percentiles.
    std = np.std(mcmc_samples.flatten())
    txt = "16-84th perc\n" + r"$\sigma={:.3f}$".format(std)
    plt.axvline(
        x=_16_50_84[0], linestyle=':', color='orange', zorder=4,
        label=txt)
    plt.axvline(x=_16_50_84[2], linestyle=':', color='orange', zorder=4)

    plt.xlim(max(-1., (_16_50_84[1]) - 4. * std), (_16_50_84[1]) + 4. * std)
    cur_ylim = ax.get_ylim()
    ax.set_ylim([0, cur_ylim[1]])
    plt.legend(fontsize=legendsize)


def traceplot(
    gs, gsx, gsy, mcmc_samples, _16_50_84, Brn_prcnt, xylabel,
        xticks=True):
    """
    Chains traceplot.
    """
    ax = plt.subplot(gs[gsy[0]:gsy[1], gsx[0]:gsx[1]])
    ax.tick_params(axis='both', which='major', labelsize=xytickssize)
    N_tot = mcmc_samples.shape[0]
    plt.plot(mcmc_samples, c='k', lw=.8, ls='-', alpha=0.5)
    Nburn = Brn_prcnt * N_tot
    plt.axvline(x=Nburn, linestyle=':', color='r', zorder=4)
    # 16th and 84th percentiles + median.
    plt.axhline(y=_16_50_84[0], linestyle=':', color='orange', zorder=4)
    plt.axhline(y=_16_50_84[1], linestyle=':', color='blue', zorder=4)
    plt.axhline(y=_16_50_84[2], linestyle=':', color='orange', zorder=4)
    if xticks is True:
        plt.xlabel("Steps", fontsize=xylabelsize)
    else:
        plt.xticks([])
    plt.ylabel(xylabel, fontsize=xylabelsize)
    # Use the last 10% of the chains.
    N = int(mcmc_samples.shape[0] * .1)
    std = np.std(mcmc_samples[-N:])
    pmin, pmax = np.min(mcmc_samples[-N:]), np.max(mcmc_samples[-N:])
    ymin, ymax = pmin - std, pmax + std
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(0, N_tot)


def autocorr(gs, gsx, gsy, Nsteps, tau_autocorr, ESS):
    """
    Mean autocorrelation time.
    """
    ax = plt.subplot(gs[gsy[0]:gsy[1], gsx[0]:gsx[1]])
    ax.tick_params(axis='both', which='major', labelsize=xytickssize)
    plt.plot(
        Nsteps * np.arange(tau_autocorr.size), tau_autocorr,
        label=r"$N_{{ESS}}\approx${:.0f}".format(ESS))
    plt.xlabel("Steps", fontsize=xylabelsize)
    plt.ylabel(r"$\hat{\tau}$", fontsize=xylabelsize)
    plt.legend(fontsize=legendsize)


def meanAF(gs, gsx, gsy, Nsteps, mean_afs):
    """
    Mean acceptance fraction.
    """
    ax = plt.subplot(gs[gsy[0]:gsy[1], gsx[0]:gsx[1]])
    ax.tick_params(axis='both', which='major', labelsize=xytickssize)
    plt.plot(Nsteps * np.arange(mean_afs.size), mean_afs)
    plt.xlabel("Steps", fontsize=xylabelsize)
    plt.ylabel(r"$MAF$", fontsize=xylabelsize)


def twoParDens(gs, gsx, gsy, KP_samples, KP_Bys_rc, KP_Bys_rt, xylabel):
    """
    """
    ax = plt.subplot(gs[gsy[0]:gsy[1], gsx[0]:gsx[1]])
    ax.tick_params(axis='both', which='major', labelsize=xytickssize)
    hist2d(ax, KP_samples[:, :, 0], KP_samples[:, :, 1])
    plt.scatter(KP_Bys_rc[1], KP_Bys_rt[1], marker='x', c='r', s=50, zorder=5)
    plt.xlabel(xylabel[0], fontsize=xylabelsize)
    plt.ylabel(xylabel[1], fontsize=xylabelsize)
    # xmed, xstd = np.median(KP_samples[:, :, 0]), np.std(KP_samples[:, :, 0])
    # ymed, ystd = np.median(KP_samples[:, :, 1]), np.std(KP_samples[:, :, 1])
    # plt.xlim(max(0.01, xmed - xstd), xmed + 2. * xstd)
    # plt.ylim(max(0.01, ymed - ystd), ymed + 2. * ystd)

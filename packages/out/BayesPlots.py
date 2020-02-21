
import numpy as np
import matplotlib.pyplot as plt
from ..aux_funcs import reject_outliers


def main(
    gs, gsx, gsy, mcmc_samples, mu_kde_x, mu_kde, _16_50_84, tau_autocorr,
        ESS, mean_afs, xylabel, Nsteps, Brn_prcnt):
    """
    """

    plt.style.use('seaborn-darkgrid')
    ax1 = plt.subplot(gs[gsy[0]:gsy[1], gsx[0]:gsx[1]])
    plt.xlabel(xylabel, fontsize=10)

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
        label=("Mean={:.3f} [mas]").format(_16_50_84[1]))
    # 16th and 84th percentiles.
    std = np.std(mcmc_samples.flatten())
    txt = "16-84th perc\n" + r"$\sigma={:.3f}$".format(std)
    plt.axvline(
        x=_16_50_84[0], linestyle=':', color='orange', zorder=4,
        label=txt)
    plt.axvline(x=_16_50_84[2], linestyle=':', color='orange', zorder=4)

    plt.xlim(max(-1., (_16_50_84[1]) - 4. * std), (_16_50_84[1]) + 4. * std)
    cur_ylim = ax1.get_ylim()
    ax1.set_ylim([0, cur_ylim[1]])
    plt.legend(fontsize='small')

    # Traceplot
    ax = plt.subplot(gs[2:3, 2:6])
    N_tot = mcmc_samples.shape[0]
    plt.plot(mcmc_samples, c='k', lw=.8, ls='-', alpha=0.5)
    Nburn = Brn_prcnt * N_tot
    plt.axvline(x=Nburn, linestyle=':', color='r', zorder=4)
    # 16th and 84th percentiles + median.
    plt.axhline(y=_16_50_84[0], linestyle=':', color='orange', zorder=4)
    plt.axhline(y=_16_50_84[1], linestyle=':', color='blue', zorder=4)
    plt.axhline(y=_16_50_84[2], linestyle=':', color='orange', zorder=4)
    plt.xlabel("Steps", fontsize=10)
    plt.ylabel(xylabel, fontsize=10)

    # Use the last 10% of the chains.
    N = int(mcmc_samples.shape[0] * .1)
    std = np.std(mcmc_samples[-N:])
    pmin, pmax = np.min(mcmc_samples[-N:]), np.max(mcmc_samples[-N:])
    ymin, ymax = pmin - std, pmax + std
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(0, N_tot)

    # Tau
    plt.subplot(gs[3:4, 2:4])
    plt.plot(
        Nsteps * np.arange(tau_autocorr.size), tau_autocorr,
        label=r"$N_{{ESS}}\approx${:.0f}".format(ESS))
    plt.xlabel("Steps", fontsize=10)
    plt.ylabel(r"$\hat{\tau}$", fontsize=10)
    plt.legend(fontsize='small')

    # MAF
    plt.subplot(gs[3:4, 4:6])
    plt.plot(10 * np.arange(mean_afs.size), mean_afs)
    plt.xlabel("Steps", fontsize=10)
    plt.ylabel(r"$MAF$", fontsize=10)

    plt.style.use('default')

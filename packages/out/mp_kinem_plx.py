
import matplotlib.pyplot as plt
import matplotlib.offsetbox as offsetbox
from matplotlib.colors import ListedColormap
import numpy as np
from astropy.stats import sigma_clipped_stats
from ..aux_funcs import reject_outliers


def plx_histo(
    gs, plx_offset, plx_clrg, plx_x_kde, kde_pl, plx_flrg,
        flag_no_fl_regs_i):
    '''
    Histogram for the distribution of parallaxes within the cluster region.
    '''
    ax = plt.subplot(gs[0:2, 0:2])
    ax.set_title('Offset applied: +{}'.format(plx_offset), fontsize=9)
    plt.xlabel('Plx [mas]', fontsize=12)
    plt.ylabel('N', fontsize=12)
    ax.minorticks_on()
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=.5,
            zorder=1)
    # Normalized histogram for cluster region.
    if len(plx_clrg) > 100:
        Nb = 100
    elif 50 < len(plx_clrg) <= 100:
        Nb = 50
    else:
        Nb = 25
    h_cl, _, _ = plt.hist(
        plx_clrg, Nb, density=True, zorder=4, color='#9aafd1',
        alpha=0.75, label=r"Cluster region ($N_{fit}$)")
    # Plot histogram for the parallaxes of the field regions.
    if not flag_no_fl_regs_i:
        plt.hist(
            plx_flrg, 100, density=True, zorder=1, color='#ef703e',
            alpha=0.5, label="Field regions")

    plt.plot(
        plx_x_kde, kde_pl / (max(kde_pl) / max(h_cl)), color='b', lw=1.,
        zorder=4)
    # Maximum KDE value.
    p_max_mas = plx_x_kde[np.argmax(kde_pl)]
    plt.axvline(x=p_max_mas, linestyle='--', color='b', lw=.7, zorder=5)
    plx_lt_zero = 100. * plx_clrg[plx_clrg < 0.].size / plx_clrg.size
    ob = offsetbox.AnchoredText(
        r"$Plx_{{max}}$={:.3f}".format(p_max_mas) + '\n' +
        r"$Plx<0 \rightarrow$ {:.1f}%".format(plx_lt_zero),
        pad=0.2, loc=1, prop=dict(size=9))
    ob.patch.set(alpha=0.85)
    ax.add_artist(ob)
    plt.xlim(
        max(-1., np.mean(plx_clrg) - 3. * np.std(plx_clrg)),
        min(3., np.mean(plx_clrg) + 3. * np.std(plx_clrg)))
    # Avoid showing the value 0.0 in the y axis.
    plt.ylim(0.01, plt.ylim()[1])
    ax.legend(fontsize='small', loc=7)


def plx_chart(gs, x_name, y_name, coord, cl_reg_fit, plx_Bys):
    '''
    Finding chart of cluster region with colors assigned according to the
    probabilities obtained and sizes according to parallaxes.
    '''
    ax = plt.subplot(gs[0:2, 4:6])
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=.5,
            zorder=1)
    # If RA is used, invert axis.
    if coord == 'deg':
        ax.invert_xaxis()
    # Set axis labels
    plt.xlabel('{} ({})'.format(x_name, coord), fontsize=12)
    plt.ylabel('{} ({})'.format(y_name, coord), fontsize=12)
    # Set minor ticks
    ax.minorticks_on()

    # Prepare data.
    x = np.array(list(zip(*cl_reg_fit))[1])
    y = np.array(list(zip(*cl_reg_fit))[2])
    mp = np.array(list(zip(*cl_reg_fit))[9])
    plx = np.array(list(zip(*list(zip(*cl_reg_fit))[7]))[0])
    msk = (~np.isnan(x)) & (~np.isnan(y)) & (~np.isnan(mp)) &\
        (~np.isnan(plx))
    x, y, mp, plx = x[msk], y[msk], mp[msk], plx[msk]
    ax.set_title('Cluster region (fit, N={})'.format(len(x)), fontsize=9)

    if plx_Bys.any():
        # Bayesian parallax value.
        p_max_mas = 1. / plx_Bys[1]
        # Distance to max value. Stars closer to the max value are larger.
        plx_d = 2. + 1. / (abs(plx - p_max_mas) + .1) ** 2.3

        # Re-arrange so stars closer to the max Plx value are on top
        plx_i = plx_d.argsort()
        x, y, mp, plx_d = x[plx_i], y[plx_i], mp[plx_i], plx_d[plx_i]

        # Color map
        cm = plt.cm.get_cmap('RdYlBu_r')
        # Get the colormap colors for my data
        my_cmap = cm(plt.Normalize(mp.min(), mp.max())(mp))
        # Set alpha
        alphas = (plx_d - plx_d.min()) / (plx_d.max() - plx_d.min())
        my_cmap[:, -1] = np.clip(alphas, a_min=.5, a_max=1.)
        # New colormap
        alpha_cmap = ListedColormap(my_cmap)
    else:
        alpha_cmap, plx_d = plt.cm.get_cmap('RdYlBu_r'), 20.

    # Plot stars selected to be used in the best fit process.
    plt.scatter(
        x, y, marker='o', c=mp, s=plx_d, edgecolors='black',
        cmap=alpha_cmap, lw=0.35, zorder=4)


def plx_vs_mag(
    gs, y_min_cmd, y_max_cmd, y_ax, mmag_plx_clp, mp_plx_clp, plx_clp,
        e_plx_clp, plx_flrg, mag_flrg, plx_Bys, plx_wa):
    '''
    Parallaxes versus main magnitude.
    '''
    ax = plt.subplot(gs[0:2, 2:4])
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=.5,
            zorder=1)

    ax.set_title("Plx clip " + r"$(med\pm2\sigma),\;N={}$".format(
        len(plx_clp)), fontsize=9)
    plt.xlabel('Plx [mas]', fontsize=12)
    plt.ylabel(y_ax, fontsize=12)
    # Set minor ticks
    ax.minorticks_on()

    cm = plt.cm.get_cmap('RdYlBu_r')
    # Plot stars selected to be used in the best fit process.
    plt.scatter(
        plx_clp, mmag_plx_clp, marker='o', c=mp_plx_clp, s=30,
        edgecolors='black', cmap=cm, lw=0.35, zorder=4)
    ax.errorbar(
        plx_clp, mmag_plx_clp, xerr=e_plx_clp, fmt='none', elinewidth=.35,
        ecolor='grey')
    if plx_flrg.any():
        # Field region values
        ax.scatter(
            plx_flrg, mag_flrg, marker='^', s=20, c='grey', linewidth=0.,
            alpha=.5, zorder=1)

    if plx_Bys.any():
        # Bayesian parallaxes in mas
        plx_16, plx_mean, plx_84 = 1. / plx_Bys
        # Bayesian parallaxes in pc
        pc_h, pc_m, pc_l = 1000. / plx_84, 1000. / plx_mean, 1000. / plx_16
        t0 = r"$Plx_{{Bay}} =$" + '\n'
        t0 += r"${:.0f}_{{{:.0f}}}^{{{:.0f}}}$ [pc]".format(
            pc_m, pc_l, pc_h)
        dm_50, dm_l, dm_h = 5. * np.log10(pc_m) - 5.,\
            5. * np.log10(pc_l) - 5., 5. * np.log10(pc_h) - 5.
        t1 = r"${:.3f}_{{{:.3f}}}^{{{:.3f}}}$ [mas]".format(
            plx_mean, plx_84, plx_16)
        t2 = r"${:.2f}_{{{:.2f}}}^{{{:.2f}}}$ [mag]".format(
            dm_50, dm_l, dm_h)
        txt = t0 + '\n' + t1 + '\n' + t2
        plt.axvline(
            x=plx_mean, linestyle='--', color='b', lw=1.2, zorder=5, label=txt)
    else:
        plt.axvline(
            x=np.mean(plx_clp), linestyle='--', color='b', lw=1.2, zorder=5,
            label=r"$Plx_{{mean}} = {:.3f}$".format(np.mean(plx_clp)))

    # Weighted average
    plt.axvline(
        x=plx_wa, linestyle='--', color='red', lw=.85, zorder=5,
        label=r"$Plx_{{wa}} = {:.3f}$".format(plx_wa))
    # Median
    plt.axvline(
        x=np.median(plx_clp), linestyle='--', color='k', lw=.85, zorder=5,
        label=r"$Plx_{{med}} = {:.3f}$".format(np.median(plx_clp)))

    cbar = plt.colorbar(pad=.01, fraction=.02, aspect=50)
    cbar.ax.tick_params(labelsize=7)
    cbar.set_label('MPs', size=8, labelpad=-20, y=1.07, rotation=0)

    ax.legend(fontsize='small', loc=0)
    min_plx, max_plx = np.min(plx_clp) - .2, np.max(plx_clp) + .2
    ax.axvspan(min_plx, 0., alpha=0.25, color='grey', zorder=1)
    plt.xlim(min_plx, max_plx)
    plt.ylim(max(mmag_plx_clp) + .5, min(mmag_plx_clp) - 1.)
    # ax.set_ylim(ax.get_ylim()[::-1])
    # plt.gca().invert_yaxis()


def pms_bys_params(
    gs, plx_bayes_flag_clp, plx_samples, plx_Bys, plx_mu_kde_x, plx_mu_kde,
        plx_tau_autocorr, mean_afs, plx_ess):
    """
    """
    if plx_bayes_flag_clp:
        plt.style.use('seaborn-darkgrid')
        ax1 = plt.subplot(gs[2:4, 0:2])
        plt.xlabel(r"$Plx_{Bay}$")

        plx_outlr = reject_outliers(plx_samples.flatten())
        # Obtain the bin values and edges using numpy
        hist, bin_edges = np.histogram(1. / plx_outlr, bins='auto')
        if len(bin_edges) > 25:
            hist, bin_edges = np.histogram(1. / plx_outlr, bins=20)
        # Plot bars with the proper positioning, height, and width.
        plt.bar(
            (bin_edges[1:] + bin_edges[:-1]) * .5, hist / float(hist.max()),
            width=(bin_edges[1] - bin_edges[0]), color='grey', alpha=0.3)
        # Plot KDE.
        plt.plot(plx_mu_kde_x, plx_mu_kde / max(plx_mu_kde), color='k', lw=1.5)
        # Mean
        plt.axvline(
            x=1. / plx_Bys[1], linestyle='--', color='blue', zorder=4,
            label=("Mean={:.3f} [mas]").format(1. / plx_Bys[1]))
        # 16th and 84th percentiles.
        std = np.std(1. / plx_samples.flatten())
        txt = "16-84th perc\n" + r"$\sigma={:.3f}$".format(std)
        plt.axvline(
            x=1. / plx_Bys[0], linestyle=':', color='orange', zorder=4,
            label=txt)
        plt.axvline(x=1. / plx_Bys[2], linestyle=':', color='orange', zorder=4)

        plt.xlim(
            max(-1., (1. / plx_Bys[1]) - 4. * std),
            (1. / plx_Bys[1]) + 4. * std)
        cur_ylim = ax1.get_ylim()
        ax1.set_ylim([0, cur_ylim[1]])
        plt.legend(fontsize='small')

        # Traceplot
        ax = plt.subplot(gs[2:3, 2:6])
        N_tot = plx_samples.shape[0]
        plt.plot(1. / plx_samples, c='k', lw=.8, ls='-', alpha=0.5)
        # HARDCODED: 25% of trace is burned
        plt.axvline(x=.25 * N_tot, linestyle=':', color='r', zorder=4)
        # 16th and 84th percentiles + median.
        plt.axhline(y=1. / plx_Bys[0], linestyle=':', color='orange', zorder=4)
        plt.axhline(y=1. / plx_Bys[1], linestyle=':', color='blue', zorder=4)
        plt.axhline(y=1. / plx_Bys[2], linestyle=':', color='orange', zorder=4)
        plt.xlabel("Steps")
        plt.ylabel(r"$Plx_{Bay}$")
        ax.set_xlim(0, N_tot)

        # Tau
        plt.subplot(gs[3:4, 2:4])
        # HARDCODED: store samples every 10 steps
        plt.plot(
            10 * np.arange(plx_tau_autocorr.size), plx_tau_autocorr,
            label=r"$N_{{ESS}}\approx${:.0f}".format(plx_ess))
        plt.xlabel("steps")
        plt.ylabel(r"$\hat{\tau}$")
        plt.legend(fontsize='small')

        # MAF
        plt.subplot(gs[3:4, 4:6])
        plt.plot(10 * np.arange(mean_afs.size), mean_afs)
        plt.xlabel("steps")
        plt.ylabel(r"$MAF$")

        plt.style.use('default')


def pms_vs_plx_mp_mag(
    gs, coord, y_ax, plx_bayes_flag_clp, plx_clp, plx_Bys, pmMP, pmRA_DE,
    pmDE, mmag_pm, pmRA_fl_DE, pmDE_fl, pm_Plx_cl, pm_Plx_fr, raPMrng,
        dePMrng):
    '''
    '''
    def axplot(
        gsi, j1, j2, yrang, xlabel, ylabel, cl_xdata, cl_ydata, cl_col,
            fr_xdata, fr_ydata, cmap, cbar_label, plx_dist):
        ax = plt.subplot(gs[gsi[0] + j1:gsi[1] + j1, gsi[2] + j2:gsi[3] + j2])
        ax.minorticks_on()
        ax.grid(b=True, which='major', color='gray', linestyle='--', lw=.5,
                zorder=1)
        ax.set_ylabel(ylabel, fontsize=12)
        ax.set_xlabel(xlabel, fontsize=12)

        plt.scatter(
            cl_xdata, cl_ydata, marker='o', c=cl_col, s=40, edgecolor='k',
            lw=.5, cmap=cmap, zorder=4, label="Cluster reg")
        cbar = plt.colorbar(pad=.01, fraction=.02, aspect=50)
        cbar.ax.tick_params(labelsize=7)
        cbar.set_label(cbar_label, size=8)
        if j2 == 0:
            cbar.ax.invert_yaxis()

        ax.scatter(
            fr_xdata, fr_ydata, marker='^', s=20, c='grey', linewidth=0.,
            alpha=.65, zorder=1, label="Field reg")

        plt.axvline(
            x=plx_dist[0], linestyle='--', color=plx_dist[1], lw=1.5,
            label=r"$Plx_{{{}}} = {:.3f}$ [mas]".format(
                plx_dist[2], plx_dist[0]))

        plx_mean, plx_median, plx_std = sigma_clipped_stats(cl_xdata)
        ax.set_xlim(
            max(-.05, plx_median - 2. * plx_std),
            min(max(cl_xdata) + .1, plx_median + 3. * plx_std))
        ax.set_ylim(yrang)

        if j1 == 0:
            plt.legend(fontsize='small')

    if plx_bayes_flag_clp:
        gsi = (4, 6, 0, 3)
        plx_dist = (1. / plx_Bys[1], 'g', 'Bay')
    else:
        gsi = (2, 4, 0, 3)
        plx_dist = (np.median(plx_clp), 'k', 'med')

    # Colored according to main mag
    if coord == 'deg':
        ylabel = r"$\mu_{{\alpha}} \, cos \delta \, \mathrm{[mas/yr]}$"
    else:
        ylabel = r"$\mu_{{\alpha}} \, \mathrm{[mas/yr]}$"
    xlabel = ""
    cmap = plt.cm.get_cmap('viridis')

    msk = np.argsort(mmag_pm)[::-1]
    pm_Plx_cl, pmRA_DE, pmDE, mmag_pm = pm_Plx_cl[msk], pmRA_DE[msk],\
        pmDE[msk], mmag_pm[msk]
    axplot(
        gsi, 0, 0, raPMrng, xlabel, ylabel, pm_Plx_cl, pmRA_DE, mmag_pm,
        pm_Plx_fr, pmRA_fl_DE, cmap, y_ax, plx_dist)
    xlabel = r"$Plx \, \mathrm{[mas]}$"
    ylabel = r"$\mu_{{\delta}} \, \mathrm{[mas/yr]}$"
    axplot(
        gsi, 2, 0, dePMrng, xlabel, ylabel, pm_Plx_cl, pmDE, mmag_pm,
        pm_Plx_fr, pmDE_fl, cmap, y_ax, plx_dist)

    # Colored according to MPs
    ylabel, xlabel = "", ""
    cmap = plt.cm.get_cmap('RdYlBu_r')
    axplot(
        gsi, 0, 3, raPMrng, xlabel, ylabel, pm_Plx_cl, pmRA_DE, pmMP,
        pm_Plx_fr, pmRA_fl_DE, cmap, 'MP', plx_dist)
    xlabel = r"$Plx \, \mathrm{[mas]}$"
    axplot(
        gsi, 2, 3, dePMrng, xlabel, ylabel, pm_Plx_cl, pmDE, pmMP, pm_Plx_fr,
        pmDE_fl, cmap, 'MP', plx_dist)


def plot(N, *args):
    '''
    Handle each plot separately.
    '''

    plt_map = {
        0: [plx_histo, 'Plx histogram'],
        1: [plx_chart, 'Plx chart'],
        2: [plx_vs_mag, 'Plx vs mag'],
        3: [pms_bys_params, 'PMs Bayes parameter'],
        4: [pms_vs_plx_mp_mag, 'PMs vs Plx, mag & MP colored']
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

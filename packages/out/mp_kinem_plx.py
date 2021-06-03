
import matplotlib.pyplot as plt
import matplotlib.offsetbox as offsetbox
from matplotlib.colors import ListedColormap
import numpy as np
# from astropy.stats import sigma_clipped_stats
from . import BayesPlots


def plx_histo(
    gs, plot_style, plx_offset, plx_clrg, plx_x_kde, kde_pl, plx_flrg,
        flag_no_fl_regs_i):
    """
    Histogram for the distribution of parallaxes within the cluster region.
    """
    ax = plt.subplot(gs[0:2, 0:2])
    ax.set_title('Offset applied: +{}'.format(plx_offset))
    plt.xlabel('Plx [mas]')
    plt.ylabel('N')
    if plot_style == 'asteca':
        ax.grid()
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
        r"$Plx_{{max}}$={:.3f} [mas]".format(p_max_mas) + '\n'
        + r"$Plx<0 \rightarrow$ {:.1f}%".format(plx_lt_zero), pad=0.2, loc=1)
    ob.patch.set(alpha=0.85)
    ax.add_artist(ob)
    plt.xlim(
        max(-1., np.mean(plx_clrg) - 3. * np.std(plx_clrg)),
        min(3., np.mean(plx_clrg) + 3. * np.std(plx_clrg)))
    # Avoid showing the value 0.0 in the y axis.
    plt.ylim(0.01, plt.ylim()[1])
    ax.legend(loc=7)


def plx_chart(gs, plot_style, x_name, y_name, coord, cl_reg_fit, plx_Bys):
    """
    Finding chart of cluster region with colors assigned according to the
    probabilities obtained and sizes according to parallaxes.
    """
    ax = plt.subplot(gs[0:2, 4:6])
    if plot_style == 'asteca':
        ax.grid()
    ax.invert_xaxis()
    # Set axis labels
    plt.xlabel('{} ({})'.format(x_name, coord))
    plt.ylabel('{} ({})'.format(y_name, coord))

    # Prepare data.
    x = np.array(list(zip(*cl_reg_fit))[1])
    y = np.array(list(zip(*cl_reg_fit))[2])
    mp = np.array(list(zip(*cl_reg_fit))[9])
    plx = np.array(list(zip(*list(zip(*cl_reg_fit))[7]))[0])
    msk = (~np.isnan(x)) & (~np.isnan(y)) & (~np.isnan(mp)) &\
        (~np.isnan(plx))
    x, y, mp, plx = x[msk], y[msk], mp[msk], plx[msk]
    ax.set_title(
        r'Cluster region ($N_{{fit}}$={})'.format(len(x)))

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
    gs, plot_style, y_min_cmd, y_max_cmd, y_ax, mmag_plx_clp,
        plx_clp, e_plx_clp, plx_flrg, mag_flrg, plx_Bys, plx_wa):
    """
    Parallaxes versus main magnitude.
    """
    ax = plt.subplot(gs[0:2, 2:4])
    if plot_style == 'asteca':
        ax.grid()

    # HARDCODED in plx_analysis: 3 sigma outlier rejection
    ax.set_title("Plx clip " + r"$(med\pm3\sigma,\;N={})$".format(
        len(plx_clp)))
    plt.xlabel('Plx [mas]')
    plt.ylabel(y_ax)

    # cm = plt.cm.get_cmap('RdYlBu_r')
    # Plot stars selected to be used in the best fit process.
    plt.scatter(
        plx_clp, mmag_plx_clp, alpha=.5, edgecolors='black', lw=0.3, zorder=4)
    # marker='o', s=30,
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
        plx_16, plx_median, plx_84, plx_mean, plx_mode = 1. / plx_Bys
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

    # cbar = plt.colorbar(pad=.01, fraction=.02, aspect=50)
    # cbar.ax.minorticks_off()
    # cbar.set_label('MPs', labelpad=-15, y=1.07, rotation=0)

    ax.legend(loc=0)
    min_plx, max_plx = np.min(plx_clp) - .2, np.max(plx_clp) + .2
    ax.axvspan(min_plx, 0., alpha=0.25, color='grey', zorder=1)
    plt.xlim(max(-1., min_plx), max_plx)
    plt.ylim(max(mmag_plx_clp) + .5, min(mmag_plx_clp) - 1.)
    # ax.set_ylim(ax.get_ylim()[::-1])
    # plt.gca().invert_yaxis()


def plx_bys_params(
    gs, plx_bayes_flag_clp, plx_burn, plx_samples, plx_Bayes_kde, plx_Bys,
        plx_tau_autocorr, mean_afs, plx_ess):
    """
    """
    if plx_bayes_flag_clp:

        dec_places = "{:.3f}"

        # Prepare data
        mcmc_samples = 1. / plx_samples
        _16_50_84_mean_mode = 1. / plx_Bys
        xylabel = "Plx"

        # HARDCODED in plx_analysis: store samples every 10 steps
        Nsteps = 10

        gsy, gsx = (2, 4), (0, 2)
        BayesPlots.histogram(
            gs, gsx, gsy, mcmc_samples, _16_50_84_mean_mode, plx_Bayes_kde,
            xylabel, dec_places)

        gsy, gsx = (2, 3), (2, 6)
        BayesPlots.traceplot(
            gs, gsx, gsy, mcmc_samples, _16_50_84_mean_mode, plx_burn,
            xylabel)

        gsy, gsx = (3, 4), (2, 4)
        BayesPlots.autocorr(
            gs, gsx, gsy, Nsteps, plx_tau_autocorr, plx_ess)

        gsy, gsx = (3, 4), (4, 6)
        BayesPlots.meanAF(gs, gsx, gsy, Nsteps, mean_afs)


# DEPRECATED 04/2021
# def pms_vs_plx_mp_mag(
#     gs, plot_style, coord, cosDE_flag, y_ax, plx_bayes_flag_clp, plx_clp,
#         plx_Bys, clreg_PMs, fregs_PMs, pm_Plx_cl, pm_Plx_fr, raPMrng, dePMrng):
#     """
#     """
#     def axplot(
#         gsi, j1, j2, yrang, xlabel, ylabel, cl_xdata, cl_ydata, cl_col,
#             fr_xdata, fr_ydata, cmap, cbar_label, plx_dist):
#         ax = plt.subplot(gs[gsi[0] + j1:gsi[1] + j1, gsi[2] + j2:gsi[3] + j2])
#         if plot_style == 'asteca':
#             ax.grid()
#         ax.set_ylabel(ylabel)
#         ax.set_xlabel(xlabel)

#         plt.scatter(
#             cl_xdata, cl_ydata, marker='o', c=cl_col, s=25, edgecolor='k',
#             lw=.2, cmap=cmap, zorder=4, label="Cluster reg")
#         cbar = plt.colorbar(pad=.01, fraction=.02, aspect=50)
#         cbar.ax.minorticks_off()
#         cbar.set_label(cbar_label)
#         if j2 == 0:
#             cbar.ax.invert_yaxis()

#         if fr_xdata.any():
#             ax.scatter(
#                 fr_xdata, fr_ydata, marker='^', s=20, c='grey', linewidth=0.,
#                 alpha=.65, zorder=1, label="Field reg")

#         plt.axvline(
#             x=plx_dist[0], linestyle='--', color=plx_dist[1], lw=1.5,
#             label=r"$Plx_{{{}}} = {:.3f}$ [mas] (median)".format(
#                 plx_dist[2], plx_dist[0]))

#         plx_median, plx_std = np.median(cl_xdata), np.std(cl_xdata)
#         min_plx = max(min(cl_xdata) - .1, plx_median - 2. * plx_std)
#         ax.set_xlim(
#             min_plx, min(max(cl_xdata) + .1, plx_median + 2. * plx_std))
#         ax.axvspan(min_plx, 0., alpha=0.25, color='grey', zorder=1)
#         ax.set_ylim(yrang)

#         if j1 == 0:
#             plt.legend(fontsize='small')

#     if plx_bayes_flag_clp:
#         gsi = (4, 6, 0, 3)
#         plx_dist = (1. / plx_Bys[1], 'g', 'Bay')
#     else:
#         gsi = (2, 4, 0, 3)
#         plx_dist = (np.median(plx_clp), 'k', 'med')

#     if cosDE_flag is False and coord == 'px':
#         ylabel = r"$\mu_{{\alpha}} \, \mathrm{[mas/yr]}$"
#     else:
#         ylabel = r"$\mu_{{\alpha}} \, cos \delta \, \mathrm{[mas/yr]}$"

#     # Colored according to main mag
#     xlabel = ""
#     cmap = plt.cm.get_cmap('viridis')
#     msk = np.argsort(clreg_PMs['mmag'])[::-1]

#     # Check field regions
#     fr_pmRA_m, pm_Plx_fr_m, fr_pmDE_m = np.array([]), np.array([]), []
#     if fregs_PMs['pmRA'].any():
#         pm_Plx_fr_m, fr_pmRA_m, fr_pmDE_m = pm_Plx_fr, fregs_PMs['pmRA'],\
#             fregs_PMs['pmDE']
#     axplot(
#         gsi, 0, 0, raPMrng, xlabel, ylabel, pm_Plx_cl[msk],
#         clreg_PMs['pmRA'][msk], clreg_PMs['mmag'][msk], pm_Plx_fr_m,
#         fr_pmRA_m, cmap, y_ax, plx_dist)
#     xlabel = r"$Plx \, \mathrm{[mas]}$"
#     ylabel = r"$\mu_{{\delta}} \, \mathrm{[mas/yr]}$"
#     axplot(
#         gsi, 2, 0, dePMrng, xlabel, ylabel, pm_Plx_cl[msk],
#         clreg_PMs['pmDE'][msk], clreg_PMs['mmag'][msk],
#         pm_Plx_fr_m, fr_pmDE_m, cmap, y_ax, plx_dist)

#     # Colored according to MPs
#     msk = np.argsort(clreg_PMs['MP'])
#     # Check field regions
#     fr_pmRA_m, pm_Plx_fr_m, fr_pmDE_m = np.array([]), np.array([]), []
#     if fregs_PMs['pmRA'].any():
#         pm_Plx_fr_m, fr_pmRA_m, fr_pmDE_m = pm_Plx_fr, fregs_PMs['pmRA'],\
#             fregs_PMs['pmDE']
#     ylabel, xlabel = "", ""
#     cmap = plt.cm.get_cmap('RdYlBu_r')
#     axplot(
#         gsi, 0, 3, raPMrng, xlabel, ylabel, pm_Plx_cl[msk],
#         clreg_PMs['pmRA'][msk], clreg_PMs['MP'][msk], pm_Plx_fr_m,
#         fr_pmRA_m, cmap, 'MP', plx_dist)
#     xlabel = r"$Plx \, \mathrm{[mas]}$"
#     axplot(
#         gsi, 2, 3, dePMrng, xlabel, ylabel, pm_Plx_cl[msk],
#         clreg_PMs['pmDE'][msk], clreg_PMs['MP'][msk], pm_Plx_fr_m,
#         fr_pmDE_m, cmap, 'MP', plx_dist)


def plot(N, *args):
    """
    Handle each plot separately.
    """

    plt_map = {
        0: [plx_histo, 'Plx histogram'],
        1: [plx_chart, 'Plx chart'],
        2: [plx_vs_mag, 'Plx vs mag'],
        3: [plx_bys_params, 'Plx Bayes parameter']
        # 4: [pms_vs_plx_mp_mag, 'PMs vs Plx, mag & MP colored']
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

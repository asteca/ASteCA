
import matplotlib.pyplot as plt
import matplotlib.offsetbox as offsetbox
from matplotlib.colors import Normalize
from matplotlib.patches import Ellipse
from matplotlib.colors import ListedColormap
import numpy as np


def plx_histo(
    gs, plx_flag, plx_clrg, plx_x_kde, kde_pl, plx_flrg,
        flag_no_fl_regs_i):
    '''
    Histogram for the distribution of parallaxes within the cluster region.
    '''
    if plx_flag:
        ax = plt.subplot(gs[0:2, 0:2])
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
            np.mean(plx_clrg) - 3. * np.std(plx_clrg),
            np.mean(plx_clrg) + 3. * np.std(plx_clrg))
        # Avoid showing the value 0.0 in the y axis.
        plt.ylim(0.01, plt.ylim()[1])
        ax.legend(fontsize='small', loc=7)


def plx_vs_mag(
    gs, y_min_cmd, y_max_cmd, y_ax, plx_flag, mmag_plx_clp, mp_plx_clp,
        plx_clp, e_plx_clp, plx_Bys, plx_wa):
    '''
    Parallaxes versus main magnitude.
    '''
    if plx_flag:
        ax = plt.subplot(gs[0:2, 2:4])
        ax.grid(b=True, which='major', color='gray', linestyle='--', lw=.5,
                zorder=1)

        ax.set_title("Plx clip " + r"$(med\pm2\sigma),\;N={}$".format(
            len(plx_clp)), fontsize=9)
        plt.xlabel('Plx [mas]', fontsize=12)
        plt.ylabel(y_ax, fontsize=12)
        # Set minor ticks
        ax.minorticks_on()

        cm = plt.cm.get_cmap('viridis')
        # Plot stars selected to be used in the best fit process.
        plt.scatter(
            plx_clp, mmag_plx_clp, marker='o', c=mp_plx_clp, s=30,
            edgecolors='black', cmap=cm, lw=0.35, zorder=4)
        ax.errorbar(
            plx_clp, mmag_plx_clp, xerr=e_plx_clp, fmt='none', elinewidth=.35,
            ecolor='grey')

        # Bayesian parallaxes in mas
        plx_16, plx_50, plx_84 = 1. / plx_Bys
        # Bayesian parallaxes in pc
        pc_h, pc_m, pc_l = 1000. / plx_84, 1000. / plx_50, 1000. / plx_16
        t0 = r"$Plx_{{Bay}} ={:.0f}_{{{:.0f}}}^{{{:.0f}}}\;pc$".format(
            pc_m, pc_l, pc_h)
        dm_50, dm_l, dm_h = 5. * np.log10(pc_m) - 5.,\
            5. * np.log10(pc_l) - 5., 5. * np.log10(pc_h) - 5.
        t0 += "\n" + r"$ {:.3f}\;[{:.2f}_{{{:.2f}}}^{{{:.2f}}}\;mag]$".format(
            plx_50, dm_50, dm_l, dm_h)
        plt.axvline(
            x=plx_50, linestyle='--', color='g', lw=1.2, zorder=5, label=t0)
        # Weighted average
        plt.axvline(
            x=plx_wa, linestyle='--', color='b', lw=.85, zorder=5,
            label=r"$Plx_{{wa}} = {:.3f}$".format(plx_wa))
        # Median
        plt.axvline(
            x=np.median(plx_clp), linestyle='--', color='k', lw=.85, zorder=5,
            label=r"$Plx_{{med}} = {:.3f}$".format(np.median(plx_clp)))

        cbar = plt.colorbar(pad=.01, fraction=.02, aspect=50)
        cbar.ax.tick_params(labelsize=7)
        # cbar.set_label('MP', size=8)

        ax.legend(fontsize='small', loc=0)
        min_plx, max_plx = np.min(plx_clp) - .2, np.max(plx_clp) + .2
        ax.axvspan(min_plx, 0., alpha=0.25, color='grey', zorder=1)
        plt.xlim(min_plx, max_plx)
        # ax.set_ylim(ax.get_ylim()[::-1])
        plt.gca().invert_yaxis()


def plx_chart(gs, plx_flag, x_name, y_name, coord, cl_reg_fit, plx_Bys):
    '''
    Finding chart of cluster region with colors assigned according to the
    probabilities obtained and sizes according to parallaxes.
    '''
    if plx_flag:
        ax = plt.subplot(gs[0:2, 4:6])
        ax.grid(b=True, which='major', color='gray', linestyle='--', lw=.5,
                zorder=1)
        ax.set_title('Cluster region (fit)'.format(), fontsize=9)

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

        if not np.isnan(plx_Bys[1]):
            # Bayesian parallax value.
            p_max_mas = 1. / plx_Bys[1]
            # Distance to max value. Stars closer to the max value are larger.
            plx_d = 2. + 1. / (abs(plx - p_max_mas) + .1) ** 2.3

            # Re-arrange so stars closer to the max Plx value are on top
            plx_i = plx_d.argsort()
            x, y, mp, plx_d = x[plx_i], y[plx_i], mp[plx_i], plx_d[plx_i]

            # Color map
            cm = plt.cm.get_cmap('viridis')
            # Get the colormap colors for my data
            my_cmap = cm(plt.Normalize(mp.min(), mp.max())(mp))
            # Set alpha
            alphas = (plx_d - plx_d.min()) / (plx_d.max() - plx_d.min())
            my_cmap[:, -1] = np.clip(alphas, a_min=.5, a_max=1.)
            # New colormap
            alpha_cmap = ListedColormap(my_cmap)
        else:
            alpha_cmap, plx_d = plt.cm.get_cmap('viridis'), 20.

        # Plot stars selected to be used in the best fit process.
        plt.scatter(
            x, y, marker='o', c=mp, s=plx_d, edgecolors='black',
            cmap=alpha_cmap, lw=0.35, zorder=4)


def pms_vpd(
    gs, coord, plx_flag, PM_flag, pmMP, pmRA_DE, e_pmRA_DE, pmDE, e_pmDE,
        pmRA_fl_DE, e_pmRA_fl_DE, pmDE_fl, e_pmDE_fl):
    '''
    '''
    if PM_flag:
        if plx_flag:
            ax = plt.subplot(gs[2:4, 0:2])
        else:
            ax = plt.subplot(gs[0:2, 0:2])
        ax.minorticks_on()

        if coord == 'deg':
            plt.xlabel(
                r"$\mu_{{\alpha}} \, cos \delta \, \mathrm{[mas/yr]}$",
                fontsize=12)
        else:
            plt.xlabel(r"$\mu_{{\alpha}} \, \mathrm{[mas/yr]}$", fontsize=12)
        plt.ylabel(r"$\mu_{{\delta}} \, \mathrm{[mas/yr]}$", fontsize=12)

        cmap = plt.cm.get_cmap('viridis')
        norm = Normalize(vmin=pmMP.min(), vmax=pmMP.max())

        ax.errorbar(
            pmRA_DE, pmDE, yerr=e_pmDE, xerr=e_pmRA_DE, fmt='none',
            elinewidth=.65, ecolor=cmap(norm(pmMP)), zorder=4)

        ax.errorbar(
            pmRA_fl_DE, pmDE_fl, yerr=e_pmDE_fl,
            xerr=e_pmRA_fl_DE, fmt='none', elinewidth=.35, ecolor='grey',
            zorder=1)

        RA_med, RA_std = np.median(pmRA_DE), np.std(pmRA_DE)
        DE_med, DE_std = np.median(pmDE), np.std(pmDE)
        plt.xlim(RA_med - 3. * RA_std, RA_med + 3. * RA_std)
        plt.ylim(DE_med - 3. * DE_std, DE_med + 3. * DE_std)


def pms_KDE_diag(
    gs, coord, plx_flag, PM_flag, pmRA_DE, pmDE, DE_pm, x_clpm, y_clpm,
        z_clpm, x_flpm, y_flpm, z_flpm, pmRA_Bys, pmDE_Bys):
    """
    KDE of PMs for cluster and field regions.
    """
    if PM_flag:
        if plx_flag:
            ax = plt.subplot(gs[2:4, 2:4])
        else:
            ax = plt.subplot(gs[0:2, 2:4])
        ax.set_title('Error-weighted KDE + Bayesian center', fontsize=9)
        ax.minorticks_on()
        ax.grid(b=True, which='major', color='gray', linestyle='--', lw=.5,
                zorder=1)

        if coord == 'deg':
            plt.xlabel(
                r"$\mu_{{\alpha}} \, cos \delta \, \mathrm{[mas/yr]}$",
                fontsize=12)
        else:
            plt.xlabel(r"$\mu_{{\alpha}} \, \mathrm{[mas/yr]}$", fontsize=12)
        plt.ylabel(r"$\mu_{{\delta}} \, \mathrm{[mas/yr]}$", fontsize=12)

        # Cluster region data
        RA_med, RA_std = np.median(pmRA_DE), np.std(pmRA_DE)
        DE_med, DE_std = np.median(pmDE), np.std(pmDE)
        ra_rang = (RA_med - 2. * RA_std, RA_med + 2. * RA_std)
        dec_rang = (DE_med - 2. * DE_std, DE_med + 2. * DE_std)

        # In place for #243
        import sys
        if sys.version_info[0] == 2:
            ax.hist2d(
                pmRA_DE, pmDE, range=(ra_rang, dec_rang), bins=25,
                normed=True, cmap=plt.cm.Greens)
        else:
            ax.hist2d(
                pmRA_DE, pmDE, range=(ra_rang, dec_rang), bins=25,
                density=True, cmap=plt.cm.Greens)

        # Cluster region contour
        max_i, max_j = np.unravel_index(z_clpm.argmax(), z_clpm.shape)
        CS = plt.contour(
            x_clpm, y_clpm, z_clpm, 5, colors='g', linewidths=1.,
            extent=(ra_rang[0], ra_rang[1], dec_rang[0], dec_rang[1]),
            zorder=3)
        CS.collections[0].set_label("Clust ({:.2f}, {:.2f})".format(
            x_clpm[max_i][max_j], y_clpm[max_i][max_j]))

        # Filed regions contours
        if z_flpm.any():
            max_i, max_j = np.unravel_index(z_flpm.argmax(), z_flpm.shape)
            CS = plt.contour(
                x_flpm, y_flpm, z_flpm, 10, colors='k', linewidths=.5,
                extent=(ra_rang[0], ra_rang[1], dec_rang[0], dec_rang[1]),
                zorder=2)
            CS.collections[0].set_label(
                "Field: ({:.2f}, {:.2f})".format(
                    x_flpm[max_i][max_j], y_flpm[max_i][max_j]))

        ellipse = Ellipse(
            xy=(pmRA_Bys[0], pmDE_Bys[0]), width=2 * pmRA_Bys[1],
            height=2. * pmDE_Bys[1], edgecolor='r', fc='None', lw=1.5,
            ls='--', zorder=4)
        ax.add_patch(ellipse)
        plt.scatter(
            pmRA_Bys[0], pmDE_Bys[0], c='r', marker='x', s=25, zorder=4,
            label=r"$[{:.3f}\pm{:.3f}, {:.3f}\pm{:.3f}]$".format(
                pmRA_Bys[0], pmRA_Bys[1], pmDE_Bys[0], pmDE_Bys[1]))

        plt.xlim(*ra_rang)
        plt.ylim(*dec_rang)
        plt.legend(fontsize='small')


def pms_vs_MP(gs, y_ax, plx_flag, PM_flag, pmMP, pm_dist_max, mmag_pm):
    '''
    '''
    if PM_flag:
        if plx_flag:
            ax = plt.subplot(gs[2:4, 4:6])
        else:
            ax = plt.subplot(gs[0:2, 4:6])
        ax.minorticks_on()
        ax.grid(b=True, which='major', color='gray', linestyle='--', lw=.5,
                zorder=1)

        ax.set_title("Distance to Bayesian median", fontsize=9)
        plt.xlabel("PM dist [mas/yr]", fontsize=12)
        plt.ylabel(y_ax, fontsize=12)

        cmap = plt.cm.get_cmap('viridis')

        plt.scatter(
            pm_dist_max, mmag_pm, marker='o', c=pmMP, s=30, edgecolors='black',
            cmap=cmap, lw=0.35, zorder=4)

        cbar = plt.colorbar(pad=.01, fraction=.02, aspect=50)
        cbar.ax.tick_params(labelsize=7)

        median_dst, std_dst = np.median(pm_dist_max), np.std(pm_dist_max)
        plt.xlim(-.05, median_dst * 2. * std_dst)
        plt.gca().invert_yaxis()


def plot(N, *args):
    '''
    Handle each plot separately.
    '''

    plt_map = {
        0: [plx_histo, 'Plx histogram'],
        1: [plx_chart, 'Plx chart'],
        2: [plx_vs_mag, 'Plx vs mag'],
        3: [pms_vpd, 'PMs vector point diagram'],
        4: [pms_KDE_diag, 'PMs KDE diagram'],
        5: [pms_vs_MP, 'PMs vs MP']
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

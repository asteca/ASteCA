
import matplotlib.pyplot as plt
import matplotlib.offsetbox as offsetbox
from matplotlib.colors import Normalize
from matplotlib.patches import Ellipse
from matplotlib.colors import ListedColormap
import numpy as np
from astropy.stats import sigma_clipped_stats


def plx_histo(gs, plx_clrg, plx_x_kde, kde_pl, plx_flrg, flag_no_fl_regs_i):
    '''
    Histogram for the distribution of parallaxes within the cluster region.
    '''
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
    gs, y_min_cmd, y_max_cmd, y_ax, mmag_plx_clp, mp_plx_clp, plx_clp,
        e_plx_clp, plx_Bys, plx_wa):
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


def plx_chart(gs, x_name, y_name, coord, cl_reg_fit, plx_Bys):
    '''
    Finding chart of cluster region with colors assigned according to the
    probabilities obtained and sizes according to parallaxes.
    '''
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
    gs, coord, plx_flag, pmMP, pmRA_DE, e_pmRA_DE, pmDE, e_pmDE, pmRA_fl_DE,
        e_pmRA_fl_DE, pmDE_fl, e_pmDE_fl, raPMrng, dePMrng):
    '''
    '''
    if plx_flag:
        ax = plt.subplot(gs[2:4, 0:2])
    else:
        ax = plt.subplot(gs[0:2, 0:2])
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

    cmap = plt.cm.get_cmap('viridis')
    norm = Normalize(vmin=pmMP.min(), vmax=pmMP.max())

    ax.errorbar(
        pmRA_DE, pmDE, yerr=e_pmDE, xerr=e_pmRA_DE, fmt='none',
        elinewidth=.65, ecolor=cmap(norm(pmMP)), zorder=4)

    ax.errorbar(
        pmRA_fl_DE, pmDE_fl, yerr=e_pmDE_fl, xerr=e_pmRA_fl_DE, fmt='none',
        elinewidth=.35, ecolor='grey', zorder=1)

    plt.xlim(raPMrng)
    plt.ylim(dePMrng)


def pms_KDE_diag(
    gs, coord, plx_flag, PM_cl_x, PM_cl_y, PM_cl_z, PM_fl_x, PM_fl_y, PM_fl_z,
    PMs_cl_cx, PMs_cl_cy, PMs_fl_cx, PMs_fl_cy, raPMrng, dePMrng, PMs_cent,
        PMs_width, PMs_height, PMs_theta):
    """
    KDE of PMs for cluster and field regions.
    """
    if plx_flag:
        ax = plt.subplot(gs[2:4, 2:4])
    else:
        ax = plt.subplot(gs[0:2, 2:4])
    ax.set_title('Error-weighted KDE + CI ellipse', fontsize=9)
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

    # Cluster region contour
    CS = plt.contour(
        PM_cl_x, PM_cl_y, PM_cl_z, 5, colors='g', linewidths=1.,
        extent=(raPMrng[0], raPMrng[1], dePMrng[0], dePMrng[1]),
        zorder=3)
    CS.collections[0].set_label(
        r"$KDE_{{max}},\, clust: ({:.3f}, {:.3f})$".format(
            PM_cl_x[PMs_cl_cx][PMs_cl_cy], PM_cl_y[PMs_cl_cx][PMs_cl_cy]))

    # Filed regions contours
    if PM_fl_z.any():
        CS = plt.contour(
            PM_fl_x, PM_fl_y, PM_fl_z, 10, colors='k', linewidths=.5,
            extent=(raPMrng[0], raPMrng[1], dePMrng[0], dePMrng[1]),
            zorder=2)
        CS.collections[0].set_label(
            r"$KDE_{{max}},\, field: ({:.3f}, {:.3f})$".format(
                PM_fl_x[PMs_fl_cx][PMs_fl_cy], PM_fl_y[PMs_fl_cx][PMs_fl_cy]))

    ellipse = Ellipse(
        xy=(PMs_cent[0], PMs_cent[1]), width=PMs_width,
        height=PMs_height, edgecolor='r', fc='None', lw=1.5,
        ls='--', zorder=4, angle=PMs_theta)
    ax.add_patch(ellipse)
    plt.scatter(
        PMs_cent[0], PMs_cent[1], c='r', marker='x', s=25, zorder=4,
        label=(r"$cent=[{:.3f},{:.3f}]$" + '\n' +
               r"$(w,h)=[{:.3f},{:.3f}]$").format(
                   PMs_cent[0], PMs_cent[1], PMs_width, PMs_height))

    plt.xlim(*raPMrng)
    plt.ylim(*dePMrng)
    plt.legend(fontsize='small')


def pms_vs_MP(gs, y_ax, plx_flag, pmMP, pm_dist_max, mmag_pm):
    '''
    '''
    if plx_flag:
        ax = plt.subplot(gs[2:4, 4:6])
    else:
        ax = plt.subplot(gs[0:2, 4:6])
    ax.minorticks_on()
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=.5,
            zorder=1)

    ax.set_title("Distance to KDE's cluster max", fontsize=9)
    plt.xlabel("PM dist [mas/yr]", fontsize=12)
    plt.ylabel(y_ax, fontsize=12)

    cmap = plt.cm.get_cmap('viridis')

    plt.scatter(
        pm_dist_max, mmag_pm, marker='o', c=pmMP, s=30, edgecolors='black',
        cmap=cmap, lw=0.35, zorder=4)

    cbar = plt.colorbar(pad=.01, fraction=.02, aspect=50)
    cbar.ax.tick_params(labelsize=7)

    plx_mean, plx_median, plx_std = sigma_clipped_stats(pm_dist_max)
    plt.xlim(-.05, plx_median + 2. * plx_std)
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

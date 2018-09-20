
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.pyplot import cm
from .prep_plots import CIEllipse


def pl_2_param_dens(_2_params, gs, min_max_p2, varIdxs, mcmc_trace):
    '''
    Parameter vs parameters density map.
    '''
    plot_dict = {
        'metal-age': [0, 2, 2, 4, 'r', 'GnBu', 0, 1],
        'metal-ext': [0, 2, 4, 6, 'b', 'YlOrBr', 0, 2],
        'metal-dist': [0, 2, 6, 8, 'r', 'GnBu', 0, 3],
        'metal-mass': [0, 2, 8, 10, 'b', 'YlOrBr', 0, 4],
        'metal-binar': [0, 2, 10, 12, 'r', 'GnBu', 0, 5],
        'age-ext': [2, 4, 4, 6, 'b', 'YlOrBr', 1, 2],
        'age-dist': [2, 4, 6, 8, 'r', 'GnBu', 1, 3],
        'age-mass': [2, 4, 8, 10, 'b', 'YlOrBr', 1, 4],
        'age-binar': [2, 4, 10, 12, 'r', 'GnBu', 1, 5],
        'ext-dist': [4, 6, 6, 8, 'r', 'GnBu', 2, 3],
        'ext-mass': [4, 6, 8, 10, 'b', 'YlOrBr', 2, 4],
        'ext-binar': [4, 6, 10, 12, 'r', 'GnBu', 2, 5],
        'dist-mass': [6, 8, 8, 10, 'b', 'YlOrBr', 3, 4],
        'dist-binar': [6, 8, 10, 12, 'r', 'GnBu', 3, 5],
        'mass-binar': [8, 10, 10, 12, 'r', 'GnBu', 4, 5]
    }

    labels = ['$z$', '$log(age)$', '$E_{(B-V)}$', '$(m-M)_o$',
              '$M\,(M_{{\odot}})$', '$b_{frac}$']

    gs_x1, gs_x2, gs_y1, gs_y2, cp, d_map, mx, my = plot_dict[_2_params]
    x_label, y_label = labels[mx], labels[my]

    ax = plt.subplot(gs[gs_y1:gs_y2, gs_x1:gs_x2])

    # To specify the number of ticks on both or any single axes
    ax.locator_params(nbins=5)
    if gs_x1 == 0:
        plt.ylabel(y_label, fontsize=16)
        plt.yticks(rotation=45)
    else:
        ax.tick_params(labelleft=False)
    if gs_y2 == 12:
        plt.xlabel(x_label, fontsize=16)
        plt.xticks(rotation=45)
    else:
        ax.tick_params(labelbottom=False)
    plt.minorticks_on()

    if mx in varIdxs and my in varIdxs:
        mx_model, my_model = varIdxs.index(mx), varIdxs.index(my)

        h2d, xbins, ybins = plt.hist2d(
            mcmc_trace[mx_model], mcmc_trace[my_model], bins=30,
            cmap=plt.get_cmap(d_map),
            range=None, zorder=2)[:-1]
        plt.contour(
            h2d.transpose(), 5,
            extent=[xbins.min(), xbins.max(), ybins.min(), ybins.max()],
            colors='#551a8b', linewidths=0.5, zorder=3)

        mean_pos, width, height, theta = CIEllipse(np.array([
            mcmc_trace[mx_model], mcmc_trace[my_model]]).T)
        # Plot 2 sigma ellipse.
        plt.scatter(
            mean_pos[0], mean_pos[1], marker='x', c=cp, s=30, linewidth=2,
            zorder=4)
        ellipse = Ellipse(xy=mean_pos, width=width, height=height, angle=theta,
                          edgecolor=cp, fc='None', lw=1., zorder=4)
        ax.add_patch(ellipse)

    xp_min, xp_max, yp_min, yp_max = min_max_p2
    ax.set_xlim([xp_min, xp_max])
    ax.set_ylim([yp_min, yp_max])
    ax.set_aspect('auto')


def pl_param_pf(
    par_name, gs, min_max_p, cp_r, cp_e, varIdxs, map_sol,
        model_done):
    '''
    Parameter posterior plot.
    '''
    plot_dict = {
        'metal': [0, 2, 0, 2, 0], 'age': [2, 4, 2, 4, 1],
        'ext': [4, 6, 4, 6, 2], 'dist': [6, 8, 6, 8, 3],
        'mass': [8, 10, 8, 10, 4], 'binar': [10, 12, 10, 12, 5]
    }

    gs_x1, gs_x2, gs_y1, gs_y2, cp = plot_dict[par_name]

    labels = [r'$z$', r'$\log(age)$', r'$E_{{(B-V)}}$', r'$(m-M)_o$',
              r'$M\,(M_{{\odot}})$', r'$b_{{frac}}$']
    frm = ["{:.4f}", "{:.2f}", "{:.2f}", "{:.2f}", "{:.0f}", "{:.2f}"]

    ld_p = labels[cp]
    p = frm[cp]

    ax = plt.subplot(gs[gs_y1:gs_y2, gs_x1:gs_x2])
    plt.title(ld_p, fontsize=10)
    # Parameter values and errors.
    xp, e_xp = map(float, [cp_r[cp], cp_e[cp]])
    # Set x axis limit.
    xp_min, xp_max = min_max_p[cp]
    ax.set_xlim(xp_min, xp_max)
    ax.locator_params(nbins=5)
    # Set minor ticks
    ax.minorticks_on()
    if cp == 5:
        plt.xlabel(ld_p, fontsize=16)
        plt.xticks(rotation=45)
    # else:
    #     ax.tick_params(labelbottom=False)
    ax.tick_params(axis='y', which='major', labelleft=False)

    if cp in varIdxs:
        c_model = varIdxs.index(cp)

        # Define KDE limits.
        x_rang = .1 * (xp_max - xp_min)
        x_kde = np.mgrid[xp_min - x_rang:xp_max + x_rang:100j]
        # Use a larger Scott bandwidth
        bw = 1.5 * len(model_done[c_model]) ** (-1. / (len(varIdxs) + 4))
        kernel_cl = stats.gaussian_kde(model_done[c_model], bw_method=bw)
        # KDE for plotting.
        try:
            kde = np.reshape(kernel_cl(x_kde).T, x_kde.shape)
            plt.plot(x_kde, kde / max(kde), color='k', lw=1.5)
            # Mode (using KDE)
            x_mode = x_kde[np.argmax(kde)]
            plt.axvline(
                x=x_mode, linestyle='--', color='cyan', zorder=4,
                label=("Mode (" + p + ")").format(x_mode))
        except (FloatingPointError, UnboundLocalError):
            pass

        # Obtain the bin values and edges using numpy
        hist, bin_edges = np.histogram(model_done[c_model], bins=25)
        # Plot bars with the proper positioning, height, and width.
        plt.bar(
            (bin_edges[1:] + bin_edges[:-1]) * .5, hist / float(hist.max()),
            width=(bin_edges[1] - bin_edges[0]), color='grey', alpha=0.3)

        # Mean
        # x_mean = np.mean(model_done[c_model])
        plt.axvline(
            x=xp, linestyle='--', color='blue', zorder=4,
            label=("Mean (" + p + ")").format(xp))
        # MAP
        plt.axvline(
            x=map_sol[cp], linestyle='--', color='red', zorder=4,
            label=("MAP (" + p + ")").format(map_sol[cp]))
        # Median
        pm = np.percentile(model_done[c_model], 50)
        plt.axvline(
            x=pm, linestyle=':', color='green', zorder=4,
            label=("Median (" + p + ")").format(pm))
        #  16th and 84th percentiles (1 sigma) around median.
        ph = np.percentile(model_done[c_model], 84) - pm
        pl = pm - np.percentile(model_done[c_model], 16)
        if ph > 0. and pm > 0.:
            # Plot error bars only if errors where assigned.
            plt.axvline(
                x=pm + ph, linestyle=':', color='orange', zorder=4,
                label=("16-84th perc\n(" + p + ", " + p + ")").format(
                    pm - pl, pm + ph))
            plt.axvline(
                x=pm - pl, linestyle=':', color='orange', zorder=4)

        cur_ylim = ax.get_ylim()
        ax.set_ylim([0, cur_ylim[1]])
        plt.legend(fontsize='small')


def pl_pdf_half(dummy, gs, mcmc_halves):
    '''
    1st and 2nd halves intersection.
    '''
    p_names = ['metal', 'age', 'ext', 'dist', 'mass', 'binar']
    ax = plt.subplot(gs[0:2, 2:4])
    ax.set_title("1 - (1st vs 2nd halves intersection)")
    ax.bar(p_names, 1. - np.array(mcmc_halves))
    # ax.set_ylim(0, 1.01)


def pl_MAP_lkl(dummy, gs, map_lkl, map_lkl_final):
    '''
    Evolution of MAP likelihood values.
    '''
    ax = plt.subplot(gs[0:2, 4:6])
    x, y = list(zip(*map_lkl))
    ax.plot(x, y, label=r"$L_{{min}}={:.1f}$".format(map_lkl_final))
    plt.xlabel("steps", fontsize=14)
    plt.ylabel("Lkl (MAP)", fontsize=14)
    ax.legend(fontsize='small', loc=0, handlelength=0.)


def pl_MAF(dummy, gs, maf_steps):
    '''
    Evolution of MAF values.
    '''
    ax = plt.subplot(gs[2:4, 4:6])
    # ax.set_title("1 - (1st vs 2nd halves intersection)")
    x, y = list(zip(*maf_steps))
    ax.plot(x, y, label=r"$MAF={:.3f}$".format(maf_steps[-1][1]))
    plt.xlabel("steps", fontsize=14)
    plt.ylabel("MAF", fontsize=14)
    plt.axhline(y=.25, color='grey', ls=':', lw=1.2, zorder=4)
    plt.axhline(y=.5, color='grey', ls=':', lw=1.2, zorder=4)
    ax.legend(fontsize='small', loc=0, handlelength=0.)


def pl_param_chain(
    par_name, gs, cp_r, min_max_p, nwalkers, nburn, nsteps, emcee_a,
        mcmc_elapsed, varIdxs, pre_bi, post_bi, autocorr_time,
        max_at_5c, min_at_5c, pymc3_ess):
    '''
    Parameter sampler chain.
    '''
    plot_dict = {
        'metal': [8, 12, 0, 1, 0], 'age': [8, 12, 1, 2, 1],
        'ext': [8, 12, 2, 3, 2], 'dist': [8, 12, 3, 4, 3],
        'mass': [8, 12, 4, 5, 4], 'binar': [8, 12, 5, 6, 5]
    }

    labels = [r'$z$', r'$\log(age)$', r'$E_{{(B-V)}}$', r'$(m-M)_o$',
              r'$M\,(M_{{\odot}})$', r'$b_{{frac}}$']

    gs_x1, gs_x2, gs_y1, gs_y2, cp = plot_dict[par_name]
    ax = plt.subplot(gs[gs_y1:gs_y2, gs_x1:gs_x2])
    if cp == 0:
        m, s = divmod(mcmc_elapsed, 60)
        h, m = divmod(m, 60)
        plt.title(
            "steps={:.0f}, chains={:.0f}, a={:.2f} | {:.0f}h{:.0f}m".format(
                nsteps, nwalkers, emcee_a, h, m), fontsize=10)
    if cp == 5:
        plt.xlabel("Steps")
    else:
        ax.tick_params(labelbottom=False)
    plt.ylabel(labels[cp])
    ax.set_xlim(0, nburn + nsteps)

    if cp in varIdxs:
        c_model = varIdxs.index(cp)

        # Worst chains
        color = iter(cm.rainbow(np.linspace(0, 1, len(max_at_5c[0]))))
        pre_bi_max_at = pre_bi[c_model][max_at_5c[c_model]]
        post_bi_max_at = post_bi[c_model][max_at_5c[c_model]]
        for w1, w2 in zip(*[pre_bi_max_at, post_bi_max_at]):
            # Burn-in stage
            plt.plot(range(nburn), w1, c='grey', lw=.5, alpha=0.5)
            # Post burn-in.
            c = next(color)
            plt.plot(np.arange(nburn, nburn + nsteps), w2, c=c, lw=.8,
                     ls='--', alpha=0.5)
        # Best chains
        color = iter(cm.rainbow(np.linspace(0, 1, len(min_at_5c[0]))))
        pre_bi_min_at = pre_bi[c_model][min_at_5c[c_model]]
        post_bi_min_at = post_bi[c_model][min_at_5c[c_model]]
        for w1, w2 in zip(*[pre_bi_min_at, post_bi_min_at]):
            # Burn-in stage
            plt.plot(range(nburn), w1, c='grey', lw=.5, alpha=0.5)
            # Post burn-in.
            c = next(color)
            plt.plot(np.arange(nburn, nburn + nsteps), w2, c=c, lw=.8,
                     alpha=0.5)

        plt.axhline(
            y=float(cp_r[cp]), color='k', ls='--', lw=1.2, zorder=4,
            label=r"$\tau={:.0f}\;(\hat{{n}}_{{eff}}={:.0f})$".format(
                autocorr_time[c_model], pymc3_ess[c_model]))
        ax.set_ylim(min_max_p[c_model][0], min_max_p[c_model][1])
        ax.legend(fontsize='small', loc=0, handlelength=0.)


def pl_mESS(dummy, gs, mESS, minESS, minESS_epsilon):
    '''
    mESS plot.
    '''
    ax = plt.subplot(gs[4:6, 6:8])
    plt.xlabel(r"$CI\;(=1-\alpha)$", fontsize=14)
    plt.ylabel(r"$\epsilon$", fontsize=16)
    ax.set_xlim(.01, 1.01)

    plt.plot(
        1. - np.array(minESS_epsilon[0]), minESS_epsilon[1],
        label="minESS ({:.0f})".format(minESS))
    plt.plot(
        1. - np.array(minESS_epsilon[0]), minESS_epsilon[2],
        label="mESS ({:.0f})".format(mESS))
    ax.legend(fontsize='small', loc=0)


def pl_tau(dummy, gs, N_steps_conv, N_conv, tol_conv, tau_index, tau_autocorr):
    '''
    Tau vs steps plot.
    '''
    ax = plt.subplot(gs[6:8, 8:10])
    plt.title(r"$N_{{conv}}={:.0f}, tol_{{conv}}={:.2f}$".format(
        N_conv, tol_conv))
    plt.xlabel("steps", fontsize=14)
    plt.ylabel(r"mean $\hat{\tau}$", fontsize=14)

    n = N_steps_conv * np.arange(1, tau_index + 1)
    plt.plot(n, n / 50, "--b", label="N/50")
    plt.plot(n, n / 100, "--g", label="N/100")
    plt.plot(n, tau_autocorr)
    plt.xlim(0, n.max())
    try:
        plt.ylim(
            0, np.nanmax(tau_autocorr) + 0.1 *
            (np.nanmax(tau_autocorr) - np.nanmin(tau_autocorr)))
    except ValueError:
        print("  WARNING: no mean autocorrelation values to plot.")
    ax.legend(fontsize='small', loc=0)


def pl_lags(dummy, gs, varIdxs, emcee_acorf):
    '''
    lags plot.
    '''
    ax = plt.subplot(gs[6:8, 10:12])
    plt.xlabel("Lag", fontsize=14)
    plt.ylabel("Autocorrelation", fontsize=14)

    plot_dict = ['metal', 'age', 'ext', 'dist', 'mass', 'binar']
    for cp, p in enumerate(emcee_acorf):
        if cp in varIdxs:
            c_model = varIdxs.index(cp)
            plt.plot(
                range(len(p)), p, lw=.8, alpha=0.5,
                label="{}".format(plot_dict[c_model]))
    ax.legend(fontsize='small', loc=0)


def pl_GW(dummy, gs, varIdxs, geweke_z):
    '''
    Geweke plot.
    '''
    ax = plt.subplot(gs[8:10, 10:12])
    ax.set_title("Geweke", fontsize=10)
    plt.xlabel("First iteration in segment", fontsize=14)
    plt.ylabel("z-score", fontsize=14)
    plt.axhline(y=2., color='grey', ls=':', lw=1.2, zorder=4)
    plt.axhline(y=-2., color='grey', ls=':', lw=1.2, zorder=4)

    plot_dict = ['metal', 'age', 'ext', 'dist', 'mass', 'binar']
    ymin, ymax = [], []
    for cp, p in enumerate(geweke_z):
        if cp in varIdxs:
            c_model = varIdxs.index(cp)
            idx, zscore = list(zip(*p))
            plt.plot(
                idx, zscore, ls="-.", linewidth=.7,
                label="{}".format(plot_dict[c_model]))
            ymin.append(np.nanmin(zscore))
            ymax.append(np.nanmax(zscore))

    if ymin and ymax:
        ax.set_ylim(max(-2.1, min(ymin)), min(2.1, max(ymax)))
    ax.legend(fontsize='small', loc=0)


def plot(N, *args):
    '''
    Handle each plot separately.
    '''
    plt.style.use('seaborn-darkgrid')
    plt_map = {
        0: [pl_2_param_dens, args[0] + ' density map'],
        1: [pl_param_pf, args[0] + ' probability function'],
        2: [pl_pdf_half, ' 1st and 2nd halfs of pdf'],
        3: [pl_MAP_lkl, ' MAP likelihood values'],
        4: [pl_MAF, ' MAF vs steps'],
        5: [pl_param_chain, args[0] + ' sampler chain'],
        6: [pl_tau, args[0]],
        7: [pl_mESS, args[0]],
        8: [pl_lags, args[0]],
        9: [pl_GW, args[0]]
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
    plt.style.use('default')

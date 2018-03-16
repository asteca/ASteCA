
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.offsetbox as offsetbox
from matplotlib.patches import Ellipse
from matplotlib.pyplot import cm


def pl_2_param_dens(_2_params, gs, min_max_p2, cp_r, cp_e, varIdxs, model_done):
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
    # Parameter values and errors.
    xp, e_xp = map(float, [cp_r[mx], cp_e[mx]])
    yp, e_yp = map(float, [cp_r[my], cp_e[my]])
    # # Axis limits.
    xp_min, xp_max = min_max_p2[mx]
    yp_min, yp_max = min_max_p2[my]

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
        # Plot best fit point (MAP, mean, median?). # TODO
        plt.scatter(xp, yp, marker='x', c=cp, s=30, linewidth=2, zorder=4)

        mx_model, my_model = varIdxs.index(mx), varIdxs.index(my)

        h2d, xbins, ybins = plt.hist2d(
            model_done[mx_model], model_done[my_model], bins=20,
            cmap=plt.get_cmap(d_map),
            range=[[xp_min, xp_max], [yp_min, yp_max]], zorder=2)[:-1]
        plt.contour(
            h2d.transpose(), 5,
            extent=[xbins.min(), xbins.max(), ybins.min(), ybins.max()],
            colors='#551a8b', linewidths=0.5, zorder=3)

        # Assume uncertainties in both dimensions are defined.
        plt.gca()
        ellipse = Ellipse(xy=(xp, yp), width=2 * e_xp, height=2 * e_yp,
                          edgecolor=cp, fc='None', lw=1., zorder=4)
        ax.add_patch(ellipse)

    ax.set_xlim(xp_min, xp_max)
    ax.set_ylim(yp_min, yp_max)
    ax.set_aspect('auto')


def pl_param_pf(par_name, gs, min_max_p2, cp_r, cp_e, varIdxs, model_done):
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
    # Parameter values and errors.
    xp, e_xp = map(float, [cp_r[cp], cp_e[cp]])
    # Set x axis limit.
    xp_min, xp_max = min_max_p2[cp]
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
    # Add textbox.
    ob = offsetbox.AnchoredText(ld_p, pad=0.1, loc=9, prop=dict(size=10))
    ob.patch.set(alpha=0.8)
    ax.add_artist(ob)

    if cp in varIdxs:
        c_model = varIdxs.index(cp)

        # Define KDE limits.
        x_rang = .1 * (xp_max - xp_min)
        x_kde = np.mgrid[xp_min - x_rang:xp_max + x_rang:100j]
        kernel_cl = stats.gaussian_kde(model_done[c_model])
        # KDE for plotting.
        kde = np.reshape(kernel_cl(x_kde).T, x_kde.shape)
        plt.plot(x_kde, kde / max(kde), color='k', lw=1.5)

        # Obtain the bin values and edges using numpy
        hist, bin_edges = np.histogram(model_done[c_model], bins=25)
        # Plot bars with the proper positioning, height, and width.
        plt.bar(
            (bin_edges[1:] + bin_edges[:-1]) * .5, hist / float(hist.max()),
            width=(bin_edges[1] - bin_edges[0]), color='grey', alpha=0.3)

        # MAP
        plt.axvline(
            x=xp, linestyle='--', color='red', zorder=4,
            label=("MAP (" + p + ")").format(xp))
        # Mode (using KDE)
        x_mode = x_kde[np.argmax(kde)]
        plt.axvline(
            x=x_mode, linestyle='--', color='cyan', zorder=4,
            label=("Mode (" + p + ")").format(x_mode))
        # Mean
        x_mean = np.mean(model_done[c_model])
        plt.axvline(
            x=x_mean, linestyle='--', color='blue', zorder=4,
            label=("Mean (" + p + ")").format(x_mean))
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


def return_intersection(hist_1, hist_2):
    """
    The ranges of both histograms must coincide for this function to work.

    Source: https://mpatacchiola.github.io/blog/2016/11/12/
            the-simplest-classifier-histogram-intersection.html
    """
    minima = np.minimum(hist_1, hist_2)
    intersection = np.true_divide(np.sum(minima), np.sum(hist_2))
    return intersection


def pl_pdf_half(par_name, gs, min_max_p, varIdxs, model_done):
    '''
    Parameter posterior plot for 1st and 2nd halfs.
    '''
    plot_dict = {
        'metal': [2, 4, 0, 2, 0], 'age': [4, 6, 0, 2, 1],
        'ext': [6, 8, 0, 2, 2], 'dist': [4, 6, 2, 4, 3],
        'mass': [6, 8, 2, 4, 4], 'binar': [6, 8, 4, 6, 5]
    }
    labels = [r'$z$', r'$\log(age)$', r'$E_{{(B-V)}}$', r'$(m-M)_o$',
              r'$M\,(M_{{\odot}})$', r'$b_{{frac}}$']

    gs_x1, gs_x2, gs_y1, gs_y2, cp = plot_dict[par_name]
    ax = plt.subplot(gs[gs_y1:gs_y2, gs_x1:gs_x2])

    ld_p = labels[cp]

    # Set x axis limit.
    xp_min, xp_max = min_max_p[cp]
    ax.locator_params(nbins=5)
    # Set minor ticks
    ax.minorticks_on()
    ax.tick_params(axis='y', which='major', labelleft=False, labelbottom=False)

    if cp in varIdxs:
        c_model = varIdxs.index(cp)
        h_min, h_max = min(model_done[c_model]), max(model_done[c_model])
        half = int(.5 * len(model_done[c_model]))
        # 1st half
        hist_1 = plt.hist(
            model_done[c_model][:half], bins=20, edgecolor='k',
            facecolor='b', linewidth=.5, alpha=0.3, range=[h_min, h_max],
            label="1st half")[0]
        # 2nd half
        hist_2 = plt.hist(
            model_done[c_model][half:], bins=20, edgecolor='k',
            facecolor='g', linewidth=.5, alpha=0.3, range=[h_min, h_max],
            label="2nd half")[0]

        # Add textbox.
        itsct = return_intersection(hist_1, hist_2)
        ob = offsetbox.AnchoredText(
            ld_p + "\n" + r"$H_{{1}} \cap H_{{2}}={:.2f}$".format(itsct),
            pad=0.1, loc=2, prop=dict(size=10))
        ob.patch.set(alpha=0.8)
        ax.add_artist(ob)
        # Legend.
        handles, labels = ax.get_legend_handles_labels()
        leg = ax.legend(handles, labels, loc='upper right', numpoints=1,
                        fontsize=10)
        leg.get_frame().set_alpha(0.7)
        cur_ylim = ax.get_ylim()
        ax.set_ylim([0, cur_ylim[1]])


def pl_param_chain(
        par_name, gs, min_max_p, cp_r, nwalkers, nsteps, nburn, m_accpt_fr,
        varIdxs, pars_chains):
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
        plt.title("Chains (walkers): {} ;  MAF: {:.3f}".format(
            nwalkers, m_accpt_fr))
    if cp == 5:
        plt.xlabel("Steps")
    else:
        ax.tick_params(labelbottom=False)
    plt.ylabel(labels[cp])
    ax.set_xlim(0, nburn + nsteps)

    if cp in varIdxs:
        c_model = varIdxs.index(cp)
        pre_bi, post_bi = pars_chains
        color = iter(cm.rainbow(np.linspace(0, 1, nwalkers)))
        for w1, w2 in zip(*[pre_bi[c_model].T, post_bi[c_model].T]):
            # Burn-in stage
            plt.plot(range(nburn), w1, c='grey', lw=.5, alpha=0.5)
            # Post burn-in.
            c = next(color)
            plt.plot(np.arange(nburn, nburn + nsteps), w2, c=c, lw=.5,
                     alpha=0.5)
        plt.axhline(y=float(cp_r[cp]), color='k', lw=1.5, zorder=4)
        ymin, ymax = np.min(post_bi[c_model]), np.max(post_bi[c_model])
        ax.set_ylim(ymin, ymax)


def plot(N, *args):
    '''
    Handle each plot separately.
    '''
    plt.style.use('seaborn-darkgrid')
    plt_map = {
        0: [pl_2_param_dens, args[0] + ' density map'],
        1: [pl_param_pf, args[0] + ' probability function'],
        2: [pl_pdf_half, args[0] + ' 1st and 2nd halfs of pdf'],
        3: [pl_param_chain, args[0] + ' sampler chain']
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

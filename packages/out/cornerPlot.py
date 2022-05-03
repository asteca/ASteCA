
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from scipy.ndimage import gaussian_filter
from matplotlib.colors import LinearSegmentedColormap, colorConverter
from .prep_plots import SigmaEllipse


def hist2d(
    ax, x, y, bins=20, range=None, weights=None, levels=None, smooth=None,
        color=None, quiet=True, plot_datapoints=True, plot_density=True,
        plot_contours=True, no_fill_contours=False, fill_contours=False,
        contour_kwargs=None, contourf_kwargs=None, data_kwargs=None, **kwargs):
    """
    Plot a 2-D histogram of samples.

    Source: https://corner.readthedocs.io

    Copyright (c) 2013-2016 Daniel Foreman-Mackey

    Parameters
    ----------
    x : array_like[nsamples,]
       The samples.
    y : array_like[nsamples,]
       The samples.
    quiet : bool
        If true, suppress warnings for small datasets.
    levels : array_like
        The contour levels to draw.
    ax : matplotlib.Axes
        A axes instance on which to add the 2-D histogram.
    plot_datapoints : bool
        Draw the individual data points.
    plot_density : bool
        Draw the density colormap.
    plot_contours : bool
        Draw the contours.
    no_fill_contours : bool
        Add no filling at all to the contours (unlike setting
        ``fill_contours=False``, which still adds a white fill at the densest
        points).
    fill_contours : bool
        Fill the contours.
    contour_kwargs : dict
        Any additional keyword arguments to pass to the `contour` method.
    contourf_kwargs : dict
        Any additional keyword arguments to pass to the `contourf` method.
    data_kwargs : dict
        Any additional keyword arguments to pass to the `plot` method when
        adding the individual data points.
    """
    # Set the default range based on the data range if not provided.
    if range is None:
        range = [[x.min(), x.max()], [y.min(), y.max()]]

    # Set up the default plotting arguments.
    if color is None:
        color = "k"

    # Choose the default "sigma" contour levels.
    if levels is None:
        levels = 1.0 - np.exp(-0.5 * np.arange(0.5, 2.1, 0.5) ** 2)

    # This is the color map for the density plot, over-plotted to indicate the
    # density of the points near the center.
    density_cmap = LinearSegmentedColormap.from_list(
        "density_cmap", [color, (1, 1, 1, 0)])

    # This color map is used to hide the points at the high density areas.
    white_cmap = LinearSegmentedColormap.from_list(
        "white_cmap", [(1, 1, 1), (1, 1, 1)], N=2)

    # This "color map" is the list of colors for the contour levels if the
    # contours are filled.
    rgba_color = colorConverter.to_rgba(color)
    contour_cmap = [list(rgba_color) for l in levels] + [rgba_color]
    for i, l in enumerate(levels):
        contour_cmap[i][-1] *= float(i) / (len(levels) + 1)

    # We'll make the 2D histogram to directly estimate the density.
    try:
        H, X, Y = np.histogram2d(x.flatten(), y.flatten(), bins=bins,
                                 range=list(map(np.sort, range)),
                                 weights=weights)
    except ValueError:
        print("  It looks like at least one of your sample columns "
              "  have no dynamic range. You could try using the "
              "  'range' argument.")
        return

    if smooth is not None:
        H = gaussian_filter(H, smooth)

    if plot_contours or plot_density:
        # Compute the density levels.
        Hflat = H.flatten()
        inds = np.argsort(Hflat)[::-1]
        Hflat = Hflat[inds]
        sm = np.cumsum(Hflat)
        sm /= sm[-1]
        V = np.empty(len(levels))
        for i, v0 in enumerate(levels):
            try:
                V[i] = Hflat[sm <= v0][-1]
            except:
                V[i] = Hflat[0]
        V.sort()
        m = np.diff(V) == 0
        if np.any(m) and not quiet:
            print("  Too few points to create valid contours")
        while np.any(m):
            V[np.where(m)[0][0]] *= 1.0 - 1e-4
            m = np.diff(V) == 0
        V.sort()

        # Compute the bin centers.
        X1, Y1 = 0.5 * (X[1:] + X[:-1]), 0.5 * (Y[1:] + Y[:-1])

        # Extend the array for the sake of the contours at the plot edges.
        H2 = H.min() + np.zeros((H.shape[0] + 4, H.shape[1] + 4))
        H2[2:-2, 2:-2] = H
        H2[2:-2, 1] = H[:, 0]
        H2[2:-2, -2] = H[:, -1]
        H2[1, 2:-2] = H[0]
        H2[-2, 2:-2] = H[-1]
        H2[1, 1] = H[0, 0]
        H2[1, -2] = H[0, -1]
        H2[-2, 1] = H[-1, 0]
        H2[-2, -2] = H[-1, -1]
        X2 = np.concatenate([
            X1[0] + np.array([-2, -1]) * np.diff(X1[:2]),
            X1,
            X1[-1] + np.array([1, 2]) * np.diff(X1[-2:]),
        ])
        Y2 = np.concatenate([
            Y1[0] + np.array([-2, -1]) * np.diff(Y1[:2]),
            Y1,
            Y1[-1] + np.array([1, 2]) * np.diff(Y1[-2:]),
        ])

    if plot_datapoints:
        if data_kwargs is None:
            data_kwargs = dict()
        data_kwargs["color"] = data_kwargs.get("color", color)
        data_kwargs["ms"] = data_kwargs.get("ms", 2.0)
        data_kwargs["mec"] = data_kwargs.get("mec", "none")
        data_kwargs["alpha"] = data_kwargs.get("alpha", 0.1)
        ax.plot(x, y, "o", zorder=-1, rasterized=True, **data_kwargs)

    # Plot the base fill to hide the densest data points.
    if (plot_contours or plot_density) and not no_fill_contours:
        ax.contourf(X2, Y2, H2.T, [V.min(), H.max()],
                    cmap=white_cmap, antialiased=False)

    if plot_contours and fill_contours:
        if contourf_kwargs is None:
            contourf_kwargs = dict()
        contourf_kwargs["colors"] = contourf_kwargs.get("colors", contour_cmap)
        contourf_kwargs["antialiased"] = contourf_kwargs.get("antialiased",
                                                             False)
        ax.contourf(
            X2, Y2, H2.T, np.concatenate([[0], V, [H.max() * (1 + 1e-4)]]),
            **contourf_kwargs)

    # Plot the density map. This can't be plotted at the same time as the
    # contour fills.
    elif plot_density:
        ax.pcolor(X, Y, H.max() - H.T, cmap=density_cmap)

    # Plot the contour edge colors.
    if plot_contours:
        if contour_kwargs is None:
            contour_kwargs = dict()
        contour_kwargs["colors"] = contour_kwargs.get("colors", color)
        ax.contour(X2, Y2, H2.T, V, **contour_kwargs)

    # ax.set_xlim(range[0])
    # ax.set_ylim(range[1])


def pl_2_param_dens(_2_params, gs, labels, min_max_p2, varIdxs, params_trace):
    '''
    Parameter vs parameters density map.
    '''
    plot_dict = {
        'metal-age': [0, 2, 2, 4, 0, 1],
        'metal-ext': [0, 2, 4, 6, 0, 2],
        'metal-dr': [0, 2, 6, 8, 0, 3],
        'metal-dist': [0, 2, 8, 10, 0, 4],
        'metal-beta': [0, 2, 10, 12, 0, 5],
        'age-ext': [2, 4, 4, 6, 1, 2],
        'age-dr': [2, 4, 6, 8, 1, 3],
        'age-dist': [2, 4, 8, 10, 1, 4],
        'age-beta': [2, 4, 10, 12, 1, 5],
        'ext-dr': [4, 6, 6, 8, 2, 3],
        'ext-dist': [4, 6, 8, 10, 2, 4],
        'ext-beta': [4, 6, 10, 12, 2, 5],
        'dr-dist': [6, 8, 8, 10, 3, 4],
        'dr-beta': [6, 8, 10, 12, 3, 5],
        'dist-beta': [8, 10, 10, 12, 4, 5]
    }

    gs_x1, gs_x2, gs_y1, gs_y2, mx, my = plot_dict[_2_params]
    x_label, y_label = labels[mx], labels[my]

    ax = plt.subplot(gs[gs_y1:gs_y2, gs_x1:gs_x2])

    # To specify the number of ticks on both or any single axes
    ax.locator_params(nbins=5)
    if gs_x1 == 0:
        plt.ylabel(y_label, fontsize=11)
        plt.yticks(rotation=45)
    else:
        ax.tick_params(labelleft=False)
    if gs_y2 == 12:
        plt.xlabel(x_label, fontsize=11)
        plt.xticks(rotation=45)
    else:
        ax.tick_params(labelbottom=False)
    plt.minorticks_on()

    if mx in varIdxs and my in varIdxs:
        mx_model, my_model = varIdxs.index(mx), varIdxs.index(my)

        ax.set_title(r"$\rho={:.2f}$".format(np.corrcoef(
            [params_trace[mx_model], params_trace[my_model]])[0][1]),
            fontsize=11)

        hist2d(ax, params_trace[mx_model], params_trace[my_model])

        mean_pos, width, height, theta = SigmaEllipse(np.array([
            params_trace[mx_model], params_trace[my_model]]).T)
        # Plot 95% confidence ellipse.
        plt.scatter(
            mean_pos[0], mean_pos[1], marker='x', c='b', s=30, linewidth=2,
            zorder=4)
        ellipse = Ellipse(xy=mean_pos, width=width, height=height, angle=theta,
                          edgecolor='r', fc='None', lw=.7, zorder=4)
        ax.add_patch(ellipse)

    xp_min, xp_max, yp_min, yp_max = min_max_p2
    ax.set_xlim([xp_min, xp_max])
    ax.set_ylim([yp_min, yp_max])

    # Grid won't respect 'zorder':
    # https://github.com/matplotlib/matplotlib/issues/5045
    # So we plot the grid behind everything else manually.
    xlocs, xlabels = plt.xticks()
    ylocs, ylabels = plt.yticks()
    for xt in xlocs:
        plt.axvline(x=xt, linestyle='-', color='w', zorder=-4)
    for yt in ylocs:
        plt.axhline(y=yt, linestyle='-', color='w', zorder=-4)


def pl_param_pf(
    par_name, gs, labels, min_max_p, varIdxs, mean_sol, map_sol, median_sol,
        mode_sol, pardist_kde, model_done):
    '''
    Parameter posterior plot.
    '''
    plot_dict = {
        'metal': [0, 2, 0, 2, 0], 'age': [2, 4, 2, 4, 1],
        'ext': [4, 6, 4, 6, 2], 'dr': [6, 8, 6, 8, 3],
        'dist': [8, 10, 8, 10, 4], 'beta': [10, 12, 10, 12, 5]
    }

    gs_x1, gs_x2, gs_y1, gs_y2, cp = plot_dict[par_name]

    frm = ["{:.4f}", "{:.3f}", "{:.3f}", "{:.3f}", "{:.0f}", "{:.2f}"]

    ld_p = labels[cp]
    p = frm[cp]

    ax = plt.subplot(gs[gs_y1:gs_y2, gs_x1:gs_x2])
    plt.title(ld_p, fontsize=11)

    # Set x axis limit.
    xp_min, xp_max = min_max_p[cp]
    ax.set_xlim(xp_min, xp_max)
    ax.locator_params(nbins=5)
    # Set minor ticks
    ax.minorticks_on()
    if cp == 5:
        plt.xlabel(ld_p, fontsize=11)
        plt.xticks(rotation=45)
    # else:
    #     ax.tick_params(labelbottom=False)
    ax.tick_params(axis='y', which='major', labelleft=False)

    if cp in varIdxs:
        c_model = varIdxs.index(cp)

        # Plot KDE.
        if pardist_kde[c_model]:
            x_kde, par_kde = pardist_kde[c_model]
            plt.plot(x_kde, par_kde / max(par_kde), color='k', lw=1.5)

        # Obtain the bin values and edges using numpy
        hist, bin_edges = np.histogram(model_done[c_model], bins='auto')
        if len(bin_edges) > 25:
            hist, bin_edges = np.histogram(model_done[c_model], bins=20)
        # Plot bars with the proper positioning, height, and width.
        plt.bar(
            (bin_edges[1:] + bin_edges[:-1]) * .5, hist / float(hist.max()),
            width=(bin_edges[1] - bin_edges[0]), color='grey', alpha=0.3)

        # Mean
        plt.axvline(
            x=mean_sol[cp], linestyle='--', color='blue', zorder=4,
            label=("Mean (" + p + ")").format(mean_sol[cp]))
        # MAP
        plt.axvline(
            x=map_sol[cp], linestyle='--', color='red', zorder=4,
            label=("MAP (" + p + ")").format(map_sol[cp]))
        # Median
        plt.axvline(
            x=median_sol[cp], linestyle=':', color='green', zorder=4,
            label=("Median (" + p + ")").format(median_sol[cp]))
        # Mode
        plt.axvline(
            x=mode_sol[cp], linestyle='--', color='cyan', zorder=4,
            label=("Mode (" + p + ")").format(mode_sol[cp]))

        # 16th and 84th percentiles (1 sigma) around median.
        ph = np.percentile(model_done[c_model], 84)
        pl = np.percentile(model_done[c_model], 16)
        plt.axvline(
            x=ph, linestyle=':', color='orange', lw=1.5, zorder=4,
            label=("16-84th perc\n(" + p + ", " + p + ")").format(pl, ph))
        plt.axvline(x=pl, linestyle=':', color='orange', lw=1.5, zorder=4)

        cur_ylim = ax.get_ylim()
        ax.set_ylim([0, cur_ylim[1]])
        plt.legend(fontsize='small')


def plot(N, *args):
    """
    Handle each plot separately.
    """
    plt_map = {
        0: [pl_2_param_dens, args[0] + ' density map'],
        1: [pl_param_pf, args[0] + ' probability function'],
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

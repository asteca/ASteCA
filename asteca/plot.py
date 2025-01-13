import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes

from .cluster import Cluster
from .synthetic import Synthetic


def radec(cluster: Cluster, ax: Axes) -> Axes:
    """Generate a (RA, DEC) plot.

    :param cluster: :py:class:`Cluster <asteca.cluster.Cluster>` object with the
        loaded data for the observed cluster
    :type cluster: Cluster
    :param ax: Matplotlib axis where to draw the plot
    :type ax: Axes

    :return: Matplotlib axis object
    :rtype: Axes
    """
    ra = cluster.ra_v
    dec = cluster.dec_v
    mag = cluster.mag_p

    msk = ~np.isnan(mag)
    ra = ra[msk]
    dec = dec[msk]
    mag = mag[msk]

    # Mag to flux
    sizes = 10 ** (mag / -2.5)
    # Normalize to 1
    sizes /= sizes.max()
    # Set min, max
    sizes = 1 + 75 * sizes

    plt.scatter(ra, dec, s=sizes, c="k", alpha=0.7)
    plt.xlabel("RA")
    plt.ylabel("DEC")
    plt.gca().invert_xaxis()

    return ax


def cluster(
    cluster: Cluster,
    ax: Axes,
    color_idx: int = 0,
    binar_probs: np.ndarray | None = None,
    prob_binar_cut: float = 0.5,
) -> Axes:
    """Generate a color-magnitude plot.

    :param cluster: :py:class:`Cluster <asteca.cluster.Cluster>` object with the
        loaded data for the observed cluster
    :type cluster: Cluster
    :param ax: Matplotlib axis where to draw the plot
    :type ax: Axes
    :param color_idx: Index of the color to plot. If ``0`` (default), plot the
        first color. If ``1`` plot the second color. Defaults to ``0``
    :type color_idx: int
    :param binar_probs: Array with probabilities of being a binary system for each
        observed star; defaults to ``None``
    :type binar_probs: np.ndarray | None
    :param prob_binar_cut: Probabilities value that separates single systems from
        binary systems; defaults to ``0.5``
    :type prob_binar_cut: float

    :raises ValueError: If ``color_idx`` is not ``0`` or ``1``

    :return: Matplotlib axis object
    :rtype: Axes
    """
    if color_idx > 1:
        raise ValueError(
            f"Wrong 'color_idx' value ({color_idx}), should be one of: [0, 1]"
        )

    if binar_probs is None:
        ax.scatter(
            cluster.colors_p[color_idx],
            cluster.mag_p,
            c="green",
            alpha=0.5,
            label=f"Observed, N={len(cluster.mag_p)}",
        )
    else:
        msk_binar = binar_probs > prob_binar_cut
        ax.scatter(
            cluster.colors_p[color_idx][~msk_binar],
            cluster.mag_p[~msk_binar],
            c="grey",
            # c=binar_probs[~msk_binar],
            marker="o",
            alpha=0.5,
            label=f"Observed (single), N={len(cluster.mag_p[~msk_binar])}",
        )
        ax.scatter(
            cluster.colors_p[color_idx][msk_binar],
            cluster.mag_p[msk_binar],
            # c="red",
            c=binar_probs[msk_binar],
            marker="s",
            alpha=0.5,
            label=f"Observed (binary), N={len(cluster.mag_p[msk_binar])}",
        )

    ax.set_ylim(max(cluster.mag_p) + 0.5, min(cluster.mag_p) - 1)

    mag_col = cluster.magnitude
    col_col = cluster.color
    if color_idx == 1:
        col_col = cluster.color2
    if col_col is not None:
        ax.set_xlabel(col_col)
    if mag_col is not None:
        ax.set_ylabel(mag_col)
    ax.legend()

    return ax


def synthetic(
    synth: Synthetic,
    ax: Axes,
    fit_params: dict,
    isoch_arr: np.ndarray | None = None,
    color_idx: int = 0,
) -> Axes:
    """Generate a color-magnitude plot for a synthetic cluster.

    The synthetic cluster is generated using the fundamental parameter values
    given in the ``fit_params`` dictionary.

    :param synth: :py:class:`Synthetic <asteca.synthetic.Synthetic>` object with the
        data required to generate synthetic clusters
    :type synth: Synthetic
    :param ax: Matplotlib axis where to draw the plot
    :type ax: Axes
    :param fit_params: Dictionary with the values for the fundamental parameters
        that were **not** included in the ``fix_params`` dictionary when the
        :py:class:`Synthetic` object was calibrated
        (:py:meth:`calibrate` method).
    :type fit_params: dict
    :param isoch_arr: Array with the isochrone data to plot; defaults to ``None``
    :type isoch_arr: np.ndarray | None
    :param color_idx: Index of the color to plot. If ``0`` (default), plot the
        first color. If ``1`` plot the second color. Defaults to ``0``
    :type color_idx: int

    :raises ValueError: If ``color_idx`` is not ``0`` or ``1``

    :return: Matplotlib axis object
    :rtype: Axes
    """
    if color_idx > 1:
        raise ValueError(
            f"Wrong 'color_idx' value ({color_idx}), should be one of: [0, 1]"
        )

    # Generate synthetic cluster.
    synth_clust = synth.generate(fit_params, full_arr_flag=True)
    if synth.binar_flag is True:
        binar_idx = ~np.isnan(synth_clust[-1])
    else:
        # No binary systems
        binar_idx = np.full(synth_clust.shape[1], False)

    y_synth = synth_clust[0]
    x_synth = synth_clust[1]
    if color_idx == 1:
        x_synth = synth_clust[2]
    # Single synthetic systems
    ax.scatter(
        x_synth[~binar_idx],
        y_synth[~binar_idx],
        marker="^",
        c="#519ddb",
        alpha=0.5,
        label=f"Synthetic (single), N={len(x_synth[~binar_idx])}",
    )
    # Binary synthetic systems
    ax.scatter(
        x_synth[binar_idx],
        y_synth[binar_idx],
        marker="v",
        c="#F34C4C",
        alpha=0.5,
        label=f"Synthetic (binary), N={len(x_synth[binar_idx])}",
    )

    plt.ylabel(synth.isochs.magnitude)
    c1, c2 = synth.isochs.color
    if color_idx == 1:
        if synth.isochs.color2 is None:
            raise ValueError("No second color available")
        c1, c2 = synth.isochs.color2
    plt.xlabel(f"{c1}-{c2}")
    ax.set_ylim(max(y_synth) + 0.5, min(y_synth) - 1)
    ax.legend()

    if isoch_arr is not None:
        c_idx = 1
        if color_idx == 1:
            c_idx = 2
        ax.plot(isoch_arr[c_idx], isoch_arr[0], c="k")

    return ax


def get_isochrone(
    synth: Synthetic,
    fit_params: dict,
    color_idx: int = 0,
) -> np.ndarray:
    """Generate an isochrone for plotting.

    The isochrone is generated using the fundamental parameter values
    given in the ``fit_params`` dictionary.

    :param synth: :py:class:`Synthetic <asteca.synthetic.Synthetic>` object with the
        data required to generate synthetic clusters
    :type synth: Synthetic
    :param fit_params: Dictionary with the values for the fundamental parameters
        that were **not** included in the ``fix_params`` dictionary when the
        :py:class:`Synthetic` object was calibrated
        (:py:meth:`calibrate` method).
    :type fit_params: dict
    :param color_idx: Index of the color to plot. If ``0`` (default), plot the
        first color. If ``1`` plot the second color. Defaults to ``0``
    :type color_idx: int

    :return: Array with the isochrone data to plot
    :rtype: np.ndarray
    """
    # Generate displaced isochrone
    fit_params_copy = dict(fit_params)

    # Generate physical synthetic cluster to extract the mas mass
    isochrone_full = synth.generate(fit_params_copy, full_arr_flag=True)
    # Extract max mass
    max_mass = isochrone_full[synth.m_ini_idx].max()

    # Generate displaced isochrone
    fit_params_copy["DR"] = 0.0
    isochrone = synth.generate(fit_params_copy, plot_flag=True)

    # Apply max mass filter to isochrone
    msk = isochrone[synth.m_ini_idx] < max_mass

    # Generate proper array for plotting
    isochrone = np.array([isochrone[0], isochrone[color_idx + 1]])[:, msk]

    return isochrone

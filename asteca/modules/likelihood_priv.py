from typing import Literal

import numpy as np
from astropy.stats import calculate_bin_edges
from fast_histogram import histogram2d
from scipy.special import gammaln

from . import HAS_NUMBA as _HAS_NUMBA
from . import njit

try:
    # C core of fast_histogram, skips the (costly) argument validation performed by
    # the public histogram2d() wrapper
    from fast_histogram._histogram_core import _histogram2d
except ImportError:  # private module, fall back to the public function

    def _histogram2d(x, y, nx, xmin, xmax, ny, ymin, ymax):
        return histogram2d(x, y, bins=[nx, ny], range=[[xmin, xmax], [ymin, ymax]])


def bin_edges_f(
    bin_method: Literal["blocks", "knuth", "scott", "freedman", "fixed"],
    mag: np.ndarray,
    colors: list[np.ndarray],
) -> tuple[list, list]:
    """Calculate bin edges for the Hess diagram.

    This function calculates the bin edges for the Hess diagram, using
    different methods.

    :param bin_method: Method used for binning the data.
    :type bin_method: Literal["blocks", "knuth", "scott", "freedman", "fixed"]
    :param mag: Array of magnitudes.
    :type mag: np.ndarray
    :param colors: List of arrays of colors.
    :type colors: list[np.ndarray]

    :return: Bin ranges and number of bins.
    :rtype: tuple[list, list]
    """

    bin_edges = []

    if bin_method == "fixed":
        N_mag, N_col = 15, 10
        # Magnitude
        mag_min, mag_max = np.nanmin(mag), np.nanmax(mag)
        bin_edges.append(np.linspace(mag_min, mag_max, N_mag))
        # Colors
        for col in colors:
            col_min, col_max = np.nanmin(col), np.nanmax(col)
            bin_edges.append(np.linspace(col_min, col_max, N_col))
    else:
        bin_edges.append(calculate_bin_edges(mag[~np.isnan(mag)], bins=bin_method))
        for col in colors:
            bin_edges.append(calculate_bin_edges(col[~np.isnan(col)], bins=bin_method))

    # Extract ranges and number of bins for each dimension (magnitude and colors),
    # used by histogram2d
    ranges, Nbins = [], []
    for be in bin_edges:
        # Add a small buffer to the edges to rejecting stars on the edges when
        # generating the Hess diagram
        ranges.append([be[0] - 0.001, be[-1] + 0.001])
        Nbins.append(len(be))

    return ranges, Nbins


def lkl_data(
    mag_v: np.ndarray, colors_v: list[np.ndarray], ranges: list, Nbins: list
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Prepare data for likelihood calculation.

    This function calculates the Hess diagram of the observed cluster and
    prepares the data for the likelihood calculation.

    :param mag_v: Array of magnitudes.
    :type mag_v: np.ndarray
    :param colors_v: List of arrays of colors.
    :type colors_v: list[np.ndarray]
    :param ranges: Per-dimension ranges.
    :type ranges: list
    :param Nbins: Per-dimension total number of bins.
    :type Nbins: list

    :return: Indexes of bins with stars, flattened histogram, and pre-computed table
        for likelihood calculation.
    :rtype: tuple[np.ndarray, np.ndarray, np.ndarray]
    """

    # Obtain histogram for observed cluster.
    hess_diag = []
    for i, col in enumerate(colors_v):
        # Fast 2D histogram
        hess_diag.append(
            histogram2d(
                mag_v,
                col,
                range=[
                    [ranges[0][0], ranges[0][1]],
                    [ranges[i + 1][0], ranges[i + 1][1]],
                ],
                bins=[Nbins[0], Nbins[i + 1]],
            )
        )

    # Flatten array
    cl_histo_f = []
    for i, diag in enumerate(hess_diag):
        cl_histo_f += list(np.array(diag).ravel())
    cl_histo_f = np.array(cl_histo_f)

    # Index of bins where stars were observed
    cl_z_idx = cl_histo_f != 0

    # Remove all bins where n_i=0 (no observed stars)
    cl_histo_f_z = cl_histo_f[cl_z_idx].astype(int)

    # table2d[n, m]: pre-computed table for likelihood calculation, where n is the
    # number of observed stars in a bin and m is the number of synthetic stars in that
    # bin. This table is used to speed up the likelihood calculation, as it avoids the
    # need to compute the gammaln function for each bin and each synthetic cluster.
    # The table is computed for all possible values of n and m up to a certain maximum,
    # which is determined based on the observed histogram.

    # The maximum number of synthetic stars in a bin should be large enough to cover
    # the range of possible values given the number of observed stars, but not too
    # large to avoid unnecessary computation. The following heuristic is used to
    # determine the maximum number of synthetic stars:
    N_obs = len(mag_v)
    if N_obs < 50:
        N_st_max = N_obs * 5
    elif N_obs < 500:
        N_st_max = N_obs * 2
    elif N_obs < 1000:
        N_st_max = N_obs * 1
    else:
        N_st_max = N_obs * 0.5
    # Create a 2D table for all combinations of n and m, where n is the number of
    # observed stars in a bin and m is the number of synthetic stars in that bin.
    # The table is computed using the gammaln function, which is used in the
    # likelihood calculation.
    n_idx, m_idx = np.arange(cl_histo_f_z.max() + 1), np.arange(N_st_max)
    N2, M2 = np.meshgrid(n_idx, m_idx, indexing="ij")
    table2d = gammaln(N2 + M2 + 0.5) - gammaln(M2 + 0.5)

    # These variables are used by the likelihood functions called by the get() method
    return cl_z_idx, cl_histo_f_z, table2d


def _median(x: np.ndarray) -> float:
    """Median of a 1D array, equal to ``np.median`` but faster as it uses
    ``np.partition`` directly"""
    N = x.shape[0]
    k = N // 2
    if N % 2:
        return np.partition(x, k)[k]
    xp = np.partition(x, (k - 1, k))
    return 0.5 * (xp[k - 1] + xp[k])


def tremmel(
    hist_args: np.ndarray,
    nz_pos: np.ndarray,
    row_off: np.ndarray,
    table_flat: np.ndarray,
    m_max: int,
    obs_mag_median: float,
    obs_col_median: float,
    cmd_widths: tuple[float, float],
    max_lkl: None | float,
    min_lkl: float,
    synth_clust: np.ndarray,
) -> float:
    r"""Poisson likelihood ratio as defined in Tremmel et al (2013), Eq 10 with
    v_{i,j}=1 (see documentation). This returns the log likelihood.

    :param hist_args: Per-color ``(nx, xmin, xmax, ny, ymin, ymax)`` arguments for
        the 2D histogram, shape ``(N_colors, 6)``
    :type hist_args: np.ndarray
    :param nz_pos: Positions in the flattened Hess diagram of the bins where the
        number of observed stars is not 0
    :type nz_pos: np.ndarray
    :param row_off: Offset in ``table_flat`` of the row associated to the number of
        observed stars in each non-empty bin
    :type row_off: np.ndarray
    :param table_flat: Flattened pre-computed table for likelihood calculation
    :type table_flat: np.ndarray
    :param m_max: Largest number of synthetic stars per bin covered by the table
        (number of columns minus 1). Larger counts are clipped to this value
    :type m_max: int
    :param obs_mag_median: Median observed magnitude
    :type obs_mag_median: float
    :param obs_col_median: Median observed color
    :type obs_col_median: float
    :param cmd_widths: Widths of the Hess diagram ranges for the magnitude and
        (first) color, used to scale the no-overlap penalty
    :type cmd_widths: tuple[float, float]
    :param max_lkl: Maximum likelihood value, used for normalization
    :type max_lkl: None | float
    :param min_lkl: Minimum likelihood value, obtained when no synthetic star falls
        in any observed bin. Used for normalization
    :type min_lkl: float
    :param synth_clust: Synthetic cluster data.
    :type synth_clust: np.ndarray

    :return: Log likelihood value.
    :rtype: float
    """
    # If synthetic cluster is empty, assign '10' indicating a large distance
    if synth_clust.size == 0:
        return 10

    mag = synth_clust[0]
    # The 1/2 term arises from using a Jeffreys prior (Eq 8 in Tremmel et al). Using a
    # flat prior instead, that 1/2 term becomes a 1
    # Equivalent to table2d[cl_histo_f_z, syn_histo_f_z].sum()
    if _HAS_NUMBA:
        tremmel_lkl = _tremmel_lkl_nb(
            mag, synth_clust[1:], hist_args, nz_pos, row_off, table_flat, m_max
        )
    else:
        # Obtain histogram for the synthetic cluster.
        syn_histo_f = np.concatenate(
            [
                _histogram2d(
                    mag,
                    synth_clust[i + 1],
                    int(h[0]),
                    h[1],
                    h[2],
                    int(h[3]),
                    h[4],
                    h[5],
                ).ravel()
                for i, h in enumerate(hist_args)
            ]
        )
        # Remove all bins where n_i = 0 (no observed stars). Counts beyond the
        # table's last column are clipped, otherwise they would read the next row
        syn_histo_f_z = np.minimum(syn_histo_f.take(nz_pos).astype(np.intp), m_max)
        tremmel_lkl = table_flat.take(row_off + syn_histo_f_z).sum()

    # If this is the call made to calibrate the likelihood
    if max_lkl is None:
        return tremmel_lkl

    # Normalize and invert the likelihood so that its values are equivalent to a
    # distance that goes from 0 (exact match) to 1 (no match). The 'no match' value
    # (min_lkl) is the likelihood obtained when no synthetic star falls in any of the
    # observed bins, i.e.: sum(table2d[cl_histo_f_z, 0]). Since table2d[n, m] grows with
    # m for n>0, tremmel_lkl >= min_lkl and thus dist <= 1 always.
    # This is only a rescaling of '1 - tremmel_lkl / max_lkl', it does not change the
    # ordering of the distances but puts the values on the same 0–1 scale for all
    # clusters (1 = no overlap)
    dist = (max_lkl - tremmel_lkl) / (max_lkl - min_lkl)

    # If m_i==0 and n_i<3, the term
    # gammaln(cl_histo_f_z + syn_histo_f_z + 0.5) - gammaln(syn_histo_f_z + 0.5)
    # will be negative, thus reducing the final likelihood value.

    # If no synthetic star falls in any observed bin (no overlap), the likelihood is
    # constant (min_lkl) no matter where the synthetic stars are, i.e.: dist=1 for all
    # these models. When most of the models sampled from the priors are in this
    # situation (e.g.: small clusters in a small region of the CMD) this tie stalls
    # the optimizer (pyABC's epsilon stays at 1). To break the tie, a penalty based on
    # the photometric distance between the synthetic and observed clusters is added.
    # This "pulls" the solution towards the region of the CMD where the observed
    # cluster is located. The penalty is only applied to the no-overlap models, so
    # that any model with some overlap (dist<1) is always better than any model with
    # no overlap (dist>=1), and the former are ranked by the likelihood alone.
    # Each dimension is scaled by the width of the Hess diagram range, so that the
    # magnitude (usually spanning a much larger range) does not dominate the color.
    # Both sums are computed identically, so they are exactly equal for no overlap
    if tremmel_lkl <= min_lkl:
        cmd_dist = ((obs_mag_median - _median(mag)) / cmd_widths[0]) ** 2 + (
            (obs_col_median - _median(synth_clust[1])) / cmd_widths[1]
        ) ** 2
        return 1 + cmd_dist

    return dist


if _HAS_NUMBA:

    @njit(cache=True)
    def _tremmel_lkl_nb(
        mag: np.ndarray,
        colors: np.ndarray,
        hist_args: np.ndarray,
        nz_pos: np.ndarray,
        row_off: np.ndarray,
        table_flat: np.ndarray,
        m_max: int,
    ) -> float:
        """Compiled (numba) version of the likelihood sum in :py:func:`tremmel`.
        A single pass bins the synthetic stars into one flattened Hess diagram
        (same binning as ``fast_histogram``: values outside ``[min, max)`` or NaN
        are dropped), then the table values of the non-empty observed bins are
        added, without the temporary arrays of the NumPy version.
        """
        N_col = hist_args.shape[0]
        nx = int(hist_args[0, 0])
        xmin, xmax = hist_args[0, 1], hist_args[0, 2]
        fx = nx / (xmax - xmin)

        # Offset of each color's Hess diagram in the flattened array
        offs = np.empty(N_col + 1, dtype=np.int64)
        offs[0] = 0
        for c in range(N_col):
            offs[c + 1] = offs[c] + nx * int(hist_args[c, 3])
        counts = np.zeros(offs[N_col], dtype=np.int64)

        for c in range(N_col):
            ny = int(hist_args[c, 3])
            ymin, ymax = hist_args[c, 4], hist_args[c, 5]
            fy = ny / (ymax - ymin)
            col = colors[c]
            for k in range(mag.shape[0]):
                x = mag[k]
                y = col[k]
                if x >= xmin and x < xmax and y >= ymin and y < ymax:
                    i = min(int((x - xmin) * fx), nx - 1)
                    j = min(int((y - ymin) * fy), ny - 1)
                    counts[offs[c] + i * ny + j] += 1

        lkl = 0.0
        for p in range(nz_pos.shape[0]):
            # Clip counts beyond the table's last column (see tremmel())
            m = min(counts[nz_pos[p]], m_max)
            lkl += table_flat[row_off[p] + m]
        return lkl


# def visual(cluster_dict, synth_clust):
#     # If synthetic cluster is empty, assign a small likelihood value.
#     if not synth_clust.any():
#         return -1.0e09

#     mag_o, colors_o = cluster_dict["mag"], cluster_dict["colors"]
#     mag_s, colors_s = synth_clust[0], synth_clust[1:]

#     N_mag, N_col = 15, 10
#     mag = list(mag_o) + list(mag_s)
#     col = list(colors_o[0]) + list(colors_s[0])
#     mag_min, mag_max = np.nanmin(mag), np.nanmax(mag)
#     bin_edges = [np.linspace(mag_min, mag_max, N_mag)]
#     col_min, col_max = np.nanmin(col), np.nanmax(col)
#     bin_edges.append(np.linspace(col_min, col_max, N_col))

#     # Obtain histogram for observed cluster.
#     cl_histo_f = []
#     for i, col_o in enumerate(colors_o):
#         hess_diag = np.histogram2d(mag_o, col_o, bins=bin_edges)[0]
#         # Flatten array
#         cl_histo_f += list(hess_diag.ravel())
#     cl_histo_f = np.array(cl_histo_f)
#     # Down sample histogram
#     # msk = cl_histo_f > 5
#     # cl_histo_f[msk] = 5

#     syn_histo_f = []
#     for i, col_s in enumerate(colors_s):
#         hess_diag = np.histogram2d(mag_s, col_s, bins=bin_edges)[0]
#         # Flatten array
#         syn_histo_f += list(hess_diag.ravel())
#     syn_histo_f = np.array(syn_histo_f)
#     # Down sample histogram
#     # msk = syn_histo_f > 5
#     # syn_histo_f[msk] = 5

#     # return -sum(abs(cluster_dict["hist_down_samp"]-syn_histo_f))

#     # create a mask for where each data set is non-zero
#     m1 = cl_histo_f != 0
#     m2 = syn_histo_f != 0
#     m1_area = m1.sum()
#     m2_area = m2.sum()
#     tot_area = m1_area + m2_area
#     # use a logical and to create a combined map where both datasets are non-zero
#     ovrlp_area = np.logical_and(m1, m2).sum()

#     # # calculate the overlapping density, where 0.5 is the bin width
#     # ol_density = np.abs((cl_histo_f - syn_histo_f) * 0.5)[ol]
#     # # calculate the total overlap percent
#     # h_overlap = ol_density.sum() * 100

#     h_overlap = ovrlp_area / tot_area

#     # ol_density = np.abs((cl_histo_f - cl_histo_f) * 0.5)[ol]
#     # h_overlap_0 = ol_density.sum() * 100
#     # print(h_overlap_0)
#     # breakpoint()

#     # # cl_histo_f = cluster_dict["hist_down_samp"]
#     # mi_cnst = np.clip(cl_histo_f, a_min=0, a_max=1)
#     # # Final chi.
#     # mig_chi = np.sum((cl_histo_f + mi_cnst - syn_histo_f)**2 / (cl_histo_f + 1.))
#     # mig_chi_0 = np.sum((cl_histo_f + mi_cnst)**2 / (cl_histo_f + 1.))
#     # mig_chi_opt = np.sum((cl_histo_f + mi_cnst - cl_histo_f)**2 / (cl_histo_f + 1.))
#     # print(mig_chi_opt, mig_chi_0, mig_chi)

#     # import matplotlib.pyplot as plt
#     # # plt.subplot(121)
#     # y_edges, x_edges = bin_edges #cluster_dict['bin_edges']
#     # for xe in x_edges:
#     #     plt.axvline(xe, c="grey", ls=":")
#     # for ye in y_edges:
#     #     plt.axhline(ye, c="grey", ls=":")
#     # plt.title(h_overlap)
#     # plt.scatter(colors_o[0], mag_o, alpha=.5)
#     # plt.scatter(colors_s[0], mag_s, alpha=.5)
#     # plt.gca().invert_yaxis()

#     # # plt.subplot(122)
#     # # plt.title(h_overlap)
#     # # plt.bar(np.arange(len(cl_histo_f)), cl_histo_f, label='obs', alpha=.5)
#     # # plt.bar(np.arange(len(syn_histo_f)), syn_histo_f, label='syn', alpha=.5)
#     # # plt.legend()
#     # plt.show()

#     return h_overlap


# def mean_dist(cluster_dict, synth_clust):
#     # If synthetic cluster is empty, assign a small likelihood value.
#     if not synth_clust.any():
#         return -1.0e09

#     # mag_o, colors_o = cluster_dict['mag'], cluster_dict['colors']
#     mag0, colors0 = cluster_dict["mag0"], cluster_dict["colors0"]
#     mag_s, colors_s = synth_clust[0], synth_clust[1:]

#     dist = (np.median(mag0) - np.median(mag_s)) ** 2 + (
#         np.median(colors0) - np.median(colors_s)
#     ) ** 2
#     return -dist

#     if len(mag_s) < 5:
#         return -1.0e09
#     import ndtest

#     P_val = ndtest.ks2d2s(mag0, colors0, mag_s, colors_s[0])

#     # import matplotlib.pyplot as plt
#     # plt.title(P_val)
#     # plt.scatter(colors_o0, mag_o, alpha=.5)
#     # plt.scatter(colors_s[0], mag_s, alpha=.5)
#     # plt.gca().invert_yaxis()
#     # plt.show()

#     return P_val


def bins_distance(
    mag_v: np.ndarray, color_v: np.ndarray, synth_clust: np.ndarray
) -> float:
    """Sum of distances to corresponding bins in the Hess diagram. Only applied
    on the first two dimensions (magnitude +  first color)

    :param mag_v: Array of magnitudes.
    :type mag_v: np.ndarray
    :param color_v: Arrays of color.
    :type color_v: np.ndarray
    :param synth_clust: Synthetic cluster data.
    :type synth_clust: np.ndarray

    :return: Sum of distances.
    :rtype: float
    """
    if not synth_clust.any():
        return 1.0e09

    # Fixed percentiles for magnitude and color
    mpercs = (0.5, 10, 20, 30, 40, 50, 60, 70, 80, 90)
    cpercs = (0.5, 10, 20, 30, 40, 50, 60, 70, 75, 80, 85, 90, 95)

    # Evaluate the magnitude and color in the defined percentiles
    perc_mag_o = np.nanpercentile(mag_v, mpercs)
    perc_colors_o = np.nanpercentile(color_v, cpercs)

    # Create a 2-dimensional array of shape: (2, len(mpercs) * len(cpercs))
    pts_o = []
    for pm in perc_mag_o:
        for pc in perc_colors_o:
            pts_o.append([pm, pc])
    pts_o = np.array(pts_o).T

    # Same for the synthetic cluster
    mag_s, colors_s = synth_clust[0], synth_clust[1]
    perc_mag_s = np.nanpercentile(mag_s, mpercs)
    perc_colors_s = np.nanpercentile(colors_s, cpercs)
    pts_s = []
    for pm in perc_mag_s:
        for pc in perc_colors_s:
            pts_s.append([pm, pc])
    pts_s = np.array(pts_s).T

    # Distance (non root squared) between the two arrays
    dist = (pts_s[0] - pts_o[0]) ** 2 + (pts_s[1] - pts_o[1]) ** 2
    # More weight to smaller magnitudes and color (top left of Hess diagram)
    weights = np.linspace(1, 0.05, len(mpercs) * len(cpercs))
    lkl = sum(dist * weights)

    # import matplotlib.pyplot as plt
    # plt.title(f"lkl={lkl:.3f}")
    # plt.scatter(colors_v[0], mag_v, alpha=0.25, c="r", label='obs')
    # plt.scatter(colors_s, mag_s, alpha=0.25, c="b", label='synth')
    # # This shows the positions where the observed Hess diagram is defined by the
    # # percentiles
    # plt.scatter(pts_o[1], pts_o[0], c="r", marker="*")
    # # Same for the synthetic Hess diagram
    # plt.scatter(pts_s[1], pts_s[0], c="b", marker="x")
    # plt.gca().invert_yaxis()
    # plt.legend()
    # plt.show()

    return lkl


def chi_square(
    ranges: list,
    Nbins: list,
    cl_z_idx: np.ndarray,
    cl_histo_f_z: np.ndarray,
    synth_clust: np.ndarray,
) -> float:
    """Calculate the chi-square value.

    This function calculates the chi-square value between the observed and
    synthetic clusters.

    :param ranges: Per-dimension ranges.
    :type ranges: list
    :param Nbins: Per-dimension total number of bins
    :type Nbins: list
    :param cl_z_idx: Index of bins where the number of stars is not 0
    :type cl_z_idx: np.ndarray
    :param cl_histo_f_z: Flattened observed Hess diagram with the empty bins removed
    :type cl_histo_f_z: np.ndarray
    :param synth_clust: Synthetic cluster data.
    :type synth_clust: np.ndarray

    :return: Chi-square value.
    :rtype: float
    """
    # If synthetic cluster is empty, assign a small likelihood value.
    if not synth_clust.any():
        return -1.0e09

    # Obtain histogram for the synthetic cluster.
    mag, colors = synth_clust[0], synth_clust[1:]
    syn_histo_f = []
    for i, col in enumerate(colors):
        hess_diag = histogram2d(
            mag,
            col,
            range=[
                [ranges[0][0], ranges[0][1]],
                [ranges[i + 1][0], ranges[i + 1][1]],
            ],
            bins=[Nbins[0], Nbins[i + 1]],
        )
        # Flatten array
        syn_histo_f += list(hess_diag.ravel())
    syn_histo_f = np.array(syn_histo_f)

    # Remove all bins where n_i = 0 (no observed stars).
    syn_histo_f_z = syn_histo_f[cl_z_idx]

    chisq = ((cl_histo_f_z - syn_histo_f_z) ** 2).sum()
    return chisq


def cmd_overlap(
    mag: np.ndarray,
    color: np.ndarray,
    synth_clust: np.ndarray,
    h_mag: float,
    h_col: float,
    chunk: int = 500,
) -> tuple[float, float]:
    """Overlap distance between the observed and synthetic clusters in the CMD
    (magnitude and first color), see :py:meth:`asteca.Likelihood.cmd_overlap`.

    :param mag: Array of observed magnitudes.
    :type mag: np.ndarray
    :param color: Array of observed colors.
    :type color: np.ndarray
    :param synth_clust: Synthetic cluster data.
    :type synth_clust: np.ndarray
    :param h_mag: Width of the Gaussian kernel in magnitude.
    :type h_mag: float
    :param h_col: Width of the Gaussian kernel in color.
    :type h_col: float
    :param chunk: Number of stars processed at once, limits the memory used.
    :type chunk: int

    :raises ValueError: If the observed cluster has fewer than two stars with valid
        photometry

    :return: Overlap distance and its approximate standard deviation for a synthetic
        cluster drawn from the same distribution as the observed one.
    :rtype: tuple[float, float]
    """
    # Scale each dimension by the width of the kernel, so that it becomes
    # k = exp(-d^2/4) with d the Euclidean distance
    obs = np.column_stack([mag / h_mag, color / h_col])
    obs = obs[~np.isnan(obs).any(axis=1)]
    N = len(obs)
    if N < 2:
        raise ValueError(
            "The observed cluster must contain at least two stars with valid photometry"
        )

    # A synthetic cluster with fewer than two stars has no overlap with the observed one
    if synth_clust.size == 0:
        return 1.0, 0.0
    syn = np.column_stack([synth_clust[0] / h_mag, synth_clust[1] / h_col])
    syn = syn[~np.isnan(syn).any(axis=1)]
    M = len(syn)
    if M < 2:
        return 1.0, 0.0

    def k_sums(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
        """Sums of k and of k^2 over all the (x, y) pairs"""
        s1 = s2 = 0.0
        for i in range(0, len(x), chunk):
            xc = x[i : i + chunk]
            d2 = (xc[:, :1] - y[:, 0]) ** 2 + (xc[:, 1:] - y[:, 1]) ** 2
            k = np.exp(-0.25 * d2)
            s1 += k.sum()
            s2 += (k * k).sum()
        return s1, s2

    if _HAS_NUMBA:
        s_oo, s2_oo = _k_sums_self_nb(obs)
        s_ss, s2_ss = _k_sums_self_nb(syn)
        s_os = _k_sums_cross_nb(obs, syn)
    else:
        s_oo, s2_oo = k_sums(obs, obs)
        s_ss, s2_ss = k_sums(syn, syn)
        s_os, _ = k_sums(obs, syn)

    # Mean kernel values over the pairs of distinct observed stars, of distinct
    # synthetic stars, and of observed-synthetic stars. Removing the pairs of a star
    # with itself (where k=1) makes these unbiased, so that a synthetic cluster
    # drawn from the same distribution as the observed one gives D=0 on average for
    # any number of stars
    k_oo = (s_oo - N) / (N * (N - 1))
    k_ss = (s_ss - M) / (M * (M - 1))
    k_os = s_os / (N * M)
    D = 1 - 2 * k_os / (k_oo + k_ss)

    # Standard deviation of D for a synthetic cluster drawn from the same
    # distribution as the observed one: leading order of the variance of the
    # unbiased MMD^2 estimator (Gretton et al. 2012), with E[k^2] in place of the
    # (smaller) second moment of the centered kernel, i.e. slightly conservative
    k2 = 0.5 * ((s2_oo - N) / (N * (N - 1)) + (s2_ss - M) / (M * (M - 1)))
    var = (2 / (N * (N - 1)) + 2 / (M * (M - 1)) + 4 / (N * M)) * k2
    D_std = np.sqrt(var) / (k_oo + k_ss)

    return float(D), float(D_std)


if _HAS_NUMBA:

    @njit(cache=True)
    def _k_sums_self_nb(x: np.ndarray) -> tuple[float, float]:
        """Compiled (numba) version of the sums of k and of k^2 over all the (x, x)
        pairs in :py:func:`cmd_overlap`, including the pairs of a star with itself
        (k=1). Since k is symmetric only the pairs i<j are evaluated, and no
        temporary arrays are created.
        """
        n = x.shape[0]
        s1 = s2 = 0.0
        for i in range(n):
            xi, yi = x[i, 0], x[i, 1]
            for j in range(i + 1, n):
                dx = xi - x[j, 0]
                dy = yi - x[j, 1]
                k = np.exp(-0.25 * (dx * dx + dy * dy))
                s1 += k
                s2 += k * k
        return n + 2 * s1, n + 2 * s2

    @njit(cache=True)
    def _k_sums_cross_nb(x: np.ndarray, y: np.ndarray) -> float:
        """Compiled (numba) version of the sum of k over all the (x, y) pairs in
        :py:func:`cmd_overlap`
        """
        s1 = 0.0
        for i in range(x.shape[0]):
            xi, yi = x[i, 0], x[i, 1]
            for j in range(y.shape[0]):
                dx = xi - y[j, 0]
                dy = yi - y[j, 1]
                s1 += np.exp(-0.25 * (dx * dx + dy * dy))
        return s1

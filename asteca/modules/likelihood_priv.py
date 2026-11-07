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


def tremmel(
    hist_args: np.ndarray,
    nz_pos: np.ndarray,
    row_off: np.ndarray,
    table_flat: np.ndarray,
    m_max: int,
    synth_clust: np.ndarray,
) -> float:
    r"""Poisson likelihood ratio as defined in Tremmel et al (2013), Eq 10 with
    v_{i,j}=1 (see documentation). This is the core shared by
    :py:meth:`asteca.Likelihood.dist` and :py:meth:`asteca.Likelihood.log_lkl`.

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
    :param synth_clust: Synthetic cluster data, ``[magnitude, color1, (color2)]``.
        Must not be empty
    :type synth_clust: np.ndarray

    :return: Log likelihood value (up to the constant terms that only depend on the
        observed cluster)
    :rtype: float
    """
    mag = synth_clust[0]
    # The 1/2 term arises from using a Jeffreys prior (Eq 8 in Tremmel et al). Using a
    # flat prior instead, that 1/2 term becomes a 1
    # Equivalent to table2d[cl_histo_f_z, syn_histo_f_z].sum()
    if _HAS_NUMBA:
        return _tremmel_lkl_nb(
            mag, synth_clust[1:], hist_args, nz_pos, row_off, table_flat, m_max
        )

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
    return table_flat.take(row_off + syn_histo_f_z).sum()


def no_overlap_dist(
    obs_mag_median: float,
    obs_col_median: float,
    cmd_widths: tuple[float, float],
    synth_clust: np.ndarray,
) -> float:
    """Distance assigned to a synthetic cluster with no stars in any of the observed
    bins.

    In that case the likelihood is constant (``min_lkl``) no matter where the
    synthetic stars are. When most of the models sampled from the priors are in this
    situation (e.g.: small clusters in a small region of the CMD) this tie stalls
    the optimizer (pyABC's epsilon stays at 1). To break the tie, a penalty based on
    the photometric distance between the synthetic and observed clusters is added to
    the 'no match' distance of ``1``. This "pulls" the solution towards the region of
    the CMD where the observed cluster is located. The penalty is only applied to the
    no-overlap models, so that any model with some overlap (distance<1) is always
    better than any model with no overlap (distance>=1), and the former are ranked by
    the likelihood alone. Each dimension is scaled by the width of the Hess diagram
    range, so that the magnitude (usually spanning a much larger range) does not
    dominate the color.

    :param obs_mag_median: Median observed magnitude
    :type obs_mag_median: float
    :param obs_col_median: Median observed color
    :type obs_col_median: float
    :param cmd_widths: Widths of the Hess diagram ranges for the magnitude and
        (first) color
    :type cmd_widths: tuple[float, float]
    :param synth_clust: Synthetic cluster data, ``[magnitude, color1, (color2)]``
    :type synth_clust: np.ndarray

    :return: Distance, ``>= 1``
    :rtype: float
    """
    cmd_dist = ((obs_mag_median - _median(synth_clust[0])) / cmd_widths[0]) ** 2 + (
        (obs_col_median - _median(synth_clust[1])) / cmd_widths[1]
    ) ** 2
    return 1 + cmd_dist


def _median(x: np.ndarray) -> float:
    """Median of a 1D array, equal to ``np.median`` but faster as it uses
    ``np.partition`` directly

    :param x: 1D array (with no ``nan`` values)
    :type x: np.ndarray

    :return: Median value
    :rtype: float
    """
    N = x.shape[0]
    k = N // 2
    if N % 2:
        return float(np.partition(x, k)[k])
    xp = np.partition(x, (k - 1, k))
    return float(0.5 * (xp[k - 1] + xp[k]))


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
        return float(lkl)


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

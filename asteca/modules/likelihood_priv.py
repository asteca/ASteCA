import numpy as np
from astropy.stats import calculate_bin_edges
from fast_histogram import histogram2d
from scipy.special import loggamma


def lkl_data(
    bin_method: str, mag_v: np.ndarray, colors_v: list[np.ndarray]
) -> tuple[list, list, np.ndarray, np.ndarray]:
    """Prepare data for likelihood calculation.

    This function calculates the Hess diagram of the observed cluster and
    prepares the data for the likelihood calculation.

    :param bin_method: Method to use for binning the data.
    :type bin_method: str
    :param mag_v: Array of magnitudes.
    :type mag_v: np.ndarray
    :param colors_v: List of arrays of colors.
    :type colors_v: list[np.ndarray]

    :return: Bin ranges, number of bins, indexes of bins with stars, and flattened histogram.
    :rtype: tuple[list, list, np.ndarray, np.ndarray]
    """
    # Obtain bin edges for each dimension, defining a grid.
    ranges, Nbins = bin_edges_f(bin_method, mag_v, colors_v)

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
    cl_histo_f_z = cl_histo_f[cl_z_idx]

    # These variables are used by the likelihood functions called by the get() method
    return ranges, Nbins, cl_z_idx, cl_histo_f_z


def bin_edges_f(
    bin_method: str, mag: np.ndarray, colors: list[np.ndarray]
) -> tuple[list, list]:
    """Calculate bin edges for the Hess diagram.

    This function calculates the bin edges for the Hess diagram, using
    different methods.

    :param bin_method: Method to use for binning the data.
    :type bin_method: str
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
        bin_edges.append(calculate_bin_edges(mag[~np.isnan(mag)], bins=bin_method))  # pyright: ignore
        for col in colors:
            bin_edges.append(
                calculate_bin_edges(col[~np.isnan(col)], bins=bin_method)  # pyright: ignore
            )

    # Extract ranges and number of bins for each dimension (magnitude and colors),
    # used by histogram2d
    ranges, Nbins = [], []
    for be in bin_edges:
        ranges.append([be[0], be[-1]])
        Nbins.append(len(be))

    return ranges, Nbins


def tremmel(
    ranges: list,
    Nbins: list,
    cl_z_idx: np.ndarray,
    cl_histo_f_z: np.ndarray,
    max_lkl: float,
    synth_clust: np.ndarray,
) -> float:
    r"""Poisson likelihood ratio as defined in Tremmel et al (2013), Eq 10 with
    v_{i,j}=1. This returns the log likelihood.

    .. math::

        p(d|\theta) = \prod_i^N \frac{\Gamma(n_i+m_i+\frac{1}{2})}
        {2^{n_i+m_i+\frac{1}{2}} n_i!\Gamma(m_i+\frac{1}{2}))}

    .. math::

        \log(p) = \sum_i^N \left[\log\Gamma(n_i+m_i+\frac{1}{2})
        - (m_i+n_i+\frac{1}{2})\log2 -\log n_i!
        - \log \Gamma(m_i+\frac{1}{2}) \right]

    Minus logarithm:

    .. math::

        \log(p) = \sum_i^N \left[\log\Gamma(n_i+m_i+\frac{1}{2})-
        \log \Gamma(m_i+\frac{1}{2}) \right]
        - 0.693  (M+N+\frac{1}{2}) - \sum_i^N \log n_i!

    .. math::

        \log(p) = SumLogGamma(n_i, m_i) -0.693 (N+\frac{1}{2}) -
        \sum_i^N \log n_i! - 0.693\,M

    .. math::

        \log(p) = f(n_i) + SumLogGamma(n_i, m_i) - 0.693\,M

    .. math::

        \log(p)\approx SumLogGamma(n_i, m_i) - 0.693\,M

    :param ranges: Per-dimension ranges.
    :type ranges: list
    :param Nbins: Per-dimension total number of bins
    :type Nbins: list
    :param cl_z_idx: Index of bins where the number of stars is not 0
    :type cl_z_idx: np.ndarray
    :param cl_histo_f_z: Flattened observed Hess diagram with the empty bins removed
    :type cl_histo_f_z: np.ndarray
    :param max_lkl: Maximum likelihood value, used for normalization
    :type max_lkl: float
    :param synth_clust: Synthetic cluster data.
    :type synth_clust: np.ndarray

    :return: Log likelihood value.
    :rtype: float
    """
    # If synthetic cluster is empty, assign a small likelihood value.
    if not synth_clust.any():
        return -1.0e09

    # Obtain histogram for the synthetic cluster.
    mag, colors = synth_clust[0], synth_clust[1:]
    syn_histo_f = []
    for i, col in enumerate(colors):
        # hess_diag = np.histogram2d(
        #     mag, col, bins=[
        #         self.bin_edges[0]] + [self.bin_edges[i + 1]])[0]
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

    tremmel_lkl = np.sum(
        loggamma(cl_histo_f_z + syn_histo_f_z + 0.5) - loggamma(syn_histo_f_z + 0.5)
    )

    # M = syn_histo_f_z.sum() <-- This is wrong and even more, this term should not
    # have been present at all
    # tremmel_lkl = SumLogGamma #- 0.693 * syn_histo_f_z.sum()  # ln(2) ~ 0.693

    return 1 - tremmel_lkl / max_lkl


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

import numpy as np
from scipy.special import loggamma
from fast_histogram import histogram2d


def tremmel(cluster_dict, synth_clust):
    r"""
    Poisson likelihood ratio as defined in Tremmel et al (2013), Eq 10 with
    v_{i,j}=1. This returns the log likelihood.

    p(d|\theta) = \prod_i^N \frac{\Gamma(n_i+m_i+\frac{1}{2})}
        {2^{n_i+m_i+\frac{1}{2}} n_i!\Gamma(m_i+\frac{1}{2}))}

    \log(p) = \sum_i^N \left[\log\Gamma(n_i+m_i+\frac{1}{2})
        - (m_i+n_i+\frac{1}{2})\log2 -\log n_i!
        - \log \Gamma(m_i+\frac{1}{2}) \right]

    Minus logarithm:

    \log(p) = \sum_i^N \left[\log\Gamma(n_i+m_i+\frac{1}{2})-
        \log \Gamma(m_i+\frac{1}{2}) \right]
        - 0.693  (M+N+\frac{1}{2}) - \sum_i^N \log n_i!

    \log(p) = SumLogGamma(n_i, m_i) -0.693 (N+\frac{1}{2}) -
        \sum_i^N \log n_i! - 0.693\,M

    \log(p) = f(n_i) + SumLogGamma(n_i, m_i) - 0.693\,M

    \log(p)\approx SumLogGamma(n_i, m_i) - 0.693\,M

    """

    # If synthetic cluster is empty, assign a small likelihood value.
    if not synth_clust.any():
        return -1.0e09

    #
    mag, colors = synth_clust[0], synth_clust[1:]

    # Obtain histogram for observed cluster.
    syn_histo_f = []
    for i, col in enumerate(colors):
        # hess_diag = np.histogram2d(
        #     mag, col, bins=[
        #         cluster_dict['bin_edges'][0]] + [cluster_dict['bin_edges'][i + 1]])[0]

        hess_diag = histogram2d(
            mag,
            col,
            range=[
                [cluster_dict["ranges"][0][0], cluster_dict["ranges"][0][1]],
                [cluster_dict["ranges"][i + 1][0], cluster_dict["ranges"][i + 1][1]],
            ],
            bins=[cluster_dict["Nbins"][0], cluster_dict["Nbins"][i + 1]],
        )

        # Flatten array
        syn_histo_f += list(hess_diag.ravel())
    syn_histo_f = np.array(syn_histo_f)

    # Remove all bins where n_i = 0 (no observed stars).
    syn_histo_f_z = syn_histo_f[cluster_dict["cl_z_idx"]]

    SumLogGamma = np.sum(
        loggamma(cluster_dict["cl_histo_f_z"] + syn_histo_f_z + 0.5)
        - loggamma(syn_histo_f_z + 0.5)
    )

    # M = syn_histo_f_z.sum()
    # ln(2) ~ 0.693
    tremmel_lkl = SumLogGamma - 0.693 * syn_histo_f_z.sum()

    return tremmel_lkl


def visual(cluster_dict, synth_clust):
    # If synthetic cluster is empty, assign a small likelihood value.
    if not synth_clust.any():
        return -1.0e09

    mag_o, colors_o = cluster_dict['mag'], cluster_dict['colors']
    mag_s, colors_s = synth_clust[0], synth_clust[1:]

    N_mag, N_col = 15, 10
    mag = list(mag_o) + list(mag_s)
    col = list(colors_o[0]) + list(colors_s[0])
    mag_min, mag_max = np.nanmin(mag), np.nanmax(mag)
    bin_edges = [np.linspace(mag_min, mag_max, N_mag)]
    col_min, col_max = np.nanmin(col), np.nanmax(col)
    bin_edges.append(np.linspace(col_min, col_max, N_col))

    # Obtain histogram for observed cluster.
    cl_histo_f = []
    for i, col_o in enumerate(colors_o):
        hess_diag = np.histogram2d(mag_o, col_o, bins=bin_edges)[0]
        # Flatten array
        cl_histo_f += list(hess_diag.ravel())
    cl_histo_f = np.array(cl_histo_f)
    # Down sample histogram
    # msk = cl_histo_f > 5
    # cl_histo_f[msk] = 5

    syn_histo_f = []
    for i, col_s in enumerate(colors_s):
        hess_diag = np.histogram2d(mag_s, col_s, bins=bin_edges)[0]
        # Flatten array
        syn_histo_f += list(hess_diag.ravel())
    syn_histo_f = np.array(syn_histo_f)
    # Down sample histogram
    # msk = syn_histo_f > 5
    # syn_histo_f[msk] = 5

    # return -sum(abs(cluster_dict["hist_down_samp"]-syn_histo_f))

    # create a mask for where each data set is non-zero
    m1 = cl_histo_f != 0
    m2 = syn_histo_f != 0
    m1_area = m1.sum()
    m2_area = m2.sum()
    tot_area = m1_area + m2_area
    # use a logical and to create a combined map where both datasets are non-zero
    ovrlp_area = np.logical_and(m1, m2).sum()

    # # calculate the overlapping density, where 0.5 is the bin width
    # ol_density = np.abs((cl_histo_f - syn_histo_f) * 0.5)[ol]
    # # calculate the total overlap percent
    # h_overlap = ol_density.sum() * 100

    h_overlap = ovrlp_area / tot_area

    # ol_density = np.abs((cl_histo_f - cl_histo_f) * 0.5)[ol]
    # h_overlap_0 = ol_density.sum() * 100
    # print(h_overlap_0)
    # breakpoint()

    # # cl_histo_f = cluster_dict["hist_down_samp"]
    # mi_cnst = np.clip(cl_histo_f, a_min=0, a_max=1)
    # # Final chi.
    # mig_chi = np.sum((cl_histo_f + mi_cnst - syn_histo_f)**2 / (cl_histo_f + 1.))
    # mig_chi_0 = np.sum((cl_histo_f + mi_cnst)**2 / (cl_histo_f + 1.))
    # mig_chi_opt = np.sum((cl_histo_f + mi_cnst - cl_histo_f)**2 / (cl_histo_f + 1.))
    # print(mig_chi_opt, mig_chi_0, mig_chi)

    # import matplotlib.pyplot as plt
    # # plt.subplot(121)
    # y_edges, x_edges = bin_edges #cluster_dict['bin_edges']
    # for xe in x_edges:
    #     plt.axvline(xe, c="grey", ls=":")
    # for ye in y_edges:
    #     plt.axhline(ye, c="grey", ls=":")
    # plt.title(h_overlap)
    # plt.scatter(colors_o[0], mag_o, alpha=.5)
    # plt.scatter(colors_s[0], mag_s, alpha=.5)
    # plt.gca().invert_yaxis()

    # # plt.subplot(122)
    # # plt.title(h_overlap)
    # # plt.bar(np.arange(len(cl_histo_f)), cl_histo_f, label='obs', alpha=.5)
    # # plt.bar(np.arange(len(syn_histo_f)), syn_histo_f, label='syn', alpha=.5)
    # # plt.legend()
    # plt.show()

    return h_overlap


def mean_dist(cluster_dict, synth_clust):
    # If synthetic cluster is empty, assign a small likelihood value.
    if not synth_clust.any():
        return -1.0e09

    # mag_o, colors_o = cluster_dict['mag'], cluster_dict['colors']
    mag0, colors0 = cluster_dict['mag0'], cluster_dict['colors0']
    mag_s, colors_s = synth_clust[0], synth_clust[1:]

    dist = (np.median(mag0) - np.median(mag_s))**2 + (np.median(colors0) - np.median(colors_s))**2
    return -dist

    if len(mag_s) < 5:
        return -1.0e09
    import ndtest
    P_val = ndtest.ks2d2s(mag0, colors0, mag_s, colors_s[0])

    # import matplotlib.pyplot as plt
    # plt.title(P_val)
    # plt.scatter(colors_o0, mag_o, alpha=.5)
    # plt.scatter(colors_s[0], mag_s, alpha=.5)
    # plt.gca().invert_yaxis()
    # plt.show()

    return P_val


def bins_distance(cluster_dict, synth_clust):
    """ Sum of distances to corresponding bins in he Hess diagram."""
    if not synth_clust.any():
        return 1.0e09

    mag0, colors0 = cluster_dict['mag'], cluster_dict['colors'][0]
    msk = np.isnan(mag0) | np.isnan(colors0)
    mag0, colors0 = mag0[~msk], colors0[~msk]
    mag_s, colors_s = synth_clust[0], synth_clust[1]

    mpercs = (.1, .5, 1, 2, 5, 10, 30, 60, 90)
    cpercs = (.5, 1, 5, 15, 30, 60, 90)

    perc_mag0 = np.percentile(mag0, mpercs)
    perc_colors0 = np.percentile(colors0, cpercs)
    pts_0 = []
    for pm in perc_mag0:
        for pc in perc_colors0:
            pts_0.append([pm, pc])
    pts_0 = np.array(pts_0).T

    # import matplotlib.pyplot as plt
    # plt.scatter(colors0, mag0)
    # plt.scatter(pts_0[1], pts_0[0], c='k')
    # plt.gca().invert_yaxis()
    # plt.show()
    # breakpoint()

    perc_mag_s = np.percentile(mag_s, mpercs)
    perc_colors_s = np.percentile(colors_s, cpercs)
    pts_s = []
    for pm in perc_mag_s:
        for pc in perc_colors_s:
            pts_s.append([pm, pc])
    pts_s = np.array(pts_s).T

    dist = np.sqrt((pts_s[0] - pts_0[0])**2 + (pts_s[1] - pts_0[1])**2)
    weights = np.linspace(1, .05, len(mpercs)*len(cpercs))
    lkl = sum(dist*weights)

    return lkl

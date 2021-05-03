
import numpy as np
from scipy.special import logsumexp, loggamma
from scipy.stats import gaussian_kde, entropy

# ############################################################
# # Timer function: http://stackoverflow.com/a/21860100/1391441
# from contextlib import contextmanager
# import time


# @contextmanager
# def timeblock(label):
#     start = time.clock()
#     try:
#         yield
#     finally:
#         end = time.clock()
#         print ('{} elapsed: {}'.format(label, end - start))
# ############################################################


def main(lkl_method, synth_clust, obs_clust):
    """
    Match the synthetic cluster to the observed cluster.
    """

    # If synthetic cluster is empty, assign a large likelihood value. This
    # assumes *all* likelihoods here need minimizing.
    if not synth_clust.any():
        return 1.e09

    # Obtain the likelihood matching the synthetic and observed clusters.
    if lkl_method == 'tremmel':
        likelihood = tremmel(synth_clust, obs_clust)
    elif lkl_method == 'dolphin':
        likelihood = dolphin(synth_clust, obs_clust)
    elif lkl_method == 'mighell':
        likelihood = mighell(synth_clust, obs_clust)
    elif lkl_method == 'tolstoy':
        likelihood = tolstoy(synth_clust, obs_clust)
    elif lkl_method == 'isochfit':
        likelihood = isochfit(synth_clust, obs_clust)
    elif lkl_method == 'dolphin_kde':
        likelihood = dolphin_kde(synth_clust, obs_clust)
    elif lkl_method == 'kdeKL':
        likelihood = kdeKL(synth_clust, obs_clust)

    return likelihood


def tremmel(synth_clust, obs_clust):
    """
    Poisson likelihood ratio as defined in Tremmel et al (2013), E1 10 with
    v_{i,j}=1. This returns the negative log likelihood.

    p(d|\theta) = \prod_i^N \frac{\Gamma(n_i+m_i+\frac{1}{2})}
        {2^{n_i+m_i+\frac{1}{2}} d_i!\Gamma(m_i+\frac{1}{2}))}

    \log(p) = \sum_i^N \left[\log\Gamma(n_i+m_i+\frac{1}{2})
        - (m_i+n_i+\frac{1}{2})\log2 -\log n_i!
        - \log \Gamma(m_i+\frac{1}{2}) \right]

    Minus logarithm:

    -\log(p) = 0.693  (M+N+\frac{1}{2}) + \sum_i^N \log n_i! -
        \sum_i^N \left[\log\Gamma(n_i+m_i+\frac{1}{2})-
        \log \Gamma(m_i+\frac{1}{2}) \right]

    -\log(p) = 0.693 (N+\frac{1}{2}) + \sum_i^N \log n_i! +
        0.693\,M - SumLogGamma(n_i, m_i)

    -\log(p) = f(n_i) + 0.693\,M - SumLogGamma(n_i, m_i)

    -\log(p)\approx 0.693\,M - SumLogGamma(n_i, m_i)

    """

    # Observed cluster's data.
    bin_edges, cl_histo_f_z, cl_z_idx = obs_clust

    # Histogram of the synthetic cluster, using the bin edges calculated
    # with the observed cluster.
    syn_histo = np.histogramdd(synth_clust, bins=bin_edges)[0]
    # Flatten N-dimensional histogram.
    syn_histo_f = syn_histo.ravel()
    # Remove all bins where n_i = 0 (no observed stars).
    syn_histo_f_z = syn_histo_f[cl_z_idx]

    SumLogGamma = np.sum(
        loggamma(cl_histo_f_z + syn_histo_f_z + .5) -
        loggamma(syn_histo_f_z + .5))

    # M = synth_clust.shape[0]
    # ln(2) ~ 0.693
    tremmel_lkl = 0.693 * synth_clust.shape[0] - SumLogGamma

    return tremmel_lkl


def dolphin(synth_clust, obs_clust):
    """
    Poisson likelihood ratio as defined in Dolphin (2002).

    -2\ln PLR = 2 \sum_i m_i - n_i + n_i \ln \frac{n_i}{m_i}
              = 2 (M- N) + 2\sum_i n_i\ln n_i - 2\sum_i n_i\ln m_i

    If the binning is made too small then  n_i, m_i --> 1 (one star per bin)
    and thus:

    -2\ln PLR --> 2*(M-N)

    In this case the likelihood will try to minimize M.

    When the number of observed stars is too small (~ N<100), this likelihood
    will tend to underestimate the total mass.
    """

    # Observed cluster's data.
    bin_edges, fill_factor, cl_histo_f_z, dolphin_cst, cl_z_idx = obs_clust

    # Histogram of the synthetic cluster, using the bin edges calculated
    # with the observed cluster.
    syn_histo = np.histogramdd(synth_clust, bins=bin_edges)[0]
    # Flatten N-dimensional histogram.
    syn_histo_f = syn_histo.ravel()
    # Remove all bins where n_i = 0 (no observed stars).
    syn_histo_f_z = syn_histo_f[cl_z_idx]

    # Assign small value to the m_i = 0 elements in 'syn_histo_f_z'.
    # The value equals 1 star divided among all empty bins.
    # If this factor --> 0, then the likelihood will try to minimize empty
    # bins, ie: M --> inf. If the factor --> 1, the likelihood has no penalty
    # for empty bins, and M --> 0.
    syn_histo_f_z[syn_histo_f_z == 0] = fill_factor

    # M = synth_phot[0].size
    # Cash's C statistic: 2 * sum(m_i - n_i * ln(m_i))
    C_cash = 2. * (synth_clust.shape[0] - np.sum(
        (cl_histo_f_z * np.log(syn_histo_f_z))))

    # Obtain (weighted) inverse logarithmic 'Poisson likelihood ratio'.
    dolph_lkl = C_cash + dolphin_cst

    return dolph_lkl


def mighell(synth_clust, obs_clust):
    """
    Chi gamma squared distribution defined in Mighell (1999)

    This likelihood is more stable than Dolphin regarding the issue of empty
    bins, but it also has less power to discriminate lower masses from the
    actual mass.

    If the number of bins is too large, it will attempt to minimize the
    synthetic cluster mass M. This is because in the infinite bins limits,
    each bin holds a single star and the chances of n_i=m_i go to zero. In
    this case, the chi-square is lowered simply lowering M.
    """

    # Observed cluster's bin edges for each dimension, flattened histogram,
    # and n_i constant.
    bin_edges, cl_histo_f, ni_cnst = obs_clust

    # Histogram of the synthetic cluster, using the bin edges calculated
    # with the observed cluster.
    syn_histo = np.histogramdd(synth_clust, bins=bin_edges)[0]
    # Flatten N-dimensional histogram.
    syn_histo_f = syn_histo.ravel()

    # Final chi.
    mig_chi = np.sum(
        (np.square(syn_histo_f) + ni_cnst * syn_histo_f) / (cl_histo_f + 1.))

    return mig_chi


def tolstoy(synth_clust, obs_clust):
    """
    Weighted (log) likelihood.
    This function follows the recipe given in Tolstoy & Saha (1996).

    Likelihood form:

    -\ln L = N\ln M -
             \sum\limits_{i=1}^{N} \ln
             \left\{
                P_i \sum\limits_{j=1}^M \frac{\exp
                    \left[
                        -\frac{1}{2}
                                \sum_{k=1}^D \left(\frac{f_{ik}-g_{jk}}
                                {\sigma_{ik}}
                            \right)^2
                    \right]}
                {\prod_{k=1}^D \sigma_{ik}}
            \right\}

    """
    # DEPRECATED 18/01/20
    # def old(obs_st, N, log_mem_probs, synth_phot, synth_errors):
    #     """
    #     -\ln L = N\ln M -
    #              \sum\limits_{i=1}^{N} \ln
    #              \left\{
    #                 P_i \sum\limits_{j=1}^M \frac{\exp
    #                     \left[
    #                         -\frac{1}{2}
    #                             \left(
    #                                 \sum_{k=1}^D \frac{(f_{ik}-g_{jk})^2}
    #                                 {\sigma_{ik}^2 + \rho_{jk}^2}
    #                             \right)
    #                     \right]}
    #                 {\prod_{k=1}^D \sqrt{\sigma_{ik}^2 + \rho_{jk}^2}}
    #             \right\}
    #     """
    #     # Square synthetic photometric errors.
    #     synth_errors = np.square(synth_errors)
    #     # Array with the proper format for the synthetic cluster.
    #     syn_st = np.dstack([np.array(synth_phot).T, synth_errors.T])

    #     # Photometric difference (observed - synthetic), for all dimensions.
    #     phot_dif = obs_st[:, None, :, 0] - syn_st[None, :, :, 0]
    #     # Sum of squared photometric errors, for all dimensions. Clip at a
    #     # minimum of 0.005 to avoid numeric issues below.
    #     sigma_sum = np.clip(
    #         obs_st[:, None, :, 1] + syn_st[None, :, :, 1], 0.005, None)

    #     # Sum for all photometric dimensions.
    #     Dsum = (np.square(phot_dif) / sigma_sum).sum(axis=-1)
    #     # Product of summed squared sigmas.
    #     sigma_prod = np.prod(sigma_sum, axis=-1)

    #     # The block below can be replaced by this line using 'logsumexp'. It
    #     # is marginally faster.
    #     sum_N = (logsumexp(-0.5 * Dsum, b=1. / np.sqrt(sigma_prod), axis=1) +
    #              log_mem_probs).sum()

    #     # # All elements inside synthetic stars summatory.
    #     # sum_M_j = np.exp(-0.5 * Dsum) / np.sqrt(sigma_prod)
    #     # # Sum for all synthetic stars.
    #     # sum_M = np.sum(sum_M_j, axis=-1)
    #     # # Multiply by membership probabilities.
    #     # sum_M_MP = np.exp(log_mem_probs) * sum_M
    #     # # Replace 0. elements before applying the logarithm below.
    #     # sum_M_MP[sum_M_MP == 0.] = 1e-7
    #     # sum_N0 = np.sum(np.log(sum_M_MP))

    #     # Final negative logarithmic likelihood
    #     tlst_lkl = N * np.log(len(syn_st)) - sum_N

    #     return tlst_lkl

    # # synthetic cluster's photometry and errors.
    # synth_phot = synth_clust.T
    # synth_errors = np.zeros((synth_phot.shape))
    # obs_st, N, log_mem_probs = obs_clust[0]
    # lkl = old(obs_st, N, log_mem_probs, synth_phot, synth_errors)

    # Observed cluster's photometry and membership probabilities.
    obs_photom, sigma, sigma_prod, N, log_mem_probs = obs_clust

    # This line takes up ~30% of the processing time
    # Sum for all photometric dimensions.
    Dsum = (np.square(
        obs_photom - synth_clust[None, :, :]) / sigma).sum(axis=-1)

    # This line takes up ~70% of the processing time
    sum_N = (
        logsumexp(-.5 * Dsum, b=1. / sigma_prod, axis=1) +
        log_mem_probs).sum()

    # Final negative logarithmic likelihood
    tlst_lkl = N * np.log(synth_clust.shape[0]) - sum_N

    return tlst_lkl


def isochfit(synth_clust, obs_clust):
    """
    In place for #358

    This is not trivial to implement.
    """

    from scipy.stats import energy_distance, wasserstein_distance
    # from scipy.spatial.distance import cdist

    clst_histo, bin_edges = obs_clust

    # Histogram of the synthetic cluster, using the bin edges calculated
    # with the observed cluster.
    syn_histo = np.histogramdd(synth_clust, bins=bin_edges)[0]
    # Flatten
    syn_histo = syn_histo.ravel()
    # Downsample
    # syn_histo[syn_histo > 5] = 5
    # The above is faster?
    # syn_histo = np.clip(syn_histo, a_min=None, a_max=1)

    lkl = energy_distance(clst_histo, syn_histo)

    # # Convert histogram into distribution of points
    # lkl = 1.e09
    # if syn_histo.sum() > 0.:
    #     idxs = np.array(np.where(syn_histo > 0)).T
    #     synth_dwnsmp = []
    #     for pt in idxs:
    #         coord = []
    #         for i, j in enumerate(pt):
    #             coord.append(bin_edges[i][j])
    #         synth_dwnsmp.append(coord)
    #     synth_dwnsmp = np.array(synth_dwnsmp)

    #     # lkl = anderson_ksamp([syn_histo_f, cl_histo_f])[0]

    #     # import matplotlib.pyplot as plt
    #     # print(np.mean(cdist(clust_dwnsmp, synth_dwnsmp)))
    #     # clust_dwnsmp = np.array(clust_dwnsmp).T
    #     # synth_dwnsmp = np.array(synth_dwnsmp).T
    #     # plt.scatter(clust_dwnsmp[1], clust_dwnsmp[0], c='g')
    #     # plt.scatter(synth_dwnsmp[1], synth_dwnsmp[0], c='r')
    #     # plt.gca().invert_yaxis()
    #     # plt.show()

    #     # Tried and none worked:
    #     # np.mean((cdist(clust_dwnsmp, synth_dwnsmp))
    #     # np.sum((cdist(clust_dwnsmp, synth_dwnsmp))
    #     # np.sum(np.mean(cdist(clust_dwnsmp, synth_dwnsmp), 1))

    #     lkl = np.sum(np.mean(cdist(clust_dwnsmp, synth_dwnsmp), 1))

    return lkl


# DEPRECATED 17/01/20
# def duong(synth_clust, obs_clust):
#     """
#     """
#     import rpy2.robjects as robjects
#     from rpy2.rinterface import RRuntimeError

#     # synthetic cluster's photometry and errors.
#     synth_phot, synth_errors = synth_clust[0]
#     # Observed cluster's photometry and membership probabilities.
#     kde_test, hpi_kfe, m_cl, hpic = obs_clust

#     # CMD for synthetic cluster.
#     matrix_f1 = np.ravel(np.column_stack((synth_phot)))
#     rows_f1 = int(len(matrix_f1) / 2)

#     m_f1 = robjects.r.matrix(robjects.FloatVector(matrix_f1),
#                              nrow=rows_f1, byrow=True)

#     try:
#         # TODO this is the second line that takes the most time.
#         # hpif1 = hpi_kfe(x=m_f1, binned=True)

#         # Call 'ks' function to obtain p_value.
#         # TODO: this line takes forever
#         # TODO: this statistic seems to select lower masses

#         # Should I:
#         # 1. explicit different bandwidths with H1,H2?
#         # res_cl = kde_test(x1=m_cl, x2=m_f1, H1=hpic, H2=hpif1)
#         # 2. use the same bandwidth defined for the cluster (hpic)?
#         # res_cl = kde_test(x1=m_cl, x2=m_f1, H1=hpic, H2=hpic)
#         # 3. not explicit any bandwidth?
#         res_cl = kde_test(x1=m_cl, x2=m_f1)
#         p_val_cl = res_cl.rx2('pvalue')
#         # Store cluster vs field p-value.
#         duong_pval = 1. - p_val_cl[0]
#     except RRuntimeError:
#         duong_pval = 1.

#     return duong_pval


def dolphin_kde(synth_clust, obs_clust):
    '''
    Poisson likelihood ratio as defined in Dolphin (2002).

    -2\ln PLR = 2 \sum_i m_i - n_i + n_i \ln \frac{n_i}{m_i}
              = 2 (M- N) + 2\sum_i n_i\ln n_i - 2\sum_i n_i\ln m_i

    If the binning is made too small then  n_i, m_i --> 1 (one star per bin)
    and thus:

    -2\ln PLR --> 2*(M-N)

    In this case the likelihood will try to minimize M.

    '''

    obs_kde, kde_pts = obs_clust

    try:
        # synth_clust.shape = (# of dims, # of data)
        kernel = gaussian_kde(synth_clust[0][0])
        synth_kde = kernel(kde_pts)

        # Assign small value to the m_i = 0 elements in 'syn_histo_f_z'.
        # The value equals 1 star divided among all empty bins.
        synth_kde[synth_kde == 0] = 1. / max(
            np.count_nonzero(synth_kde == 0), 1.)

        # Cash's C statistic.
        # M = synth_clust[0][0][0].size
        # C_cash = M - np.sum(N * obs_kde * np.log(M * synth_kde))
        C_cash = np.sum(synth_kde - obs_kde * np.log(synth_kde))

        # Obtain (weighted) inverse logarithmic 'Poisson likelihood ratio'.
        # dolph_lkl = min(C_cash, 1e09)  # + dolphin_cst
        dolph_lkl = C_cash\
            if not np.isnan(C_cash) and not np.isinf(C_cash) else 1e09
    except:
        dolph_lkl = 1e09

    return dolph_lkl


def kdeKL(synth_clust, obs_clust):
    """
    Kullback-Leibler divergence between the N-dimensional KDE of the observed
    cluster versus the synthetic cluster.
    """
    # TODO seems to work, but it is terribly slow.
    obs_kde, kde_pts = obs_clust

    try:
        # synth_clust.shape = (# of dims, # of data)
        kernel = gaussian_kde(synth_clust[0][0])
        synth_kde = kernel(kde_pts)
        # if (synth_kde == 0.).all():
        #     print('zero')
        # if np.isnan(synth_kde).all():
        #     print('nan')
        # if np.isinf(synth_kde).all():
        #     print("inf")
        kl = min(entropy(obs_kde, synth_kde), 1000.)
        kl = kl if not np.isnan(kl) else 1000.
    except:
        kl = 1e09

    print(kl)
    return kl


# def entropy(pk, qk=None, base=None):
#     """Calculate the entropy of a distribution for given probability values.
#     If only probabilities `pk` are given, the entropy is calculated as
#     ``S = -sum(pk * log(pk), axis=0)``.
#     If `qk` is not None, then compute the Kullback-Leibler divergence
#     ``S = sum(pk * log(pk / qk), axis=0)``.
#     This routine will normalize `pk` and `qk` if they don't sum to 1.
#     Parameters
#     ----------
#     pk : sequence
#         Defines the (discrete) distribution. ``pk[i]`` is the (possibly
#         unnormalized) probability of event ``i``.
#     qk : sequence, optional
#         Sequence against which the relative entropy is computed. Should be in
#         the same format as `pk`.
#     base : float, optional
#         The logarithmic base to use, defaults to ``e`` (natural logarithm).
#     Returns
#     -------
#     S : float
#         The calculated entropy.
#     """
#     from scipy.special import rel_entr
#     # pk = np.asarray(pk)
#     pk = 1. * pk / np.sum(pk, axis=0)
#     # if qk is None:
#     #     vec = entr(pk)
#     # else:
#     # qk = np.asarray(qk)
#     # if len(qk) != len(pk):
#     #     raise ValueError("qk and pk must have same length.")
#     qk0 = 1. * qk / np.sum(qk, axis=0)
#     vec = rel_entr(pk, qk0)
#     S = np.sum(vec, axis=0)
#     # if base is not None:
#     #     S /= log(base)
#     return S

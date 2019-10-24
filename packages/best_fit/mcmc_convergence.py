
import numpy as np
from scipy.stats import chi2
from scipy.special import gammaln
# from scipy.signal import fftconvolve
import logging
import warnings
from emcee import autocorr
from .ptemcee import util


def fminESS(p, alpha=.05, eps=.05, ess=None):
    """
    Minimum effective sample size
    """

    crit = chi2.ppf(1 - alpha, p)
    foo = 2. / p

    if ess is None:
        logminESS = foo * np.log(2.) + np.log(np.pi) - foo * np.log(p) -\
            foo * gammaln(p / 2.) - 2. * np.log(eps) + np.log(crit)
        return np.round(np.exp(logminESS))
    else:
        if isinstance(ess, str):
            raise ValueError("Only numeric entry allowed for ess")
        try:
            logEPS = .5 * foo * np.log(2.) + .5 * np.log(np.pi) -\
                .5 * foo * np.log(p) - .5 * foo * gammaln(p / 2.) -\
                .5 * np.log(ess) + .5 * np.log(crit)
        except FloatingPointError:
            return np.nan
        return np.exp(logEPS)


def multiESS(X, b='sqroot', Noffsets=10, Nb=None):
    """
    Compute multivariate effective sample size of a single Markov chain X,
    using the multivariate dependence structure of the process.

    X: MCMC samples of shape (n, p)
    n: number of samples
    p: number of parameters

    b: specifies the batch size for estimation of the covariance matrix in
       Markov chain CLT. It can take a numeric value between 1 and n/2, or a
       char value between:

    'sqroot'    b=floor(n^(1/2)) (for chains with slow mixing time; default)
    'cuberoot'  b=floor(n^(1/3)) (for chains with fast mixing time)
    'lESS'      pick the b that produces the lowest effective sample size
                for a number of b ranging from n^(1/4) to n/max(20,p); this
                is a conservative choice

    If n is not divisible by b Sigma is recomputed for up to Noffsets subsets
    of the data with different offsets, and the output mESS is the average over
    the effective sample sizes obtained for different offsets.

    Nb specifies the number of values of b to test when b='less'
    (default NB=200). This option is unused for other choices of b.

    Original source: https://github.com/lacerbi/multiESS

    Reference:
    Vats, D., Flegal, J. M., & Jones, G. L. "Multivariate Output Analysis
    for Markov chain Monte Carlo", arXiv preprint arXiv:1512.07713 (2015).

    """

    # Trying to implement the changes mentioned in
    # https://stats.stackexchange.com/a/359474/10416
    # to work with non-independent chains. I'm getting an ESS that is
    # significantly larger than both PyMC3's and (nsteps*nwalkers)/acorr_t

    # MCMC samples and parameters
    w = X.shape[1]
    X_w_avrg = np.mean(X, axis=1)
    n, p = X_w_avrg.shape
    X_concat = X.reshape(-1, p)

    if p > n:
        raise ValueError(
            "More dimensions than data points, cannot compute effective "
            "sample size.")

    # Input check for batch size B
    if isinstance(b, str):
        if b not in ['sqroot', 'cuberoot', 'less']:
            raise ValueError(
                "Unknown string for batch size. Allowed arguments are "
                "'sqroot', 'cuberoot' and 'lESS'.")
        if b != 'less' and Nb is not None:
            raise Warning(
                "Nonempty parameter NB will be ignored (NB is used "
                "only with 'lESS' batch size B).")
    else:
        if not 1. < b < (n / 2):
            raise ValueError(
                "The batch size B needs to be between 1 and N/2.")

    # Compute multiESS for the chain
    mESS = multiESS_chain(X_w_avrg, X_concat, n, w, p, b, Noffsets, Nb)

    return mESS


def multiESS_chain(X_w_avrg, X_concat, n, w, p, b, Noffsets, Nb):
    """
    Compute multiESS for a MCMC chain.
    """

    # Determinant of sample covariance matrix
    if p == 1:
        detLambda = np.cov(X_concat.T)
    else:
        detLambda = np.linalg.det(np.cov(X_concat.T))

    # Parameters sample mean
    theta = np.mean(X_w_avrg, axis=0)

    if b == 'sqroot':
        b = [int(np.floor(n ** (1. / 2)))]
    elif b == 'cuberoot':
        b = [int(np.floor(n ** (1. / 3)))]
    elif b == 'less':
        b_min = np.floor(n ** (1. / 4))
        b_max = max(np.floor(n / max(p, 20)), np.floor(np.sqrt(n)))
        if Nb is None:
            Nb = 200
        # Try NB log-spaced values of B from B_MIN to B_MAX
        b = set(map(int, np.round(np.exp(
            np.linspace(np.log(b_min), np.log(b_max), Nb)))))

    # Compute mESS
    mESS_i = []
    for bi in b:
        mESS_i.append(multiESS_batch(
            X_w_avrg, n, w, p, theta, detLambda, bi, Noffsets))

    # Return lowest mESS
    if np.isnan(np.array(mESS_i)).all():
        mESS = np.nan
    else:
        mESS = np.nanmin(mESS_i)

    return mESS


def multiESS_batch(X_w_avrg, n, w, p, theta, detLambda, bi, Noffsets):
    """
    Compute multiESS for a given batch size B.
    """

    # Compute batch estimator for SIGMA
    a = int(np.floor(n / bi))
    Sigma = np.zeros((p, p))
    offsets = np.sort(list(set(map(int, np.round(
        np.linspace(0, n - np.dot(a, bi), Noffsets))))))

    for j in offsets:
        # Swapped a, bi in reshape compared to the original code.
        Y = X_w_avrg[j + np.arange(a * bi), :].reshape((a, bi, p))
        Ybar = np.squeeze(np.mean(Y, axis=1))
        Z = Ybar - theta
        for i in range(a):
            if p == 1:
                Sigma += Z[i] ** 2
            else:
                Sigma += Z[i][np.newaxis, :].T * Z[i]

    Sigma = (Sigma * bi) / (a - 1) / len(offsets)

    detSigma = np.linalg.det(w * Sigma)
    if detLambda > 1e-10 and detSigma > 1.e-10:
        mESS = n * w * (detLambda / detSigma) ** (1. / p)
    else:
        mESS = np.nan

    return mESS


def geweke(x, first=.1, last=.5, intervals=20):
    R"""
    Source: https://github.com/pymc-devs/pymc3/blob/master/pymc3/diagnostics.py

    Return z-scores for convergence diagnostics.
    Compare the mean of the first % of series with the mean of the last % of
    series. x is divided into a number of segments for which this difference is
    computed. If the series is converged, this score should oscillate between
    -1 and 1.
    Parameters
    ----------
    x : array-like
      The trace of some stochastic parameter.
    first : float
      The fraction of series at the beginning of the trace.
    last : float
      The fraction of series at the end to be compared with the section
      at the beginning.
    intervals : int
      The number of segments.
    Returns
    -------
    scores : list [[]]
      Return a list of [i, score], where i is the starting index for each
      interval and score the Geweke score on the interval.
    Notes
    -----
    The Geweke score on some series x is computed by:
      .. math:: \frac{E[x_s] - E[x_e]}{\sqrt{V[x_s] + V[x_e]}}
    where :math:`E` stands for the mean, :math:`V` the variance,
    :math:`x_s` a section at the start of the series and
    :math:`x_e` a section at the end of the series.
    References
    ----------
    Geweke (1992)
    """

    if np.ndim(x) > 1:
        return [geweke(y, first, last, intervals) for y in np.transpose(x)]

    # Filter out invalid intervals
    for interval in (first, last):
        if interval <= 0 or interval >= 1:
            raise ValueError(
                "Invalid intervals for Geweke convergence analysis",
                (first,
                 last))
    if first + last >= 1:
        raise ValueError(
            "Invalid intervals for Geweke convergence analysis",
            (first,
             last))

    # Initialize list of z-scores
    zscores = []

    # Last index value
    end = len(x) - 1

    # Start intervals going up to the <last>% of the chain
    last_start_idx = (1 - last) * end

    # Calculate starting indices
    start_indices = np.arange(0, int(last_start_idx), step=int(
        (last_start_idx) / (intervals - 1)))

    # Loop over start indices
    for start in start_indices:
        # Calculate slices
        first_slice = x[start: start + int(first * (end - start))]
        last_slice = x[int(end - last * (end - start)):]

        z = first_slice.mean() - last_slice.mean()
        try:
            z /= np.sqrt(first_slice.var() + last_slice.var())
        except FloatingPointError:
            z = np.nan

        zscores.append([start, z])

    if intervals is None:
        return np.array(zscores[0])
    else:
        return np.array(zscores)


# def return_intersection(hist_1, hist_2):
#     """
#     The ranges of both histograms must coincide for this function to work.

#     Source: https://mpatacchiola.github.io/blog/2016/11/12/
#             the-simplest-classifier-histogram-intersection.html
#     """
#     minima = np.minimum(hist_1, hist_2)
#     intersection = np.true_divide(np.sum(minima), np.sum(hist_2))
#     return intersection
#
#
# def pdfHalfves(varIdxs, mcmc_trace):
#     """
#     Estimate the difference between the first and second halves of a 20 bins
#     histogram of the flat trace, for each parameter.
#     """
#     mcmc_halves = []
#     for cp in [0, 1, 2, 3, 4, 5]:
#         if cp in varIdxs:
#             c_model = varIdxs.index(cp)
#             h_min, h_max = min(mcmc_trace[c_model]), max(mcmc_trace[c_model])
#             half = int(.5 * len(mcmc_trace[c_model]))
#             # 1st half
#             hist_1 = np.histogram(
#                 mcmc_trace[c_model][:half], bins=20, range=[h_min, h_max])[0]
#             # 2nd half
#             hist_2 = np.histogram(
#                 mcmc_trace[c_model][half:], bins=20, range=[h_min, h_max])[0]

#             mcmc_halves.append(return_intersection(hist_1, hist_2))
#         else:
#             mcmc_halves.append(1.)

#     return mcmc_halves

def convergenceVals(algor, ndim, varIdxs, N_conv, chains_nruns, bi_steps):
    """
    Convergence statistics.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        # Mean Tau across chains, shape: (post-bi steps, ndims)
        x = np.mean(chains_nruns.T, axis=1).T
        tau_autocorr = []
        for j in np.arange(50, x.shape[0], 50):
            # tau.shape: ndim
            tau = util.autocorr_integrated_time(x[:j])
            # Autocorrelation time. Mean across dimensions.
            tau_autocorr.append([bi_steps + j, np.mean(tau)])
        # Add one last point with the entire chain.
        if j < x.shape[0]:
            tau = util.autocorr_integrated_time(x)
            tau_autocorr.append([bi_steps + x.shape[0], np.mean(tau)])
        tau_autocorr = np.array(tau_autocorr).T

        # Autocorrelation time for each parameter, mean across chains.
        if algor == 'emcee':
            acorr_t = autocorr.integrated_time(
                chains_nruns, tol=N_conv, quiet=True)
        elif algor == 'ptemcee':
            x = np.mean(chains_nruns.transpose(1, 0, 2), axis=0)
            acorr_t = util.autocorr_integrated_time(x)
        elif algor == 'abc':
            acorr_t = np.array([np.nan] * ndim)

        # Autocorrelation time for each chain for each parameter.
        logger = logging.getLogger()
        logger.disabled = True
        at = []
        # For each parameter/dimension
        for p in chains_nruns.T:
            at_p = []
            # For each chain for this parameter/dimension
            for c in p:
                # at_p.append(autocorr.integrated_time(c, quiet=True)[0])
                at_p.append(util.autocorr_integrated_time(c))
            at.append(at_p)
        logger.disabled = False

        # IAT for all chains and all parameters.
        all_taus = [item for subl in at for item in subl]

        # # Worst chain: chain with the largest acorr time.
        # max_at_c = [np.argmax(a) for a in at]
        # # Best chain: chain with the smallest acorr time.
        # min_at_c = [np.argmin(a) for a in at]
        # Chain with the closest to the median IAT
        med_at_c = [np.argmin(np.abs(np.median(a) - a)) for a in at]

        # Mean Geweke z-scores and autocorrelation functions for all chains.
        geweke_z, acorr_function = [[] for _ in range(ndim)],\
            [[] for _ in range(ndim)]
        for i, p in enumerate(chains_nruns.T):
            for c in p:
                try:
                    geweke_z[i].append(geweke(c))
                except ZeroDivisionError:
                    geweke_z[i].append([np.nan, np.nan])
                try:
                    # acorr_function[i].append(autocorr.function_1d(c))
                    acorr_function[i].append(util.autocorr_function(c))
                except FloatingPointError:
                    acorr_function[i].append([np.nan])
        # Mean across chains
        geweke_z = np.nanmean(geweke_z, axis=1)
        acorr_function = np.nanmean(acorr_function, axis=1)

        # # Cut the autocorrelation function just after *all* the parameters
        # # have crossed the zero line.
        # try:
        #     lag_zero = max([np.where(_ < 0)[0][0] for _ in acorr_function])
        # except IndexError:
        #     # Could not obtain zero lag
        #     lag_zero = acorr_function.shape[-1]
        # acorr_function = acorr_function[:, :int(lag_zero + .2 * lag_zero)]

        # # Approx IAT
        # lag_iat = 1. + 2. * np.sum(acorr_function, axis=1)
        # print("  Approx (zero lag) IAT: ", lag_iat)

        # Effective Sample Size (per param) = (nsteps / tau) * nchains
        mcmc_ess = (chains_nruns.shape[0] / acorr_t) * chains_nruns.shape[1]

        # TODO fix this function
        # # Minimum effective sample size (ESS), and multi-variable ESS.
        # minESS, mESS = fminESS(ndim), multiESS(chains_nruns)
        # # print("mESS: {}".format(mESS))
        # mESS_epsilon = [[], [], []]
        # for alpha in [.01, .05, .1, .2, .3, .4, .5, .6, .7, .8, .9, .95]:
        #     mESS_epsilon[0].append(alpha)
        #     mESS_epsilon[1].append(fminESS(ndim, alpha=alpha, ess=minESS))
        #     mESS_epsilon[2].append(fminESS(ndim, alpha=alpha, ess=mESS))

    return tau_autocorr, acorr_t, med_at_c, all_taus, geweke_z,\
        acorr_function, mcmc_ess

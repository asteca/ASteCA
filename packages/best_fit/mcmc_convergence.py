
import numpy as np
from scipy.stats import chi2
from scipy.special import gammaln
from scipy.signal import fftconvolve
import logging
from .emcee3rc2 import autocorr
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

    # MCMC samples and parameters
    n, p = X.shape

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
    mESS = multiESS_chain(X, n, p, b, Noffsets, Nb)

    return mESS


def multiESS_chain(Xi, n, p, b, Noffsets, Nb):
    """
    Compute multiESS for a MCMC chain.
    """

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

    # Sample mean
    theta = np.mean(Xi, axis=0)
    # Determinant of sample covariance matrix
    if p == 1:
        detLambda = np.cov(Xi.T)
    else:
        detLambda = np.linalg.det(np.cov(Xi.T))

    # Compute mESS
    mESS_i = []
    for bi in b:
        mESS_i.append(multiESS_batch(Xi, n, p, theta, detLambda, bi, Noffsets))
    # Return lowest mESS
    if np.isnan(np.array(mESS_i)).all():
        mESS = np.nan
    else:
        mESS = np.nanmin(mESS_i)

    return mESS


def multiESS_batch(Xi, n, p, theta, detLambda, b, Noffsets):
    """
    Compute multiESS for a given batch size B.
    """

    # Compute batch estimator for SIGMA
    a = int(np.floor(n / b))
    Sigma = np.zeros((p, p))
    offsets = np.sort(list(set(map(int, np.round(
        np.linspace(0, n - np.dot(a, b), Noffsets))))))

    for j in offsets:
        # Swapped a, b in reshape compared to the original code.
        Y = Xi[j + np.arange(a * b), :].reshape((a, b, p))
        Ybar = np.squeeze(np.mean(Y, axis=1))
        Z = Ybar - theta
        for i in range(a):
            if p == 1:
                Sigma += Z[i] ** 2
            else:
                Sigma += Z[i][np.newaxis, :].T * Z[i]

    Sigma = (Sigma * b) / (a - 1) / len(offsets)

    detSigma = np.linalg.det(Sigma)
    # print(detLambda, detSigma)
    if detLambda > 1e-10 and detSigma > 1.e-10:
        mESS = n * (detLambda / detSigma) ** (1. / p)
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


def effective_n(mtrace, varnames=None, include_transformed=False):
    R"""
    https://docs.pymc.io/api/diagnostics.html#pymc3.diagnostics.effective_n

    Returns estimate of the effective sample size of a set of traces.
    Parameters
    ----------
    mtrace : MultiTrace or trace object
      A MultiTrace object containing parallel traces (minimum 2)
      of one or more stochastic parameters.
    varnames : list
      Names of variables to include in the effective_n report
    include_transformed : bool
      Flag for reporting automatically transformed variables in addition
      to original variables (defaults to False).
    Returns
    -------
    n_eff : dictionary of floats (MultiTrace) or float (trace object)
        Return the effective sample size, :math:`\hat{n}_{eff}`
    Notes
    -----
    The diagnostic is computed by:
    .. math:: \hat{n}_{eff} = \frac{mn}{1 + 2 \sum_{t=1}^T \hat{\rho}_t}
    where :math:`\hat{\rho}_t` is the estimated autocorrelation at lag t, and T
    is the first odd positive integer for which the sum
    :math:`\hat{\rho}_{T+1} + \hat{\rho}_{T+1}` is negative.
    The current implementation is similar to Stan, which uses Geyer's initial
    monotone sequence criterion (Geyer, 1992; Geyer, 2011).
    References
    ----------
    Gelman et al. BDA (2014)"""

    def autocorr(x, lag=None):
        """
        Compute autocorrelation using FFT for every lag for the input array
        https://en.wikipedia.org/wiki/Autocorrelation#Efficient_computation
        Parameters
        ----------
        x : Numpy array
            An array containing MCMC samples
        Returns
        -------
        acorr: Numpy array same size as the input array
        """
        y = x - x.mean()
        n = len(y)
        result = fftconvolve(y, y[::-1])
        acorr = result[len(result) // 2:]
        acorr /= np.arange(n, 0, -1)
        acorr /= acorr[0]
        return acorr

    def autocov(x, lag=None):
        """Compute autocovariance estimates for every lag for the input array
        Parameters
        ----------
        x : Numpy array
            An array containing MCMC samples
        Returns
        -------
        acov: Numpy array same size as the input array
        """
        acorr = autocorr(x)
        varx = np.var(x, ddof=1) * (len(x) - 1) / len(x)
        acov = acorr * varx
        return acov

    def get_neff(x):
        """Compute the effective sample size for a 2D array
        """
        trace_value = x.T
        nchain, n_samples = trace_value.shape

        acov = np.asarray([
            autocov(trace_value[chain]) for chain in range(nchain)])

        chain_mean = trace_value.mean(axis=1)
        chain_var = acov[:, 0] * n_samples / (n_samples - 1.)
        acov_t = acov[:, 1] * n_samples / (n_samples - 1.)
        mean_var = np.mean(chain_var)
        var_plus = mean_var * (n_samples - 1.) / n_samples
        var_plus += np.var(chain_mean, ddof=1)

        rho_hat_t = np.zeros(n_samples)
        rho_hat_even = 1.
        rho_hat_t[0] = rho_hat_even
        rho_hat_odd = 1. - (mean_var - np.mean(acov_t)) / var_plus
        rho_hat_t[1] = rho_hat_odd
        # Geyer's initial positive sequence
        max_t = 1
        t = 1
        while t < (n_samples - 2) and (rho_hat_even + rho_hat_odd) >= 0.:
            rho_hat_even = 1. - (mean_var - np.mean(acov[:, t + 1])) / var_plus
            rho_hat_odd = 1. - (mean_var - np.mean(acov[:, t + 2])) / var_plus
            if (rho_hat_even + rho_hat_odd) >= 0:
                rho_hat_t[t + 1] = rho_hat_even
                rho_hat_t[t + 2] = rho_hat_odd
            max_t = t + 2
            t += 2

        # Geyer's initial monotone sequence
        t = 3
        while t <= max_t - 2:
            if (rho_hat_t[t + 1] + rho_hat_t[t + 2]) >\
                    (rho_hat_t[t - 1] + rho_hat_t[t]):
                rho_hat_t[t + 1] = (rho_hat_t[t - 1] + rho_hat_t[t]) / 2.
                rho_hat_t[t + 2] = rho_hat_t[t + 1]
            t += 2
        ess = nchain * n_samples
        ess = ess / (-1. + 2. * np.sum(rho_hat_t))
        return ess

    def generate_neff(trace_values):
        x = np.array(trace_values)
        shape = x.shape

        # Make sure to handle scalars correctly, adding extra dimensions if
        # needed. We could use np.squeeze here, but we don't want to squeeze
        # out dummy dimensions that a user inputs.
        if len(shape) == 2:
            x = np.atleast_3d(trace_values)

        # Transpose all dimensions, which makes the loop below
        # easier by moving the axes of the variable to the front instead
        # of the chain and sample axes.
        x = x.transpose()

        # Get an array the same shape as the var
        _n_eff = np.zeros(x.shape[:-2])

        # Iterate over tuples of indices of the shape of var
        for tup in np.ndindex(*list(x.shape[:-2])):
            _n_eff[tup] = get_neff(x[tup])

        if len(shape) == 2:
            return _n_eff[0]

        return np.transpose(_n_eff)

    # if not isinstance(mtrace, MultiTrace):
    #     # Return neff for non-multitrace array
    #     return generate_neff(mtrace)

    # if mtrace.nchains < 2:
    #     raise ValueError(
    #         'Calculation of effective sample size requires multiple chains '
    #         'of the same length.')

    # if varnames is None:
    #     varnames = get_default_varnames(
    #         mtrace.varnames,include_transformed=include_transformed)

    # n_eff = {}
    # for var in varnames:
    #     n_eff[var] = generate_neff(mtrace.get_values(var, combine=False))

    n_eff = generate_neff(mtrace)

    return n_eff


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

def convergenceVals(algor, ndim, varIdxs, N_conv, chains_nruns, mcmc_trace):
    """
    Convergence statistics.
    """
    # Autocorrelation time for each parameter.
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
    for p in chains_nruns.T:
        at_p = []
        for c in p:
            at_p.append(autocorr.integrated_time(c, quiet=True)[0])
        at.append(at_p)
    logger.disabled = False

    # Select the indexes of the 5 chains with the largest acorr times, and
    # the 5 chains with the smallest acorr times, for each parameter.
    if len(at[0]) >= 5:
        max_at_5c = [np.argpartition(a, -5)[-5:] for a in at]
        min_at_5c = [np.argpartition(a, 5)[:5] for a in at]
    else:
        max_at_5c, min_at_5c = [np.array([0])] * ndim, [np.array([0])] * ndim

    # Worst chain: chain with the largest acorr time.
    # max_at_c = [np.argmax(a) for a in at]

    # Mean Geweke z-scores and autocorrelation functions for all chains.
    geweke_z, emcee_acorf = [[] for _ in range(ndim)],\
        [[] for _ in range(ndim)]
    for i, p in enumerate(chains_nruns.T):
        for c in p:
            try:
                geweke_z[i].append(geweke(c))
            except ZeroDivisionError:
                geweke_z[i].append([np.nan, np.nan])
            try:
                emcee_acorf[i].append(autocorr.function_1d(c))
            except FloatingPointError:
                emcee_acorf[i].append([np.nan])
    geweke_z = np.nanmean(geweke_z, axis=1)
    emcee_acorf = np.nanmean(emcee_acorf, axis=1)

    # # PyMC3 effective sample size.
    # try:
    #     # Change shape to (nchains, nstesp, ndim)
    #     pymc3_ess = effective_n(chains_nruns.transpose(1, 0, 2))
    # except FloatingPointError:
    #     pymc3_ess = np.array([np.nan] * ndim)
    # N_steps / tau effective sample size
    mcmc_ess = mcmc_trace.shape[-1] / acorr_t

    # TODO fix this function
    # Minimum effective sample size (ESS), and multi-variable ESS.
    minESS, mESS = fminESS(ndim), multiESS(mcmc_trace.T)
    mESS_epsilon = [[], [], []]
    for alpha in [.01, .05, .1, .2, .3, .4, .5, .6, .7, .8, .9, .95]:
        mESS_epsilon[0].append(alpha)
        mESS_epsilon[1].append(fminESS(ndim, alpha=alpha, ess=minESS))
        mESS_epsilon[2].append(fminESS(ndim, alpha=alpha, ess=mESS))

    return acorr_t, max_at_5c, min_at_5c, geweke_z, emcee_acorf, mcmc_ess,\
        minESS, mESS, mESS_epsilon

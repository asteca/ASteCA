
import numpy as np
from scipy.stats import chi2
from scipy.special import gammaln
from scipy.signal import fftconvolve


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
    R"""Returns estimate of the effective sample size of a set of traces.
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
            if (
                rho_hat_t[t + 1] + rho_hat_t[t + 2]) > (rho_hat_t[t - 1] +
                    rho_hat_t[t]):
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

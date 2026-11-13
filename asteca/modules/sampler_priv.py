"""Sequential Monte Carlo core used by :py:class:`asteca.Sampler`: adaptive
tempering with persistent sampling and t-pCN MCMC moves.

This module is adapted from pocoMC v1.2.6 (https://github.com/minaskar/pocomc),
Copyright 2022-Now Minas Karamanis and contributors, distributed under the GNU General
Public License v3 (GPL-3.0).

Only the parts of pocoMC used by ASteCA are kept: bounded (uniform) priors mapped to
an unbounded space with a probit transform and a diagonal affine scaling, persistent
sampling (importance weights that combine the particles of all the iterations), the
choice of the next inverse temperature by bisection on the effective sample size
(ESS) with the dynamic adjustment of the target ESS, multinomial resampling, and
t-preconditioned Crank-Nicolson (t-pCN) MCMC moves, whose proposals use a Student-t
fitted to the particles in each iteration. The normalizing flow (and thus PyTorch),
i.e.: the preconditioning that gives Preconditioned Monte Carlo (PMC) its name, the
evidence estimation with the flow, the blobs, the checkpointing and the
multiprocessing pools are not included. Hence this is not PMC but the Sequential
Monte Carlo sampler that pocoMC runs with ``precondition=False``: with these settings
the results are those of pocoMC for the same random seed.

Changes with respect to pocoMC's code, which do not change the results:

- :py:func:`trim_weights` (discards the samples with negligible importance weights,
  in each iteration and for the posterior) obtains the ESS of all the candidate
  thresholds at once from the cumulative sums of the sorted weights, instead of
  computing them one at a time with a ``np.percentile()`` call over all the weights
  for each of the (up to 1000) candidates. The selected threshold is then checked
  with pocoMC's computation, so the kept samples and weights are identical. In
  pocoMC this function dominated the run time of the sampler itself (excluding the
  likelihood).

References:

- Karamanis et al. (2022), MNRAS, 516, 1644: "Accelerating astronomical and
  cosmological inference with preconditioned Monte Carlo"
  https://academic.oup.com/mnras/article/516/2/1644/6670808
- Karamanis et al. (2022), JOSS, 7, 79: "pocoMC: A Python package for
  accelerated Bayesian inference in astronomy and cosmology"
  https://joss.theoj.org/papers/10.21105/joss.04634
- Karamanis & Seljak (2025), Stat. Comput., 35, 144: "Persistent sampling: enhancing
  the efficiency of sequential Monte Carlo"
  https://link.springer.com/article/10.1007/s11222-025-10682-y
"""

from collections.abc import Callable

import numpy as np
from scipy import optimize, special

SQRTEPS = float(np.sqrt(np.finfo(np.float64).eps))


class PersistentSMC:
    """Sequential Monte Carlo sampler with adaptive tempering, persistent sampling,
    and t-pCN moves (pocoMC without the normalizing flow preconditioning).

    :param log_like: Vectorized log-likelihood, takes an array of shape
        ``(n_particles, n_dim)`` and returns an array of shape ``(n_particles,)``
    :type log_like: Callable[[np.ndarray], np.ndarray]
    :param bounds: Array of shape ``(n_dim, 2)`` with the ``[min, max]`` limits of the
        uniform priors
    :type bounds: np.ndarray
    :param n_effective: Target ESS used to select the next inverse temperature
    :type n_effective: int
    :param n_active: Number of particles moved with MCMC in each iteration
    :type n_active: int
    :param rng: Random number generator. A legacy ``RandomState`` is used so that the
        stream of random numbers is the same as pocoMC's for the same seed
    :type rng: np.random.RandomState
    :param dynamic: Adjust the target ESS based on the number of unique particles
    :type dynamic: bool
    :param callback: Function called at the end of each iteration with the dictionary
        of the current particles (used to print the progress)
    :type callback: Callable[[dict], None] | None

    :raises ValueError: If the prior bounds are not finite or not ``min < max``
    """

    def __init__(
        self,
        log_like: Callable[[np.ndarray], np.ndarray],
        bounds: np.ndarray,
        n_effective: int,
        n_active: int,
        rng: np.random.RandomState,
        dynamic: bool = True,
        callback: Callable[[dict], None] | None = None,
    ) -> None:
        self.log_like = log_like
        self.rng = rng
        self.callback = callback

        # Uniform priors
        self.bounds = np.asarray(bounds, dtype=float)
        self.low, self.high = self.bounds.T
        if not (np.all(np.isfinite(self.bounds)) and np.all(self.high > self.low)):
            raise ValueError("Prior bounds must be finite and with min < max")
        self.n_dim = len(self.bounds)

        self.n_active = int(n_active)
        self.n_effective = int(n_effective)
        # Number of MCMC steps after the logP plateau, and maximum number of steps
        self.n_steps = int(self.n_dim // 2)
        self.n_max_steps = 10 * self.n_steps
        # Number of prior samples drawn in the warm-up
        self.n_prior = int(
            2 * max(self.n_effective // self.n_active, 1) * self.n_active
        )

        # Dynamic ESS: keep the original ratio of unique active to effective particles
        self.dynamic = dynamic
        self.dynamic_ratio = (
            unique_sample_size(np.ones(self.n_effective), k=self.n_active)
            / self.n_active
        )

        # Proposal scale of the t-pCN moves
        self.proposal_scale = 2.38 / self.n_dim**0.5

        # Probit + affine reparameterization (fitted on the prior samples)
        self.scaler = Reparameterize(self.low, self.high)
        # Multivariate Student-t fitted to the particles in the unbounded space
        self.geometry = Geometry(rng)

        self.particles = Particles()
        self.current = {}
        self.t = 0
        self.calls = 0
        self.warmup = True
        self.logz: float

    # ---------------------------------------------------------------------------
    # Prior
    # ---------------------------------------------------------------------------
    def log_prior(self, x: np.ndarray) -> np.ndarray:
        """Log-density of the uniform priors (``-inf`` outside the limits)"""
        logp = np.zeros(len(x))
        for i in range(self.n_dim):
            inside = (x[:, i] >= self.low[i]) & (x[:, i] <= self.high[i])
            logp += np.where(inside, -np.log(self.high[i] - self.low[i]), -np.inf)
        return logp

    def sample_prior(self, size: int) -> np.ndarray:
        """Samples of the uniform priors, of shape ``(size, n_dim)``"""
        samples = [
            self.rng.uniform(0.0, 1.0, size) * (hi - lo) + lo
            for lo, hi in zip(self.low, self.high)
        ]
        return np.transpose(samples)

    # ---------------------------------------------------------------------------
    # Main loop
    # ---------------------------------------------------------------------------
    def run(self, n_total: int) -> None:
        """Run the sampler until the inverse temperature reaches ``1`` and the
        ESS of the persistent (importance weighted) samples is at least ``n_total``.
        Calling it again with a larger ``n_total`` continues the run.
        """
        self.n_total = int(n_total)

        if self.warmup:
            self._warmup()

        while self._not_termination():
            self._reweight()
            self._train()
            self._resample()
            self._mutate()
            self.particles.update(self.current)
            if self.callback is not None:
                self.callback(self.current)

        _, self.logz = self.particles.compute_logw_and_logz(1.0)

    def _warmup(self) -> None:
        """Draw the initial particles from the prior, in batches of ``n_active``"""
        prior_samples = self.sample_prior(self.n_prior)
        self.scaler.fit(prior_samples)

        for i in range(self.n_prior // self.n_active):
            x = prior_samples[i * self.n_active : (i + 1) * self.n_active]
            u = self.scaler.forward(x)
            logdetj = self.scaler.inverse(u)[1]
            logp = self.log_prior(x)
            logl = self.log_like(x)
            self.calls += self.n_active

            # Replace prior particles with infinite likelihoods by finite ones
            inf_msk = np.isinf(logl)
            if np.any(inf_msk):
                finite_idx = np.flatnonzero(~inf_msk)
                if len(finite_idx) == 0:
                    raise ValueError(
                        "All the prior samples in a batch returned an infinite "
                        + "log-likelihood"
                    )
                infinite_idx = np.flatnonzero(inf_msk)
                idx = self.rng.choice(finite_idx, size=len(infinite_idx), replace=True)
                for arr in (x, u, logdetj, logp, logl):
                    arr[infinite_idx] = arr[idx]

            self.current = {
                "u": u,
                "x": x,
                "logl": logl,
                "logp": logp,
                "logdetj": logdetj,
                "iter": self.t,
                "calls": self.calls,
                "steps": 1,
                "efficiency": 1.0,
                "ess": self.n_effective,
                "accept": 1.0,
                "beta": 0.0,
                "logz": 0.0,
            }
            self.particles.update(self.current)
            if self.callback is not None:
                self.callback(self.current)
            self.t += 1

        self.warmup = False

    def _not_termination(self) -> bool:
        """True while the termination criterion is not satisfied"""
        logw, _ = self.particles.compute_logw_and_logz(1.0)
        ess = effective_sample_size(np.exp(logw - np.max(logw)))
        return 1.0 - self.current["beta"] >= 1e-4 or ess < self.n_total

    def _reweight(self) -> None:
        """Select the next inverse temperature (bisection on the ESS of the
        persistent weights), and get the weighted particles of all iterations"""
        self.t += 1
        beta_prev = float(self.particles.get("beta", index=-1))
        beta_max = 1.0
        beta_min = beta_prev

        def weights_and_ess(beta):
            logw, _ = self.particles.compute_logw_and_logz(beta)
            weights = np.exp(logw - np.max(logw))
            return weights, effective_sample_size(weights)

        _, ess_prev = weights_and_ess(beta_prev)
        _, ess_max = weights_and_ess(beta_max)

        if ess_prev <= self.n_effective:
            beta, ess_est = beta_prev, ess_prev
            logz = self.particles.get("logz", index=-1)
        elif ess_max >= self.n_effective:
            beta, ess_est = beta_max, ess_max
            _, logz = self.particles.compute_logw_and_logz(beta)
        else:
            while True:
                beta = (beta_max + beta_min) * 0.5
                _, ess_est = weights_and_ess(beta)
                if (
                    np.abs(ess_est - self.n_effective) < 0.01 * self.n_effective
                    or beta == 1.0
                ):
                    _, logz = self.particles.compute_logw_and_logz(beta)
                    break
                elif ess_est < self.n_effective:
                    beta_max = beta
                else:
                    beta_min = beta

        logw, _ = self.particles.compute_logw_and_logz(beta)
        weights = np.exp(logw - np.max(logw))
        weights /= np.sum(weights)

        if self.dynamic:
            # Adjust the target ESS based on the expected number of unique particles
            n_unique_active = unique_sample_size(weights, k=self.n_active)
            if n_unique_active < self.n_active * (0.95 * self.dynamic_ratio):
                self.n_effective = int(
                    self.n_active / n_unique_active * self.n_effective
                )
            elif n_unique_active > self.n_active * min(1.05 * self.dynamic_ratio, 1.0):
                self.n_effective = int(
                    n_unique_active / self.n_active * self.n_effective
                )

        idx, weights = trim_weights(np.arange(len(weights)), weights)
        for key in ("u", "x", "logdetj", "logl", "logp"):
            self.current[key] = self.particles.get(key, flat=True)[idx]
        self.current["logz"] = logz
        self.current["beta"] = beta
        self.current["weights"] = weights
        self.current["ess"] = ess_est
        self.current["iter"] = self.t

    def _train(self) -> None:
        """Fit the Student-t distribution used by the t-pCN proposals"""
        self.geometry.fit(self.current["u"], weights=self.current["weights"])

    def _resample(self) -> None:
        """Multinomial resampling of ``n_active`` particles"""
        weights = self.current["weights"]
        idx = self.rng.choice(
            np.arange(len(weights)), size=self.n_active, replace=True, p=weights
        )
        for key in ("u", "x", "logdetj", "logl", "logp"):
            self.current[key] = self.current[key][idx]

    def _mutate(self) -> None:
        """Move the particles with t-preconditioned Crank-Nicolson MCMC"""
        results = self._tpcn()
        for key in ("u", "x", "logdetj", "logl", "logp"):
            self.current[key] = results[key]
        self.current["efficiency"] = results["efficiency"] / (2.38 / self.n_dim**0.5)
        self.current["steps"] = results["steps"]
        self.current["accept"] = results["accept"]
        self.calls += results["calls"]
        self.current["calls"] = self.calls
        self.proposal_scale = results["proposal_scale"]

    def _tpcn(self) -> dict:
        """t-preconditioned Crank-Nicolson moves in the unbounded space"""
        cur = self.current
        u, x = cur["u"].copy(), cur["x"].copy()
        logdetj, logl, logp = (
            cur["logdetj"].copy(),
            cur["logl"].copy(),
            cur["logp"].copy(),
        )
        beta = cur["beta"]

        n_calls = 0
        sigma = min(self.proposal_scale, 0.99)
        n_walkers, n_dim = x.shape
        mu, cov, nu = self.geometry.t_mean, self.geometry.t_cov, self.geometry.t_nu
        inv_cov = np.linalg.inv(cov)
        chol_cov = np.linalg.cholesky(cov)

        logp2_val = np.mean(logl + logp)
        cnt, i = 0, 0
        while True:
            i += 1
            diff = u - mu
            # Mahalanobis distances of the current points
            maha = np.einsum("ij,jk,ik->i", diff, inv_cov, diff)
            # Scale mixture variables (one per walker)
            s = 1.0 / self.rng.gamma((n_dim + nu) / 2, 2.0 / (nu + maha))
            # Propose new points in the unbounded space
            z = self.rng.randn(n_walkers, n_dim) @ chol_cov.T
            u_prime = (
                mu + (1.0 - sigma**2) ** 0.5 * diff + sigma * np.sqrt(s)[:, None] * z
            )

            # Transform to the parameter space
            x_prime, logdetj_prime = self.scaler.inverse(u_prime)

            finite_msk = np.isfinite(logdetj_prime) & np.isfinite(x_prime).all(axis=1)
            logp_prime = np.full(n_walkers, -np.inf)
            logp_prime[finite_msk] = self.log_prior(x_prime[finite_msk])
            finite_msk &= np.isfinite(logp_prime)

            logl_prime = np.full(n_walkers, -np.inf)
            if finite_msk.any():
                logl_prime[finite_msk] = self.log_like(x_prime[finite_msk])
            n_calls += int(np.sum(finite_msk))

            # Metropolis-Hastings acceptance
            diff_prime = u_prime - mu
            maha_prime = np.einsum("ij,jk,ik->i", diff_prime, inv_cov, diff_prime)
            A = -(n_dim + nu) / 2 * np.log(1 + maha_prime / nu)
            B = -(n_dim + nu) / 2 * np.log(1 + maha / nu)
            with np.errstate(invalid="ignore", over="ignore"):
                alpha = np.minimum(
                    1.0,
                    np.exp(
                        logl_prime * beta
                        - logl * beta
                        + logp_prime
                        - logp
                        + logdetj_prime
                        - logdetj
                        - A
                        + B
                    ),
                )
            alpha[np.isnan(alpha)] = 0.0

            msk = self.rng.rand(n_walkers) < alpha
            u[msk] = u_prime[msk]
            x[msk] = x_prime[msk]
            logdetj[msk] = logdetj_prime[msk]
            logl[msk] = logl_prime[msk]
            logp[msk] = logp_prime[msk]

            # Adapt the scale towards an acceptance rate of 0.234 (diminishing
            # adaptation)
            sigma = np.abs(
                np.minimum(
                    sigma + 1 / (i + 1) ** 0.75 * (np.mean(alpha) - 0.234),
                    np.minimum(2.38 / n_dim**0.5, 0.99),
                )
            )

            # Stop once the mean logP has not increased for long enough
            logp2_val_new = np.mean(logl + logp)
            if logp2_val_new > logp2_val:
                cnt = 0
                logp2_val = logp2_val_new
            else:
                cnt += 1
                if cnt >= self.n_steps * ((2.38 / n_dim**0.5) / sigma) ** 2.0:
                    break
            if i >= self.n_max_steps:
                break

        return {
            "u": u,
            "x": x,
            "logdetj": logdetj,
            "logl": logl,
            "logp": logp,
            "efficiency": sigma,
            "accept": np.mean(alpha),
            "steps": i,
            "calls": n_calls,
            "proposal_scale": sigma,
        }

    # ---------------------------------------------------------------------------
    # Results
    # ---------------------------------------------------------------------------
    def posterior(
        self, resample: bool = False
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Posterior samples from all the iterations (persistent sampling).

        :param resample: If ``True``, return equally weighted samples (multinomial
            resampling of the trimmed importance weighted samples)
        :type resample: bool

        :return: ``(samples, weights, logl)``; ``weights`` are all equal to
            ``1/N`` if ``resample=True``
        :rtype: tuple[np.ndarray, np.ndarray, np.ndarray]
        """
        samples = self.particles.get("x", flat=True)
        logl = self.particles.get("logl", flat=True)
        logw, _ = self.particles.compute_logw_and_logz(1.0)
        idx, weights = trim_weights(np.arange(len(samples)), np.exp(logw))
        samples, logl = samples[idx], logl[idx]

        if resample:
            idx_r = self.rng.choice(
                np.arange(len(weights)), size=len(samples), replace=True, p=weights
            )
            N = len(samples)
            return samples[idx_r], np.full(N, 1.0 / N), logl[idx_r]

        return samples, weights, logl


class Particles:
    """Store the particles of all the iterations (persistent sampling)"""

    keys = (
        "u",
        "x",
        "logdetj",
        "logl",
        "logp",
        "iter",
        "logz",
        "calls",
        "steps",
        "efficiency",
        "ess",
        "accept",
        "beta",
    )

    def __init__(self) -> None:
        self.past = {k: [] for k in self.keys}

    def update(self, data: dict) -> None:
        for key in self.keys:
            if key in data:
                self.past[key].append(data[key])

    def get(self, key: str, index: int | None = None, flat: bool = False):
        if index is not None:
            return self.past[key][index]
        if flat:
            return np.concatenate(self.past[key])
        return np.asarray(self.past[key])

    def compute_logw_and_logz(
        self, beta_final: float = 1.0
    ) -> tuple[np.ndarray, float]:
        """Normalized log importance weights of all the particles at the inverse
        temperature ``beta_final``, and the log-evidence estimate.

        The particles of iteration ``i`` are treated as samples of the mixture of all
        the tempered posteriors ``L^beta_i * prior / Z_i`` (persistent sampling).
        """
        logz = self.get("logz")
        logl = self.get("logl")
        beta = self.get("beta")

        A = logl * beta_final
        b = logl[None, :, :] * beta[:, None, None] - logz[:, None, None]
        B = np.logaddexp.reduce(b, axis=0) - np.log(len(beta))
        logw = np.concatenate(A - B)
        logz_new = np.logaddexp.reduce(logw) - np.log(len(logw))
        logw -= np.logaddexp.reduce(logw)
        return logw, float(logz_new)


class Reparameterize:
    """Map the bounded parameters to an unbounded space (probit transform),
    followed by a diagonal affine scaling learned from the prior samples"""

    def __init__(self, low: np.ndarray, high: np.ndarray) -> None:
        self.low = low
        self.high = high
        self.width = high - low
        self.mu: np.ndarray
        self.sigma: np.ndarray

    def fit(self, x: np.ndarray) -> None:
        u = self._forward(x)
        self.mu = np.mean(u, axis=0)
        self.sigma = np.std(u, axis=0)

    def forward(self, x: np.ndarray) -> np.ndarray:
        return (self._forward(x) - self.mu) / self.sigma

    def inverse(self, u: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Return the parameters and the log-determinant of the Jacobian"""
        v = self.mu + self.sigma * u
        log_det_J = np.sum(np.log(self.sigma)) * np.ones(len(u))
        p = (special.erf(v / np.sqrt(2.0)) + 1.0) / 2.0
        x = p * self.width + self.low
        J = np.log(self.width) + (-(v**2.0) / 2.0) - np.log(np.sqrt(2.0 * np.pi))
        return x, log_det_J + np.sum(J, axis=1)

    def _forward(self, x: np.ndarray) -> np.ndarray:
        if np.any(x < self.low) or np.any(x > self.high):
            raise ValueError("Input outside of the prior limits")
        p = np.clip((x - self.low) / self.width, 1e-13, 1.0 - 1e-13)
        return np.sqrt(2.0) * special.erfinv(2.0 * p - 1.0)


class Geometry:
    """Multivariate Student-t fitted to the (resampled) weighted particles"""

    def __init__(self, rng: np.random.RandomState) -> None:
        self.rng = rng
        self.t_mean: np.ndarray
        self.t_cov: np.ndarray
        self.t_nu: float

    def fit(self, theta: np.ndarray, weights: np.ndarray) -> None:
        idx = systematic_resample(len(theta), weights, self.rng)
        self.t_mean, self.t_cov, self.t_nu = fit_mvstud(theta[idx])
        if not np.isfinite(self.t_nu):
            self.t_nu = 1e6


def fit_mvstud(
    data: np.ndarray, tolerance: float = 1e-6, max_iter: int = 100
) -> tuple[np.ndarray, np.ndarray, float]:
    """Fit a multivariate Student-t distribution to ``data`` (shape ``(n, dim)``)
    with the EM algorithm. Returns the mean, scale matrix, and degrees of freedom.
    """

    def opt_nu(delta_iobs):
        def func0(nu):
            w_iobs = (nu + dim) / (nu + delta_iobs)
            return (
                -special.psi(nu / 2)
                + np.log(nu / 2)
                + np.sum(np.log(w_iobs)) / n
                - np.sum(w_iobs) / n
                + 1
                + special.psi((nu + dim) / 2)
                - np.log((nu + dim) / 2)
            )

        if func0(1e300) >= 0:
            return np.inf
        return optimize.bisect(func0, 1e-300, 1e300)

    data = data.T
    dim, n = data.shape
    mu = np.array([np.median(data, 1)]).T
    Sigma = np.cov(data) * (n - 1) / n + (1 / n) * np.diag(np.var(data, axis=1))
    nu = 20

    last_nu = 0
    i = 0
    while np.abs(last_nu - nu) > tolerance and i < max_iter:
        i += 1
        diffs = data - mu
        delta_iobs = np.sum(diffs * np.linalg.solve(Sigma, diffs), 0)

        last_nu = nu
        nu = opt_nu(delta_iobs)
        if nu == np.inf:
            return mu.T[0], Sigma, nu

        w_iobs = (nu + dim) / (nu + delta_iobs)
        Sigma = np.dot(w_iobs * diffs, diffs.T) / n

        mu = np.sum(w_iobs * data, 1) / sum(w_iobs)
        mu = np.array([mu]).T

    return mu.T[0], Sigma, nu


def trim_weights(
    samples: np.ndarray, weights: np.ndarray, ess: float = 0.99, bins: int = 1000
) -> tuple[np.ndarray, np.ndarray]:
    """Discard the samples with the smallest weights while keeping a fraction
    ``ess`` of the ESS. Returns the kept samples and their normalized weights.

    The threshold is the largest of the ``bins`` percentiles of the weights (from the
    99th down to the 0th) whose kept samples (``weights >= threshold``) retain the
    fraction ``ess`` of the ESS, as in pocoMC. pocoMC tests the percentiles one at a
    time, each with its own ``np.percentile()`` call over all the weights, which
    dominated the run time of the sampler. Here the ESS of all the candidates is
    obtained at once from the cumulative sums of the sorted weights, and the
    selected percentile is then checked with pocoMC's exact computation, so that the
    result is identical. Since the samples are added from the largest weight down,
    and each added weight is not larger than the mean (weighted by the weights) of
    those already kept, the ESS can only grow as the threshold decreases: the first
    percentile (from the top) that passes the test is the boundary.
    """
    weights = weights / np.sum(weights)
    ess_total = 1.0 / np.sum(weights**2.0)
    percentiles = np.linspace(0, 99, bins)

    def trimmed(i: int) -> tuple[np.ndarray, np.ndarray, bool]:
        """pocoMC's computation for the i-th percentile"""
        threshold = np.percentile(weights, percentiles[i])
        mask = weights >= threshold
        weights_trimmed = weights[mask]
        weights_trimmed /= np.sum(weights_trimmed)
        ess_trimmed = 1.0 / np.sum(weights_trimmed**2.0)
        return mask, weights_trimmed, bool(ess_trimmed / ess_total >= ess)

    # ESS ratio of every candidate threshold, from the cumulative sums of the weights
    # sorted in decreasing order
    w_asc = np.sort(weights)
    w_desc = w_asc[::-1]
    c1 = np.cumsum(w_desc)
    c2 = np.cumsum(w_desc**2.0)
    # Number of weights >= each threshold (always >= 1)
    n_kept = len(weights) - np.searchsorted(
        w_asc, np.percentile(weights, percentiles), side="left"
    )
    ratio = (c1[n_kept - 1] ** 2 / c2[n_kept - 1]) / ess_total
    # Largest percentile that keeps the required ESS (the 0th keeps all the samples)
    passing = np.flatnonzero(ratio >= ess)
    i = int(passing[-1]) if passing.size else 0

    # Check with the exact computation, moving down (or up) if the cumulative sums
    # put the boundary one percentile off due to rounding
    mask, weights_trimmed, ok = trimmed(i)
    while not ok and i > 0:
        i -= 1
        mask, weights_trimmed, ok = trimmed(i)
    while i < bins - 1:
        mask_up, weights_up, ok_up = trimmed(i + 1)
        if not ok_up:
            break
        i, mask, weights_trimmed = i + 1, mask_up, weights_up

    return samples[mask], weights_trimmed


def effective_sample_size(weights: np.ndarray) -> float:
    """Kish effective sample size of (unnormalized) weights"""
    weights = weights / np.sum(weights)
    return float(1.0 / np.sum(weights**2.0))


def unique_sample_size(weights: np.ndarray, k: int | None = None) -> float:
    """Expected number of unique samples when ``k`` samples are drawn"""
    if k is None:
        k = len(weights)
    weights = weights / np.sum(weights)
    return float(np.sum(1.0 - (1.0 - weights) ** k))


def systematic_resample(
    size: int, weights: np.ndarray, rng: np.random.RandomState
) -> np.ndarray:
    """Indexes of a systematic resampling of ``size`` samples"""
    if abs(np.sum(weights) - 1.0) > SQRTEPS:
        weights = np.asarray(weights) / np.sum(weights)
    positions = (rng.random_sample() + np.arange(size)) / size
    cumulative_sum = np.cumsum(weights)
    idx = np.searchsorted(cumulative_sum, positions, side="left")
    return np.minimum(idx, len(weights) - 1)

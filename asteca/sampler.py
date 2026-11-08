import numpy as np

from .likelihood import Likelihood
from .synthetic import Synthetic


class Sampler:
    """Define a :py:class:`Sampler` object.

    This object estimates the fundamental parameters of an observed cluster with
    Preconditioned Monte Carlo (PMC), a Sequential Monte Carlo method that moves a set
    of particles from the prior to the posterior by slowly increasing the weight given
    to the likelihood. The implementation is adapted from
    `pocoMC <https://pocomc.readthedocs.io/>`__ (Karamanis et al.
    `2022a <https://ui.adsabs.harvard.edu/abs/2022MNRAS.516.1644K>`__,
    `2022b <https://arxiv.org/abs/2207.05660>`__), without the normalizing flow
    preconditioning, and with persistent sampling
    (`Karamanis & Seljak 2025 <https://arxiv.org/abs/2407.20722>`__).

    The synthetic clusters are generated with the :py:meth:`asteca.Synthetic.generate`
    method of the (calibrated) ``synthcl`` object, and compared to the observed
    cluster with the :py:meth:`asteca.Likelihood.log_lkl` method of the
    ``likelihood`` object. The parameters that are not included in ``priors`` keep
    the default values of the ``synthcl`` object.

    :param synthcl: :py:class:`asteca.Synthetic` object, calibrated with the observed
        cluster via :py:meth:`asteca.Synthetic.calibrate`
    :type synthcl: Synthetic
    :param likelihood: :py:class:`asteca.Likelihood` object defined with the observed
        cluster
    :type likelihood: Likelihood
    :param priors: Dictionary with the ``(min, max)`` limits of the uniform priors of
        the fitted parameters, e.g.: ``{"loga": (7.0, 9.5), "dm": (8.0, 10.5)}``. The
        keys must be parameters of the ``synthcl`` object (``met, loga, alpha, beta,
        Av, DR, Rv, dm``)
    :type priors: dict[str, tuple[float, float]]
    :param temperature: Temperature passed to :py:meth:`asteca.Likelihood.log_lkl`.
        Values larger than ``1`` widen the posterior (by a factor of about
        ``sqrt(temperature)``)
    :type temperature: float
    :param n_effective: Target effective sample size used to select the next
        inverse temperature in each iteration. It is adjusted during the run based
        on the number of unique particles
    :type n_effective: int
    :param n_active: Number of particles moved with MCMC in each iteration
    :type n_active: int
    :param seed: Random seed. If ``None`` a random integer will be generated and used
    :type seed: int | None
    :param verbose: Verbose level. A value of ``0`` hides all output
    :type verbose: int

    :raises ValueError: If any of the attributes is not recognized as a valid option
    """

    def __init__(
        self,
        synthcl: Synthetic,
        likelihood: Likelihood,
        priors: dict[str, tuple[float, float]],
        temperature: float = 1.0,
        n_effective: int = 256,
        n_active: int = 128,
        seed: int | None = None,
        verbose: int = 1,
    ) -> None:
        self.synthcl = synthcl
        self.likelihood = likelihood
        self.temperature = float(temperature)
        self.n_effective = int(n_effective)
        self.n_active = int(n_active)
        self.seed = seed
        self.verbose = verbose

        self.names = list(priors)
        self.bounds = np.array([priors[k] for k in self.names], dtype=float)
        self._check_args()

        if self.seed is None:
            self.seed = int(np.random.randint(0, 100_000))
        # A legacy RandomState is used so that the results match those of pocoMC
        # (with precondition=False) for the same seed
        self.rng = np.random.RandomState(self.seed)

        # Best solution found during the run (largest log-likelihood)
        self.best_lnl = -np.inf
        self.best_params = {}

        from .modules.sampler_priv import PMC

        self._pmc = PMC(
            self._log_like,
            self.bounds,
            self.n_effective,
            self.n_active,
            self.rng,
            callback=self._progress if self.verbose > 1 else None,
        )

        self._vp("\nSampler object generated")
        self._vp(f"Fitted params  : {self.names}", 1)
        self._vp(f"Temperature    : {self.temperature:g}", 1)
        self._vp(f"n_effective    : {self.n_effective}", 1)
        self._vp(f"n_active       : {self.n_active}", 1)
        self._vp(f"Random seed    : {self.seed}", 1)

    def _vp(self, mssg: str, level: int = 0) -> None:
        """Verbose print method"""
        if self.verbose > level:
            print(mssg)

    def _check_args(self) -> None:
        """Validate the input arguments"""
        if len(self.names) == 0:
            raise ValueError("At least one parameter must be included in 'priors'")
        for k, (pmin, pmax) in zip(self.names, self.bounds):
            if k not in self.synthcl.def_params:
                raise ValueError(
                    f"Parameter '{k}' not recognized. Should be one of "
                    + f"{list(self.synthcl.def_params)}"
                )
            if not (np.isfinite(pmin) and np.isfinite(pmax) and pmin < pmax):
                raise ValueError(f"Prior of '{k}' must be finite and with min < max")

        isochs = self.synthcl.isochs
        lims = {
            "met": (isochs.zmin, isochs.zmax),
            "loga": (isochs.amin, isochs.amax),
        }
        for k, (pmin, pmax) in zip(self.names, self.bounds):
            if k in lims:
                lmin, lmax = lims[k]
                if pmin < lmin or pmax > lmax:
                    raise ValueError(
                        f"Prior of '{k}' {pmin, pmax} is outside of the range of "
                        + f"the loaded isochrones [{lmin}, {lmax}]"
                    )
            elif k == "Rv" and pmin <= 0.0:
                raise ValueError("Prior of 'Rv' must be larger than 0")
            elif k != "Rv" and pmin < 0.0:
                raise ValueError(f"Prior of '{k}' must not be negative")

        if self.temperature <= 0.0:
            raise ValueError("'temperature' must be larger than 0")
        if self.n_active < 2 or self.n_effective < self.n_active:
            raise ValueError("Must have 'n_effective' >= 'n_active' >= 2")

    def _log_like(self, params_arr: np.ndarray) -> np.ndarray:
        """Vectorized log-likelihood, ``params_arr`` has shape
        ``(n_particles, n_dim)``"""
        lnl = np.array(
            [
                self.likelihood.log_lkl(
                    self.synthcl.generate(dict(zip(self.names, p))), self.temperature
                )
                for p in params_arr
            ],
            dtype=float,
        )
        lnl[~np.isfinite(lnl)] = -np.inf

        i = int(np.argmax(lnl))
        if lnl[i] > self.best_lnl:
            self.best_lnl = float(lnl[i])
            self.best_params = dict(zip(self.names, map(float, params_arr[i])))
        return lnl

    def _progress(self, particles: dict) -> None:
        """Print the progress of the run"""
        print(
            f"iter {particles['iter']:>3d} | beta {particles['beta']:.4f} | "
            + f"ESS {particles['ess']:>6.0f} | calls {particles['calls']:>7d} | "
            + f"accept {particles['accept']:.3f} | steps {particles['steps']:>3d}"
        )

    def run(self, n_total: int = 512) -> None:
        """Run the sampler.

        The run ends when the inverse temperature of the likelihood reaches ``1``
        (i.e.: the particles sample the posterior) and the effective sample size of
        the importance weighted samples of all the iterations is at least
        ``n_total``. Calling this method again with a larger ``n_total`` continues the
        run.

        :param n_total: Effective sample size of the posterior samples required to
            stop the run
        :type n_total: int
        """
        self._vp("\nRunning sampler...", 1)
        self._pmc.run(n_total)

        self.n_calls = int(self._pmc.calls)
        self.n_iter = int(self._pmc.t)
        self.logz = float(self._pmc.logz)
        _, weights, _ = self._pmc.posterior()
        self.ess = float(weights.sum() ** 2 / (weights**2).sum())
        # Distance of the best solution (see Likelihood.log_lkl)
        self.best_dist = float(
            -self.best_lnl * self.temperature / self.likelihood.lkl_range
        )

        self._vp(f"Iterations     : {self.n_iter}", 1)
        self._vp(f"Likelihood calls: {self.n_calls}", 1)
        self._vp(f"ESS            : {self.ess:.0f}", 1)
        self._vp(f"Best distance  : {self.best_dist:.3f}", 1)
        self._vp("Sampler run completed")

    def posterior(self, resample: bool = False) -> tuple[np.ndarray, np.ndarray]:
        """Posterior samples, after :py:meth:`run`.

        The samples of all the iterations are combined with importance weights
        (persistent sampling), and those with negligible weights are discarded.

        :param resample: If ``True``, return equally weighted samples (resampled
            with replacement according to their importance weights)
        :type resample: bool

        :return: ``(samples, weights)``, where ``samples`` has shape
            ``(N, n_params)`` with the columns in the order of the ``priors``
            dictionary (stored in the ``names`` attribute), and ``weights`` are the
            normalized importance weights (all equal to ``1/N`` if ``resample=True``)
        :rtype: tuple[np.ndarray, np.ndarray]

        :raises ValueError: If the sampler has not been run
        """
        if self._pmc.warmup:
            raise ValueError("The sampler must be run first, see the 'run()' method")
        samples, weights, _ = self._pmc.posterior(resample)
        return samples, weights

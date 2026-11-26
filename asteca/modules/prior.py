"""Priors of the parameters fitted by :py:class:`asteca.Sampler`: uniform ``(min,
max)`` limits, or ``(dist, min, max)`` frozen ``scipy.stats`` distributions truncated
to ``[min, max]`` and renormalized. Used by the sampling engine
(:py:class:`asteca.modules.smc.PersistentSMC`) for the log-density of the particles,
the draws of the initial ones, and the (finite) limits of its probit transform."""

from collections.abc import Callable

import numpy as np


class Prior:
    """Independent priors of the parameters, one entry per parameter:

    - ``(min, max)``: uniform prior in ``[min, max]``
    - ``(dist, min, max)``: frozen ``scipy.stats`` continuous distribution truncated
      to ``[min, max]`` (intersected with the support of ``dist``) and renormalized to
      the probability mass in that range, so that the evidence is correct

    The uniform entries reproduce the original (uniform only) code exactly: same
    log-density and same random draws, so the results are unchanged for the same
    seed. The truncated distributions are drawn with the inverse of their CDF (one
    array of uniform values per parameter, as for the uniform priors).

    The log-density of a distribution is evaluated in every MCMC step, and the
    public ``dist.logpdf()`` spends ~0.1 ms per call validating its arguments
    (independently of the number of particles). Its ``loc, scale`` and shape
    parameters are hence unpacked once, and the distribution's ``_logpdf()`` method is
    called directly (6-9x faster). This uses a private ``scipy`` API: it is checked
    against ``dist.logpdf()`` once, and the public method is used if it fails or
    gives different values.

    :param names: Names of the parameters (used in the error messages)
    :type names: list[str]
    :param entries: Prior of each parameter, ``(min, max)`` or ``(dist, min, max)``
    :type entries: list

    :raises ValueError: If an entry is not one of the two forms, if ``dist`` is not a
        frozen continuous ``scipy.stats`` distribution, if its support does not
        overlap ``[min, max]``, or if it has no probability mass in that range
    """

    def __init__(self, names: list[str], entries: list) -> None:
        self.names = list(names)
        # (kind, low, high, log_norm, logpdf, cdf_lo, cdf_hi, upper_tail, dist)
        self._terms: list[tuple] = []
        for name, entry in zip(self.names, entries):
            self._terms.append(self._parse(name, entry))
        self.bounds = np.array([[t[1], t[2]] for t in self._terms], dtype=float)

    @staticmethod
    def _parse(name: str, entry) -> tuple:
        """Validate one entry and precompute its constants"""
        msg = (
            f"Prior of '{name}' must be '(min, max)' or '(dist, min, max)', with "
            + "'dist' a frozen scipy.stats continuous distribution"
        )
        if hasattr(entry, "logpdf"):
            raise ValueError(msg + ", got a distribution without its '(min, max)'")
        try:
            n_el = len(entry)
        except TypeError:
            raise ValueError(msg) from None

        if n_el == 2 and not any(hasattr(e, "logpdf") for e in entry):
            try:
                lo, hi = float(entry[0]), float(entry[1])
            except (TypeError, ValueError):
                raise ValueError(msg) from None
            # Infinite or inverted limits are reported by asteca.Sampler
            with np.errstate(all="ignore"):
                log_norm = -np.log(hi - lo) if hi > lo else -np.inf
            return ("uniform", lo, hi, log_norm, None, 0.0, 1.0, False, None)

        if n_el != 3:
            raise ValueError(msg)
        dist = entry[0]
        if not all(hasattr(dist, a) for a in ("logpdf", "cdf", "ppf", "support")):
            raise ValueError(msg)
        try:
            lo, hi = float(entry[1]), float(entry[2])
        except (TypeError, ValueError):
            raise ValueError(msg) from None
        s_lo, s_hi = (float(v) for v in dist.support())
        lo, hi = max(lo, s_lo), min(hi, s_hi)
        if not lo < hi:
            raise ValueError(
                f"Prior of '{name}': the support of the distribution does not overlap "
                + "its '(min, max)' range"
            )
        # Probability mass in [lo, hi]. If the range is in the upper tail of the
        # distribution the survival function is used, to keep the precision of the
        # draws with the inverse CDF
        c_lo, c_hi = float(dist.cdf(lo)), float(dist.cdf(hi))
        upper_tail = c_lo > 0.5
        if upper_tail:
            c_lo, c_hi = float(dist.sf(lo)), float(dist.sf(hi))
        mass = abs(c_hi - c_lo)
        if not (np.isfinite(mass) and mass > 0.0):
            raise ValueError(
                f"Prior of '{name}': the distribution has no probability mass in "
                + f"[{lo:g}, {hi:g}]"
            )
        logpdf = Prior._fast_logpdf(dist, lo, hi)
        return ("dist", lo, hi, -np.log(mass), logpdf, c_lo, c_hi, upper_tail, dist)

    @staticmethod
    def _fast_logpdf(dist, lo: float, hi: float) -> Callable[[np.ndarray], np.ndarray]:
        """``logpdf`` of ``dist`` without ``scipy``'s per-call argument checks, or the
        public ``dist.logpdf`` if the private API fails or differs (checked at points
        inside ``[lo, hi]``)"""
        try:
            gen = dist.dist
            shapes, loc, scale = gen._parse_args(*dist.args, **dist.kwds)
            raw, log_scale = gen._logpdf, float(np.log(scale))

            def logpdf(x: np.ndarray) -> np.ndarray:
                return raw((x - loc) / scale, *shapes) - log_scale

            xs = lo + (hi - lo) * np.linspace(0.05, 0.95, 7)
            with np.errstate(all="ignore"):
                ok = np.allclose(logpdf(xs), dist.logpdf(xs), rtol=1e-10, atol=1e-12)
            if ok:
                return logpdf
        # Any failure of the private API (e.g.: changed in a scipy version) falls
        # back to the public method
        except Exception:  # noqa: BLE001, S110
            pass
        return dist.logpdf

    def logpdf(self, x: np.ndarray) -> np.ndarray:
        """Log-density of the priors at the points ``x`` (shape ``(N, n_dim)``),
        ``-inf`` outside of the limits"""
        logp = np.zeros(len(x))
        for i, (kind, lo, hi, log_norm, logpdf, *_) in enumerate(self._terms):
            inside = (x[:, i] >= lo) & (x[:, i] <= hi)
            if kind == "uniform":
                logp += np.where(inside, log_norm, -np.inf)
            else:
                lp = np.full(len(x), -np.inf)
                if inside.any():
                    with np.errstate(divide="ignore"):
                        lp[inside] = logpdf(x[inside, i]) + log_norm
                logp += lp
        return logp

    def rvs(self, size: int, rng: np.random.RandomState) -> np.ndarray:
        """Random draws of the priors, shape ``(size, n_dim)``"""
        samples = []
        for kind, lo, hi, _, _, c_lo, c_hi, upper_tail, dist in self._terms:
            # One array of uniform values per parameter
            u = rng.uniform(0.0, 1.0, size)
            if kind == "uniform":
                samples.append(u * (hi - lo) + lo)
                continue
            c = c_lo + u * (c_hi - c_lo)
            x = dist.isf(c) if upper_tail else dist.ppf(c)
            samples.append(np.clip(x, lo, hi))
        return np.transpose(samples)

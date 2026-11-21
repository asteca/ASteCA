"""Data used to generate synthetic clusters and to evaluate their likelihood.

:py:class:`SynthEngine` holds everything a :py:class:`asteca.Synthetic` object needs
to generate a synthetic cluster, and is the single entry point used to generate them
(one or many, :py:meth:`SynthEngine.synth_one` and :py:meth:`SynthEngine.synth_batch`)
and to evaluate their likelihood (:py:meth:`SynthEngine.log_lkl_batch`, used by
:py:class:`asteca.Sampler`), with the compiled functions of :py:mod:`compiled` if
numba is installed, else with :py:func:`generate.synth_one`. Its parameters are
validated by :py:meth:`SynthEngine.validate_params`. :py:class:`ObservedData` holds the
observed cluster stored by :py:meth:`asteca.Synthetic.calibrate`. Both are created
once and not modified afterwards.
"""

from dataclasses import dataclass, field
from functools import cached_property
from typing import TYPE_CHECKING

import numpy as np

from .. import HAS_NUMBA
from . import compiled
from . import generate as gen

if TYPE_CHECKING:
    from .likelihood import Tremmel

# Order of the full parameter vector of a model (see SynthEngine.validate_params())
FULL_PARAMS = ("met", "loga", "alpha", "beta", "Av", "DR", "Rv", "dm")
_FULL_PARAMS_SET = frozenset(FULL_PARAMS)
# Order in which the parameters are checked (it defines which error is raised when
# several parameters are invalid): met, loga, alpha, beta, Av, DR, dm, Rv
_CHECK_ORDER = (0, 1, 2, 3, 4, 5, 7, 6)
# Types accepted as single values by the fast path of SynthEngine.validate_params()
_SCALARS = (int, float, np.integer, np.floating)


@dataclass(frozen=True, eq=False)
class SynthEngine:
    """Data that defines the synthetic clusters of a :py:class:`asteca.Synthetic`
    object, and the methods that validate their parameters, generate them and
    evaluate their likelihood. It is fixed once created:
    :py:meth:`asteca.Synthetic.calibrate` replaces it with a copy with the random
    values, maximum magnitude, error distribution and number of stars of the
    observed cluster (an object not calibrated uses a copy sized for the requested
    number of stars, with no magnitude cut and no errors). The
    differential reddening model is not included, as it can be changed at any time: it
    is passed to each method, which selects the corresponding random values in
    ``rand_floats``.

    :ivar met_age_dict: Metallicity and age grids of the isochrones
    :vartype met_age_dict: dict
    :ivar met_tol: Tolerance for metallicity values matching
    :vartype met_tol: float
    :ivar age_tol: Tolerance for log(age) values matching
    :vartype age_tol: float
    :ivar def_params: Default values of the fundamental parameters
    :vartype def_params: dict
    :ivar m_ini_idx: Index of the initial mass row
    :vartype m_ini_idx: int
    :ivar theor_tracks: Theoretical isochrones, with the binary systems added
    :vartype theor_tracks: np.ndarray
    :ivar ext_law: Extinction law
    :vartype ext_law: str
    :ivar ext_coefs: Extinction coefficients, see
        :py:func:`precompute.ccmo_filter_coefs` (no rows for ``GAIADR3``)
    :vartype ext_coefs: np.ndarray
    :ivar rand_floats: Random values drawn in advance, see
        :py:func:`precompute.draw_rand_floats` (empty for the engine of an object not yet
        sized for a number of stars)
    :vartype rand_floats: dict
    :ivar imf_knots: Knots ``(x, y)`` (CDF values, masses) of the inverse CDF of the
        IMF, see :py:func:`precompute.imf_cdf_tables`
    :vartype imf_knots: tuple
    :ivar imf_cdf: IMF CDF at the initial masses of each grid node of the isochrones
    :vartype imf_cdf: np.ndarray
    :ivar imf_lut: Bucket lookup tables of ``imf_cdf``
    :vartype imf_lut: np.ndarray
    :ivar max_mag: Maximum magnitude of the synthetic clusters
    :vartype max_mag: float
    :ivar err_dist: Error distributions added to the synthetic photometry, one row
        per magnitude and color (no rows: no errors), see
        :py:func:`precompute.error_distribution`
    :vartype err_dist: np.ndarray
    :ivar N_synth: Number of stars of the synthetic clusters
    :vartype N_synth: int
    """

    met_age_dict: dict
    met_tol: float
    age_tol: float
    def_params: dict
    m_ini_idx: int
    theor_tracks: np.ndarray
    ext_law: str
    ext_coefs: np.ndarray
    rand_floats: dict
    imf_knots: tuple
    imf_cdf: np.ndarray
    imf_lut: np.ndarray
    max_mag: float = np.inf
    err_dist: np.ndarray = field(default_factory=lambda: np.empty((0, 0)))
    N_synth: int = 0

    def validate_params(self, params: dict | None) -> np.ndarray:
        """Validated values of all the parameters of one or more models, in the
        order of ``FULL_PARAMS``. The parameters not in ``params`` take their default
        values. Used by both the single-model and the batched paths, so that they
        accept and reject the same values.

        :param params: Dictionary with a single value, or an array of values (one per
            model), for each parameter. Single values are shared by all the models.
            If ``None``, the default values are used
        :type params: dict | None

        :raises ValueError: If a parameter is not recognized, if the arrays have
            different lengths or are not 1D, if (z, a) are outside of the isochrones'
            ranges, if any of ``alpha, beta, Av, DR, dm`` is negative, if
            ``Rv <= 0``, or if any value is ``nan``

        :return: Array of shape ``(N_models, 8)``
        :rtype: np.ndarray
        """
        params = {} if params is None else params
        unknown = params.keys() - _FULL_PARAMS_SET
        if unknown:
            raise ValueError(f"Parameter(s) {sorted(unknown)} not recognized")

        values = [params[k] if k in params else self.def_params[k] for k in FULL_PARAMS]
        # Fast path for a single model, checked with the same comparisons in Python
        # (much faster than with arrays this small)
        lo, hi = self._bounds
        if all(isinstance(v, _SCALARS) for v in values):
            vals = [float(v) for v in values]
            if all(lo[i] <= v <= hi[i] for i, v in enumerate(vals)):
                return np.array([vals])

        P = _broadcast_params(values)
        # 'nan' values are rejected, as they fail both comparisons
        bad = ~((P >= lo) & (P <= hi))
        if bad.any():
            self._raise_invalid(P, bad)
        return P

    @cached_property
    def _bounds(self) -> tuple[list[float], list[float]]:
        """Valid ranges of the parameters. ``met, loga`` with a single value in the
        grid must be equal to it, within their tolerance"""
        lo = [0.0] * len(FULL_PARAMS)
        hi = [float("inf")] * len(FULL_PARAMS)
        # Rv must be > 0
        lo[FULL_PARAMS.index("Rv")] = float(np.nextafter(0.0, 1.0))
        for i, (k, tol) in enumerate((("met", self.met_tol), ("loga", self.age_tol))):
            grid = self.met_age_dict[k]
            if len(grid) > 1:
                lo[i], hi[i] = float(grid[0]), float(grid[-1])
            else:
                lo[i], hi[i] = float(grid[0]) - tol, float(grid[0]) + tol
        return lo, hi

    @cached_property
    def sampling_ranges(self) -> dict[str, tuple[float, float]]:
        """Valid range ``(low, high)`` of each parameter, used to sample parameter
        values. ``met, loga`` with a single value in the grid have ``low = high``
        (equal to that value)"""
        lo, hi = self._bounds
        lims = {k: (lo[i], hi[i]) for i, k in enumerate(FULL_PARAMS)}
        for k in ("met", "loga"):
            grid = self.met_age_dict[k]
            if len(grid) == 1:
                lims[k] = (float(grid[0]), float(grid[0]))
        return lims

    def _raise_invalid(self, P: np.ndarray, bad: np.ndarray) -> None:
        """Raise the error for the first invalid parameter

        :param P: Array of parameters
        :type P: np.ndarray
        :param bad: Mask of the invalid values in ``P``
        :type bad: np.ndarray

        :raises ValueError: Always
        """
        for i in _CHECK_ORDER:
            if not bad[:, i].any():
                continue
            k, v = FULL_PARAMS[i], P[bad[:, i], i][0]
            if k in ("met", "loga"):
                grid = self.met_age_dict[k]
                if len(grid) == 1:
                    raise ValueError(f"Parameter {k}={v:.5f} != fixed value {grid[0]}")
                raise ValueError(
                    f"Parameter {k}={v:.5f} outside isochrone"
                    + f" range [{grid[0]}, {grid[-1]}]"
                )
            if k == "Rv":
                raise ValueError(f"Parameter Rv={v:.5f} must be > 0")
            raise ValueError(f"Parameter {k}={v:.5f} must be >= 0")

    def synth_one(
        self,
        p: np.ndarray,
        dr_model: str,
        show_errors: bool = False,
        isoch_flag: bool = False,
        isoch_dr_frac: float = 0.0,
        mass_info: dict | None = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        """One synthetic cluster (and its isochrone), identical to
        :py:func:`generate.synth_one` (which is used if the compiled functions
        are not available). The synthetic cluster is generated by the compiled
        :py:func:`compiled.synth_one_nb`, the isochrone and ``mass_info`` are
        obtained from the data it returns.

        :param p: Values of the fundamental parameters: a row of the array
            returned by :py:meth:`SynthEngine.validate_params
            <engine.SynthEngine.validate_params>`
        :type p: np.ndarray
        :param dr_model: Differential reddening model
        :type dr_model: str
        :param show_errors: Print the reason when an empty cluster is returned
        :type show_errors: bool
        :param isoch_flag: Also return the isochrone
        :type isoch_flag: bool
        :param isoch_dr_frac: Fraction of the differential reddening applied to the
            returned isochrone
        :type isoch_dr_frac: float
        :param mass_info: Filled with the data required to estimate the total mass,
            see :py:func:`generate.synth_one`
        :type mass_info: dict | None

        :return: Isochrone (empty if ``isoch_flag=False``) and synthetic cluster (an
            empty 1D array if the cluster is empty)
        :rtype: tuple[np.ndarray, np.ndarray]
        """
        if not HAS_NUMBA:
            return gen.synth_one(
                self,
                p,
                dr_model,
                show_errors=show_errors,
                isoch_flag=isoch_flag,
                isoch_dr_frac=isoch_dr_frac,
                mass_info=mass_info,
            )

        out = compiled.synth_one_nb(p, *compiled.synth_args(self, dr_model))
        isoch = gen.isochrone_at_dr_frac(
            self, p, dr_model, isoch_flag, isoch_dr_frac, out.isoch_moved, out.isoch_av
        )
        if out.status != 0:
            if show_errors:
                print(gen.empty_message(out.status, p))
            return isoch, np.array([])

        if mass_info is not None:
            gen.fill_mass_info(mass_info, self, p, out)
        return isoch, out.synth

    def synth_batch(self, P: np.ndarray, dr_model: str) -> list[np.ndarray]:
        """Synthetic clusters of a batch of models, identical to those generated
        one at a time by :py:func:`generate.synth_one`

        :param P: Array of parameters, see :py:meth:`SynthEngine.validate_params
            <engine.SynthEngine.validate_params>`
        :type P: np.ndarray
        :param dr_model: Differential reddening model
        :type dr_model: str

        :return: One synthetic cluster per model (an empty 1D array if the cluster
            is empty)
        :rtype: list[np.ndarray]
        """
        if not HAS_NUMBA:
            return [gen.synth_one(self, p, dr_model)[1] for p in P]
        out = compiled.synth_batch_nb(P, *compiled.synth_args(self, dr_model))
        # An empty synthetic cluster is returned as an empty 1D array, as with
        # synth_one()
        return [a if a.shape[1] else np.array([]) for a in out]

    def log_lkl_batch(
        self,
        P: np.ndarray,
        dr_model: str,
        lkl: "Tremmel",
        temperature: float,
        parallel: bool = False,
        n_threads: int | None = None,
    ) -> np.ndarray:
        """Log-likelihoods of a batch of models, identical to
        ``lkl.log_lkl(generate.synth_one(self, p, dr_model)[1], temperature)`` for
        each row ``p`` of ``P``. With the compiled functions, each synthetic cluster
        is generated and evaluated in the same call (the synthetic clusters never go
        back to Python, which is where the gain comes from), optionally distributing
        the models among parallel threads.

        :param P: Array of parameters, see :py:meth:`SynthEngine.validate_params
            <engine.SynthEngine.validate_params>`
        :type P: np.ndarray
        :param dr_model: Differential reddening model
        :type dr_model: str
        :param lkl: Likelihood of the observed cluster
        :type lkl: Tremmel
        :param temperature: Likelihood temperature, assumed positive (validated by
            :py:class:`asteca.Sampler`)
        :type temperature: float
        :param parallel: Evaluate the models in parallel threads (only used if the
            compiled functions are available)
        :type parallel: bool
        :param n_threads: Number of threads used if ``parallel=True``, assumed valid
            (validated by :py:class:`asteca.Sampler`). If ``None``, numba's current
            setting is used
        :type n_threads: int | None

        :return: Array with the log-likelihood of each model, ``0`` is the best value
        :rtype: np.ndarray
        """
        if not HAS_NUMBA:
            return np.array(
                [
                    lkl.log_lkl(gen.synth_one(self, p, dr_model)[1], temperature)
                    for p in P
                ]
            )

        args = (
            P,
            *compiled.synth_args(self, dr_model),
            lkl.kernel_args,
            float(temperature),
        )
        if not parallel:
            return compiled.log_lkl_batch_nb(*args)
        if n_threads is None:
            return compiled.log_lkl_batch_parallel_nb(*args)

        import numba

        # Set the number of threads for this call only
        prev_threads = numba.get_num_threads()
        numba.set_num_threads(int(n_threads))
        try:
            return compiled.log_lkl_batch_parallel_nb(*args)
        finally:
            numba.set_num_threads(prev_threads)

    @cached_property
    def kernel_data(self) -> dict:
        """Arrays used by the compiled batch functions, obtained the first time they
        are needed (see :py:func:`compiled.kernel_data`). Requires numba"""
        return compiled.kernel_data(self)


@dataclass(frozen=True, eq=False)
class ObservedData:
    """Data of the observed cluster stored by :py:meth:`asteca.Synthetic.calibrate`.
    The data that defines its synthetic clusters (maximum magnitude, error
    distribution, number of stars) is stored in the :py:class:`SynthEngine`

    :ivar mag: Magnitudes
    :vartype mag: np.ndarray
    :ivar e_mag: Magnitude uncertainties
    :vartype e_mag: np.ndarray
    :ivar colors: ``[color, color2]``, with ``color2=None`` if not defined
    :vartype colors: list
    :ivar e_colors: ``[e_color, e_color2]``, with ``e_color2=None`` if not defined
    :vartype e_colors: list
    :ivar star_id: Identifiers of the stars (can be ``None``)
    :vartype star_id: np.ndarray | None
    :ivar ra: Right ascensions (can be ``None``)
    :vartype ra: np.ndarray | None
    :ivar dec: Declinations (can be ``None``)
    :vartype dec: np.ndarray | None
    :ivar lkl: Likelihood of synthetic clusters given the observed one
    :vartype lkl: Tremmel
    """

    mag: np.ndarray
    e_mag: np.ndarray
    colors: list
    e_colors: list
    star_id: np.ndarray | None
    ra: np.ndarray | None
    dec: np.ndarray | None
    lkl: "Tremmel"


def _broadcast_params(values: list) -> np.ndarray:
    """Array of parameters from a mix of single values and arrays

    :param values: Value(s) of each parameter, in the order of ``FULL_PARAMS``
    :type values: list

    :raises ValueError: If an array is not 1D, or if the arrays have different
        lengths

    :return: Array of shape ``(N_models, 8)``
    :rtype: np.ndarray
    """
    arrays = [np.asarray(v, dtype=float) for v in values]
    for k, v in zip(FULL_PARAMS, arrays):
        if v.ndim > 1:
            raise ValueError(f"Parameter {k} must be a single value or a 1D array")
    lengths = {v.size for v in arrays if v.ndim == 1 and v.size != 1}
    if len(lengths) > 1:
        raise ValueError("All the parameter arrays must have the same length")
    N = lengths.pop() if lengths else 1
    # Filled by rows (contiguous) and transposed once, faster than by columns
    P = np.empty((len(FULL_PARAMS), N))
    for i, v in enumerate(arrays):
        P[i] = v
    return P.T.copy()

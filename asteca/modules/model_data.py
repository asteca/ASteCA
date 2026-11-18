"""Data used to generate synthetic clusters and to evaluate their likelihood.

:py:class:`ModelEngine` holds everything a :py:class:`asteca.Synthetic` object needs
to generate a synthetic cluster, and is the single entry point used to generate them
(one or many, :py:meth:`ModelEngine.get_model` and :py:meth:`ModelEngine.generate`)
and to evaluate their likelihood (:py:meth:`ModelEngine.log_lkl`, used by
:py:class:`asteca.Sampler`), with the compiled functions of :py:mod:`synth_batch` if
numba is installed, else with :py:func:`synth_model.get_model`. Its parameters are
validated by :py:meth:`ModelEngine.params_array`. :py:class:`ObservedData` holds what
:py:meth:`asteca.Synthetic.calibrate` obtains from the observed cluster. Both are
created once and not modified afterwards.
"""

from dataclasses import dataclass
from functools import cached_property
from typing import TYPE_CHECKING

import numpy as np

from . import HAS_NUMBA
from . import synth_batch as sb
from . import synth_model as sm

if TYPE_CHECKING:
    from .likelihood_priv import Tremmel

# Order of the full parameter vector of a model (see ModelEngine.params_array())
FULL_PARAMS = ("met", "loga", "alpha", "beta", "Av", "DR", "Rv", "dm")
_FULL_PARAMS_SET = frozenset(FULL_PARAMS)
# Order in which the parameters are checked (it defines which error is raised when
# several parameters are invalid): met, loga, alpha, beta, Av, DR, dm, Rv
_CHECK_ORDER = (0, 1, 2, 3, 4, 5, 7, 6)
# Types accepted as single values by the fast path of ModelEngine.params_array()
_SCALARS = (int, float, np.integer, np.floating)


@dataclass(frozen=True, eq=False)
class ModelEngine:
    """Data that defines the synthetic clusters of a :py:class:`asteca.Synthetic`
    object (fixed once the object is created), and the methods that validate their
    parameters, generate them and evaluate their likelihood. The differential
    reddening model is not included, as it can be changed at any time: it is passed to
    each method, which selects the corresponding random values in ``rand_floats``.

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
    :ivar ext_coefs: Extinction coefficients (empty for ``GAIADR3``)
    :vartype ext_coefs: list
    :ivar rand_floats: Random values drawn in advance
    :vartype rand_floats: dict
    :ivar st_dist_mass: Sampled IMF masses for each (z, a) box
    :vartype st_dist_mass: list
    :ivar mass_idx: Search indexes of the sampled IMF masses into the mass grid of
        the isochrones, see :py:func:`synth_cluster_priv.mass_idx_lookup`
    :vartype mass_idx: dict
    """

    met_age_dict: dict
    met_tol: float
    age_tol: float
    def_params: dict
    m_ini_idx: int
    theor_tracks: np.ndarray
    ext_law: str
    ext_coefs: list
    rand_floats: dict
    st_dist_mass: list
    mass_idx: dict

    def params_array(self, params: dict | None) -> np.ndarray:
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
            ok = all(lo[i] <= v <= hi[i] for i, v in enumerate(vals)) and all(
                abs(vals[i] - g) <= tol for i, g, tol in self._fixed
            )
            if ok:
                return np.array([vals])

        P = _broadcast_params(values)
        # 'nan' values are rejected, as they fail both comparisons
        bad = ~((P >= lo) & (P <= hi))
        for i, g, tol in self._fixed:
            bad[:, i] = ~(np.abs(P[:, i] - g) <= tol)
        if bad.any():
            self._raise_invalid(P, bad)
        return P

    @cached_property
    def _bounds(self) -> tuple[list[float], list[float]]:
        """Valid ranges of the parameters (``met, loga`` with a single value in the
        grid are checked separately, see :py:attr:`_fixed`)"""
        lo = [0.0] * len(FULL_PARAMS)
        hi = [float("inf")] * len(FULL_PARAMS)
        # Rv must be > 0
        lo[FULL_PARAMS.index("Rv")] = float(np.nextafter(0.0, 1.0))
        for i, k in enumerate(("met", "loga")):
            grid = self.met_age_dict[k]
            if len(grid) > 1:
                lo[i], hi[i] = float(grid[0]), float(grid[-1])
            else:
                lo[i], hi[i] = -float("inf"), float("inf")
        return lo, hi

    @cached_property
    def _fixed(self) -> tuple[tuple[int, float, float], ...]:
        """Index, grid value and tolerance of ``met, loga`` if their grid has a
        single value: the parameter must then be equal to it, within the tolerance"""
        return tuple(
            (i, float(self.met_age_dict[k][0]), tol)
            for i, (k, tol) in enumerate(
                (("met", self.met_tol), ("loga", self.age_tol))
            )
            if len(self.met_age_dict[k]) == 1
        )

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

    @property
    def compiled(self) -> bool:
        """``True`` if the compiled functions of :py:mod:`synth_batch` are used (numba
        installed), else the Python implementation :py:func:`synth_model.get_model`
        is used. The results are identical"""
        return HAS_NUMBA

    def check_threads(self, parallel: bool, n_threads: int | None) -> None:
        """Validate ``n_threads``, only if the threads are used (compiled functions
        and ``parallel=True``)

        :param parallel: Evaluate the log-likelihoods in parallel threads
        :type parallel: bool
        :param n_threads: Number of threads, or ``None`` for numba's current setting
        :type n_threads: int | None

        :raises ValueError: If ``n_threads`` is not an integer between ``1`` and
            ``NUMBA_NUM_THREADS``
        """
        if not (self.compiled and parallel and n_threads is not None):
            return
        import numba

        if (
            isinstance(n_threads, bool)
            or not isinstance(n_threads, (int, np.integer))
            or n_threads < 1
        ):
            raise ValueError(
                f"'n_threads' must be a positive integer, got {n_threads!r}"
            )
        # NUMBA_NUM_THREADS is set dynamically in numba.config (unknown to type
        # checkers)
        max_threads = numba.config.NUMBA_NUM_THREADS  # pyright: ignore[reportPrivateImportUsage, reportAttributeAccessIssue]
        if n_threads > max_threads:
            raise ValueError(
                f"'n_threads={n_threads}' can not be larger than NUMBA_NUM_THREADS="
                + f"{max_threads}"
            )

    def get_model(
        self,
        p: np.ndarray,
        dr_model: str,
        calib: tuple[float, list[np.ndarray], int],
        show_errors: bool = False,
        isoch_flag: bool = False,
        isoch_dr_frac: float = 0.0,
        mass_info: dict | None = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        """One synthetic cluster (and its isochrone), identical to
        :py:func:`synth_model.get_model` (which is used if the compiled functions
        are not available). The synthetic cluster is generated by the compiled
        :py:func:`synth_batch._synth_one`, the isochrone and ``mass_info`` are
        obtained from the data it returns.

        :param p: Values of the fundamental parameters: a row of the array
            returned by :py:meth:`ModelEngine.params_array
            <model_data.ModelEngine.params_array>`
        :type p: np.ndarray
        :param dr_model: Differential reddening model
        :type dr_model: str
        :param calib: Maximum magnitude, error distribution and number of stars of
            the synthetic cluster, as given by ``Synthetic._calib_values()``
        :type calib: tuple[float, list[np.ndarray], int]
        :param show_errors: Print the reason when an empty cluster is returned
        :type show_errors: bool
        :param isoch_flag: Also return the isochrone
        :type isoch_flag: bool
        :param isoch_dr_frac: Fraction of the differential reddening applied to the
            returned isochrone
        :type isoch_dr_frac: float
        :param mass_info: Filled with the data required to estimate the total mass,
            see :py:func:`synth_model.get_model`
        :type mass_info: dict | None

        :return: Isochrone (empty if ``isoch_flag=False``) and synthetic cluster (an
            empty 1D array if the cluster is empty)
        :rtype: tuple[np.ndarray, np.ndarray]
        """
        if not self.compiled:
            return sm.get_model(
                self,
                p,
                dr_model,
                calib,
                show_errors=show_errors,
                isoch_flag=isoch_flag,
                isoch_dr_frac=isoch_dr_frac,
                mass_info=mass_info,
            )

        out = sb._synth_one(p, *sb.batch_args(self, dr_model, calib))
        isoch = sm.output_isochrone(
            self, p, dr_model, isoch_flag, isoch_dr_frac, out.isoch_moved, out.isoch_av
        )
        if out.status != 0:
            if show_errors:
                print(sm.empty_message(out.status, p))
            return isoch, np.array([])

        if mass_info is not None:
            sm.fill_mass_info(
                mass_info,
                self,
                p,
                calib[2],
                out.isoch_moved,
                out.idx_cut,
                out.N_pool,
                out.keep,
            )
        return isoch, out.synth

    def generate(
        self,
        P: np.ndarray,
        dr_model: str,
        calib: tuple[float, list[np.ndarray], int],
    ) -> list[np.ndarray]:
        """Synthetic clusters of a batch of models, identical to those generated
        one at a time by :py:func:`synth_model.get_model`

        :param P: Array of parameters, see :py:meth:`ModelEngine.params_array
            <model_data.ModelEngine.params_array>`
        :type P: np.ndarray
        :param dr_model: Differential reddening model
        :type dr_model: str
        :param calib: Maximum magnitude, error distribution and number of stars of
            the synthetic clusters, as given by ``Synthetic._calib_values()``
        :type calib: tuple[float, list[np.ndarray], int]

        :return: One synthetic cluster per model (an empty 1D array if the cluster
            is empty)
        :rtype: list[np.ndarray]
        """
        if not self.compiled:
            return [sm.get_model(self, p, dr_model, calib)[1] for p in P]
        out = sb._generate_batch(P, *sb.batch_args(self, dr_model, calib))
        # An empty synthetic cluster is returned as an empty 1D array, as with
        # get_model()
        return [a if a.shape[1] else np.array([]) for a in out]

    def log_lkl(
        self,
        P: np.ndarray,
        dr_model: str,
        obs: "ObservedData",
        temperature: float,
        parallel: bool = False,
        n_threads: int | None = None,
    ) -> np.ndarray:
        """Log-likelihoods of a batch of models, identical to
        ``obs.lkl.log_lkl(synth_model.get_model(self, p, dr_model, obs.calib_values())[1],
        temperature)`` for each row ``p`` of ``P``. With the compiled functions,
        each synthetic cluster is generated and evaluated in the same call (the
        synthetic clusters never go back to Python, which is where the gain comes
        from), optionally distributing the models among parallel threads.

        :param P: Array of parameters, see :py:meth:`ModelEngine.params_array
            <model_data.ModelEngine.params_array>`
        :type P: np.ndarray
        :param dr_model: Differential reddening model
        :type dr_model: str
        :param obs: Data of the observed cluster
        :type obs: ObservedData
        :param temperature: Likelihood temperature, assumed positive (validated by
            :py:class:`asteca.Sampler`)
        :type temperature: float
        :param parallel: Evaluate the models in parallel threads (only used if the
            compiled functions are available)
        :type parallel: bool
        :param n_threads: Number of threads used if ``parallel=True``, assumed valid
            (see :py:meth:`check_threads`). If ``None``, numba's current setting is
            used
        :type n_threads: int | None

        :return: Array with the log-likelihood of each model, ``0`` is the best value
        :rtype: np.ndarray
        """
        lk = obs.lkl
        calib = obs.calib_values()
        if not self.compiled:
            return np.array(
                [
                    lk.log_lkl(sm.get_model(self, p, dr_model, calib)[1], temperature)
                    for p in P
                ]
            )

        args: tuple = (
            P,
            *sb.batch_args(self, dr_model, calib),
            lk.m_ini_idx,
            lk.hist_args,
            lk.nz_pos,
            lk.row_off,
            lk.table_flat,
            lk.m_max,
            lk.min_lkl,
            lk.max_lkl,
            lk.lkl_range,
            float(lk.obs_mag_median),
            float(lk.obs_col_median),
            float(lk.cmd_widths[0]),
            float(lk.cmd_widths[1]),
            float(temperature),
        )
        if not parallel:
            return sb._lnl_batch_serial(*args)
        if n_threads is None:
            return sb._lnl_batch_parallel(*args)

        import numba

        # Set the number of threads for this call only
        prev_threads = numba.get_num_threads()
        numba.set_num_threads(int(n_threads))
        try:
            return sb._lnl_batch_parallel(*args)
        finally:
            numba.set_num_threads(prev_threads)

    @cached_property
    def batch_data(self) -> dict:
        """Arrays used by the compiled batch functions, obtained the first time they
        are needed (see :py:func:`synth_batch.batch_data`). Requires numba"""
        return sb.batch_data(self)


@dataclass(frozen=True, eq=False)
class ObservedData:
    """Data of the observed cluster stored by :py:meth:`asteca.Synthetic.calibrate`

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
    :ivar max_mag: Maximum magnitude of the synthetic clusters
    :vartype max_mag: float
    :ivar N_stars: Number of stars of the synthetic clusters
    :vartype N_stars: int
    :ivar err_dist: Error distribution of the synthetic clusters
    :vartype err_dist: list[np.ndarray]
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
    max_mag: float
    N_stars: int
    err_dist: list[np.ndarray]
    lkl: "Tremmel"

    def calib_values(self) -> tuple[float, list[np.ndarray], int]:
        """Maximum magnitude, error distribution and number of stars of the
        synthetic clusters

        :return: Maximum magnitude, error distribution and number of stars
        :rtype: tuple[float, list[np.ndarray], int]
        """
        return self.max_mag, self.err_dist, self.N_stars


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

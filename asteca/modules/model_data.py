"""Data used to generate synthetic clusters and to evaluate their likelihood.

:py:class:`ModelContext` holds everything a :py:class:`asteca.Synthetic` object needs
to generate a synthetic cluster, and :py:class:`ObservedData` holds what
:py:meth:`asteca.Synthetic.calibrate` obtains from the observed cluster. Both are
created once and not modified afterwards; the single-model
(:py:func:`synth_model.get_model`) and the batched (:py:mod:`synth_batch`) paths
read their data from them.
"""

from dataclasses import dataclass
from functools import cached_property
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from .likelihood_priv import Tremmel


@dataclass(frozen=True, eq=False)
class ModelContext:
    """Data that defines the synthetic clusters of a :py:class:`asteca.Synthetic`
    object, fixed once the object is created. The differential reddening model is
    not included, as it can be changed at any time: it is passed to
    :py:func:`synth_model.get_model` and :py:func:`synth_batch.batch_args`, which
    select the corresponding random values in ``rand_floats``.

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

    @cached_property
    def batch_data(self) -> dict:
        """Arrays used by the compiled batch functions, obtained the first time they
        are needed (see :py:func:`synth_batch.batch_data`). Requires numba"""
        from .synth_batch import batch_data

        return batch_data(self)


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

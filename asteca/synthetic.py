import warnings
from collections.abc import Sequence

import numpy as np

from .cluster import Cluster
from .isochrones import Isochrones
from .modules.synthetic.engine import ModelEngine, ObservedData

# Default values of the fundamental parameters, used for those not in 'def_params'
DEFAULT_PARAMS = {
    "met": 0.0152,
    "loga": 8.0,
    "alpha": 0.09,
    "beta": 0.94,
    "Rv": 3.1,
    "DR": 0.0,
    "Av": 0.2,
    "dm": 9.0,
}


class Synthetic:
    """Define a :py:class:`Synthetic` object.

    Use the isochrones loaded in the
    :py:class:`Isochrones <asteca.isochrones.Isochrones>` object
    to generate a :py:class:`Synthetic` object. This object is used
    to generate synthetic clusters given a :py:class:`Cluster <asteca.cluster.Cluster>`
    object and a set of input fundamental parameters (metallicity, age, distance,
    extinction, etc.).

    See the :ref:`synthetic_module` section for more details.

    :param isochs: :py:class:`Isochrones <asteca.isochrones.Isochrones>` object with
        the loaded files for the theoretical isochrones
    :type isochs: Isochrones
    :param def_params: Default values of the fundamental parameters, used for the
        parameters not given when generating a synthetic cluster. The parameters not
        included here (all of them if ``None`` is given) take the following values::

            {
                "met": 0.0152,
                "loga": 8.0,
                "alpha": 0.09,
                "beta": 0.94,
                "Rv": 3.1,
                "DR": 0.0,
                "Av": 0.2,
                "dm": 9.0,
            }

        A warning is shown if a default value is not valid (e.g.: a metallicity or age
        outside the range of the isochrones); the parameter must then be given each
        time a synthetic cluster is generated
    :type def_params: dict[str, float] | None
    :param IMF_name: Name of the initial mass function used to populate the isochrones,
        one of ``salpeter_1955, kroupa_2001, chabrier_2014``
    :type IMF_name: str
    :param max_mass: Maximum total initial mass. Should be large enough to allow
        generating as many synthetic stars as observed stars
    :type max_mass: int
    :param gamma: Distribution function for the mass ratio of the binary systems,
        float or one of ``D&K, fisher_stepped, fisher_peaked, raghavan``
    :type gamma: float | str
    :param ext_law: Extinction law, one of ``CCMO, GAIADR3``.
        If ``GAIADR3`` is selected, the magnitude and first
        color defined in the :py:class:`Isochrones <asteca.isochrones.Isochrones>` and
        :py:class:`Cluster <asteca.cluster.Cluster>` objects are assumed to be Gaia's
        (E)DR3 **G** and **(BP-RP)** respectively. A second color can not be used with
        this law (a ``ValueError`` is raised if one is defined); use ``CCMO`` instead.
        The ``Rv`` parameter has no effect with this law (its coefficients were derived
        for a single extinction curve, ``Rv~3.1``); a warning is shown if the default
        ``Rv`` value in ``def_params`` is not ``3.1``
    :type ext_law: str
    :param DR_distribution: Model for the differential reddening, one of
        ``foreground, internal``. ``foreground``: dust-free cluster behind a patchy
        foreground screen. ``internal``: dust mixed with the cluster stars, behind a
        uniform foreground screen. In both models the absorption of each star is
        ``Av + f * DR`` with ``f`` in ``[0, 1]``, i.e.: ``Av`` is the minimum and
        ``Av + DR`` the maximum absorption. The model can be changed at any time
        by setting the ``DR_distribution`` attribute of the object (e.g.:
        ``synthcl.DR_distribution = "internal"``), without re-creating or
        re-calibrating it. It is used by all the methods that generate synthetic
        clusters and by :py:class:`asteca.Sampler`. All the models use the same
        random values, so they can be compared without random noise
    :type DR_distribution: str
    :param seed: Random seed. If ``None`` a random integer will be generated and used
    :type seed: int | None
    :param verbose: Verbose level. A value of ``0`` hides all output
    :type verbose: int

    :raises ValueError: If any of the attributes is not recognized as a valid option
    """

    def __init__(
        self,
        isochs: Isochrones,
        def_params: dict[str, float] | None = None,
        IMF_name: str = "chabrier_2014",
        max_mass: int = 10_000,
        gamma: float | str = "D&K",
        ext_law: str = "CCMO",
        DR_distribution: str = "foreground",
        seed: int | None = None,
        verbose: int = 1,
    ) -> None:

        self.isochs = isochs
        self.def_params = {**DEFAULT_PARAMS, **(def_params or {})}
        unknown = self.def_params.keys() - DEFAULT_PARAMS.keys()
        if unknown:
            raise ValueError(
                f"Parameter(s) {sorted(unknown)} in 'def_params' not recognized"
            )
        self.IMF_name = IMF_name
        self.max_mass = max_mass
        self.gamma = gamma
        self.ext_law = ext_law
        # Validated by the property setter
        self.DR_distribution = DR_distribution
        # Set seed
        self.seed: int = (
            int(np.random.randint(0, 100_000)) if seed is None else int(seed)
        )
        self.verbose = verbose
        self.rng = np.random.default_rng(self.seed)

        # Check gamma distribution
        gammas = ("D&K", "fisher_stepped", "fisher_peaked", "raghavan")
        if isinstance(self.gamma, str) and self.gamma not in gammas:
            raise ValueError(
                f"gamma '{self.gamma}' not recognized. Should be one of {gammas}"
            )

        from .modules.synthetic import precompute as pre

        # Check IMF function
        imfs = ("salpeter_1955", "kroupa_2001", "chabrier_2014")
        if self.IMF_name not in imfs:
            raise ValueError(
                f"IMF '{self.IMF_name}' not recognized. Should be one of {imfs}"
            )

        # Check extinction law
        ext_laws = ("CCMO", "GAIADR3")
        if self.ext_law not in ext_laws:
            raise ValueError(
                f"Extinction law '{self.ext_law}' not recognized. Should be "
                + f"one of {ext_laws}"
            )

        self._vp("\nInstantiating synthetic")

        # Sample the selected IMF
        Nmets, Nages = self.isochs.theor_tracks.shape[:2]
        self._vp("(Sampling IMF)", 1)
        st_dist_mass = pre.sample_imf(
            self.rng, self.IMF_name, self.max_mass, Nmets, Nages
        )

        # Add binary systems
        theor_tracks = pre.add_binarity(
            self.rng,
            self.gamma,
            self.isochs.color,
            self.isochs.color2,
            self.isochs.theor_tracks,
            self.isochs.color_filters,
        )

        # Index of the initial mass row of the synthetic clusters
        m_ini_idx = 2  # (0->mag, 1->color, 2->mass_ini)
        if self.isochs.color2_effl is not None:
            m_ini_idx = 3  # (0->mag, 1->color, 2->color2, 3->mass_ini)
        # Set tolerances for the ranges of the isochrones in `generate()`
        age_tol = 0.01
        met_tol = 0.001  # z
        if self.isochs.z_to_FeH is not None:
            met_tol = 0.01  # FeH

        if self.ext_law in ("GAIADR3",):
            if (
                self.isochs.magnitude_effl is not None
                or self.isochs.color_effl is not None
            ):
                warnings.warn(
                    f"\nExtinction law '{self.ext_law}' does not require effective lambda"
                    + " values for the\nmagnitude and first color (assumed to be 'G' "
                    + "and 'BP-RP', respectively)."
                )

            if self.isochs.color2_effl is not None:
                raise ValueError(
                    f"\nExtinction law '{self.ext_law}' does not apply to the "
                    "second color.\nUse the 'CCMO' extinction law if you need to "
                    "process a second color"
                )

            # The Rv parameter is not used by this law
            if self.def_params.get("Rv", 3.1) != 3.1:
                warnings.warn(
                    f"\nExtinction law '{self.ext_law}' does not use the 'Rv' parameter: "
                    + f"the default value Rv={self.def_params['Rv']} will have no effect."
                )

        # Get extinction coefficients for these filters
        ext_coefs = []  # It is important to pass an empty list because the
        # function `extinction()` checks its length later on
        if self.ext_law == "CCMO":
            if self.isochs.magnitude_effl is None or self.isochs.color_effl is None:
                raise ValueError(
                    f"Extinction law '{self.ext_law}' requires effective lambda\n"
                    + "values for the magnitude and first color."
                )
            ext_coefs = pre.ccmo_ext_coeffs(
                self.isochs.magnitude_effl,
                self.isochs.color_effl,
                self.isochs.color2_effl,
            )

        # Data that defines the synthetic clusters, and the single entry point used to
        # generate them (one or many) and to evaluate their likelihood (Sampler)
        self._engine = ModelEngine(
            met_age_dict=self.isochs.met_age_dict,
            met_tol=met_tol,
            age_tol=age_tol,
            def_params=self.def_params,
            m_ini_idx=m_ini_idx,
            theor_tracks=theor_tracks,
            ext_law=self.ext_law,
            ext_coefs=ext_coefs,
            # Random floats used by `generate()`
            rand_floats=pre.randVals(self.rng, theor_tracks, st_dist_mass),
            st_dist_mass=st_dist_mass,
            # Pre-computed search indexes used by mass_interp()
            mass_idx=pre.mass_idx_lookup(theor_tracks, st_dist_mass, m_ini_idx),
        )
        # Observed cluster data, set by calibrate()
        self._obs: ObservedData | None = None

        # Turn-off points used to define the BSS region, stored by the last call to
        # cluster_properties()
        self.turn_off_points: dict | None = None

        # Validate the default values with the same checks as the parameters of the
        # synthetic clusters. An invalid value is only an error if it is used
        try:
            self._engine.params_array(None)
        except ValueError as err:
            warnings.warn(
                f"Invalid 'def_params': {err}. The parameter must be included in "
                + "'params' when generating synthetic clusters"
            )

        self._vp(f"Default params : {self.def_params}", 1)
        self._vp(f"Extinction law : {self.ext_law}", 1)
        self._vp(f"Diff reddening : {self.DR_distribution}", 1)
        self._vp(f"IMF            : {self.IMF_name}", 1)
        self._vp(f"Max init mass  : {self.max_mass}", 1)
        self._vp(f"Gamma dist     : {self.gamma}", 1)
        self._vp(f"Random seed    : {self.seed}", 1)
        self._vp("Synthetic clusters object generated")

    @property
    def DR_distribution(self) -> str:
        """Differential reddening model used to generate the synthetic clusters, one
        of ``foreground, internal`` (see :py:class:`Synthetic`). It can be changed at
        any time, e.g.: ``synthcl.DR_distribution = "internal"``. Setting a model that
        is not recognized raises a ``ValueError``
        """
        return self._DR_distribution

    @DR_distribution.setter
    def DR_distribution(self, name: str) -> None:
        from .modules.synthetic import precompute as pre

        if name not in pre.DR_MODELS:
            raise ValueError(
                f"Differential reddening model '{name}' not recognized. Should be "
                + f"one of {pre.DR_MODELS}"
            )
        self._DR_distribution = name

    def _vp(self, mssg: str, level: int = 0) -> None:
        """Verbose print method"""
        if self.verbose > level:
            print(mssg)

    def calibrate(self, cluster: Cluster):
        """Calibrate a :py:class:`Synthetic` object based on a
        :py:class:`Cluster` object .

        Use the data obtained from your observed cluster, stored in the
        :py:class:`Cluster` object, to calibrate a :py:class:`Synthetic` object.
        The observed cluster is also binned into a Hess diagram (with Knuth's rule
        for the number of bins), used by the likelihood of the
        :py:class:`asteca.Sampler` (see :ref:`sampler_module`).

        See the :ref:`synthetic_module` section for more details.

        :param cluster: :py:class:`Cluster` object with the
            processed data from your observed cluster
        :type cluster: Cluster

        :raises ValueError:
            - If required 'Cluster' attributes are missing
            - If the number of colors defined in the :py:class:`Cluster` and
              :py:class:`Synthetic` objects do not match
        """
        from .modules.synthetic import likelihood as lh
        from .modules.synthetic import precompute as pre

        required = ("mag", "e_mag", "color", "e_color")
        values = [getattr(cluster, attr) for attr in required]
        for attr, v in zip(required, values):
            if v is None:
                raise ValueError(
                    f"Attribute '{attr}' is required to calibrate 'Cluster' object."
                )
        mag, e_mag, color, e_color = values

        # Check that the number of colors match
        if self.isochs.color2_effl is not None and cluster.color2 is None:
            raise ValueError(
                "Two colors were defined in 'Synthetic' but a single color\n"
                + "was defined in 'cluster'."
            )
        if self.isochs.color2_effl is None and cluster.color2 is not None:
            raise ValueError(
                "Two colors were defined in 'cluster' but a single color\n"
                + "was defined in 'synthetic'."
            )

        # Warning with hard coded values
        f_stars_mass = 7
        if mag.size * f_stars_mass > self.max_mass:
            warnings.warn(
                (
                    f"\nThe number of stars in this `Cluster` object "
                    f"({mag.size}) is large.\nThis may lead to synthetic "
                    "clusters with insufficient stars to\n"
                    "properly reproduce the observations.\n"
                    "Consider increasing the 'max_mass' parameter in the Synthetic "
                    f"cluster \n object (current value: {self.max_mass}). As a rule "
                    f"of thumb, set 'max_mass' \nto at least {f_stars_mass} times "
                    "the number of observed stars."
                ),
                UserWarning,
            )

        # All the data is obtained before it is stored, so that the object is left
        # unchanged if any step fails
        obs = ObservedData(
            mag=mag,
            e_mag=e_mag,
            colors=[color, cluster.color2],
            e_colors=[e_color, cluster.e_color2],
            star_id=cluster.star_id,
            ra=cluster.ra,
            dec=cluster.dec,
            # Calibrate parameters using the observed cluster
            max_mag=float(np.nanmax(mag)),
            N_stars=cluster.N_stars - cluster.N_nan_phot,
            err_dist=pre.error_distribution(
                mag,
                e_mag,
                e_color,
                cluster.e_color2,
                self._engine.rand_floats["norm"][1],
            ),
            # Likelihood used by asteca.Sampler, with the defined colors
            lkl=lh.Tremmel(mag, [c for c in (color, cluster.color2) if c is not None]),
        )
        self._obs = obs

        self._vp("\nCalibrated observed cluster")
        self._vp(f"N_stars_obs    : {obs.N_stars}", 1)
        self._vp(f"Max magnitude  : {obs.max_mag:.2f}", 1)
        self._vp("Error distribution loaded", 1)

    def generate(
        self,
        params: dict | None = None,
        N_stars: int = 100,
        batch: bool = False,
        show_errors: bool = False,
    ) -> np.ndarray | list[np.ndarray]:
        """Generate a synthetic cluster, or a batch of synthetic clusters.

        The synthetic cluster is generated according to the parameters given in
        the ``params`` dictionary and the already calibrated
        :py:class:`Synthetic` object.

        The ``params`` dictionary can include values for all or some of the parameters.
        E.g. of a dictionary with some fundamental parameters:
        ``params = {met: 0.0152, loga: 8.1, Av: 0.2, DR: 0., Rv: 3.1, dm: 9.7}``
        If no ``params`` dictionary is given, an empty one is used and the default values
        defined in the :py:class:`Synthetic` object are used to generate the synthetic
        cluster.

        If the object was calibrated with :py:meth:`calibrate`, the synthetic cluster
        contains up to the same number of stars as the observed cluster (with valid
        photometry), and is cut at its maximum magnitude and assigned its
        uncertainties. Otherwise, ``N_stars`` stars are generated with no magnitude
        cut and no uncertainties.

        With ``batch=True`` many synthetic clusters are generated in a single call.
        The ``params`` dictionary then holds an array of values for each parameter,
        all with the same length (the number of models), or a single value shared by
        all the models. Parameters not included take the default values defined in the
        :py:class:`Synthetic` object. E.g., for an array ``x`` of shape
        ``(N_models, N_params)`` with the values of the parameters in ``names``:
        ``params = dict(zip(names, x.T))``. The results are identical to calling this
        method with ``batch=False`` for each set of parameters, but faster: with
        `numba <https://numba.pydata.org/>`__ installed the whole batch is generated
        in a single compiled call, which removes the Python overhead of each call (the
        gain is largest for clusters with few stars). **Without numba there is no
        speedup**: each model is generated as with ``batch=False``, which is
        equivalent to (and as fast as) a loop.

        This method is the entry point to combine the synthetic clusters with a
        likelihood or distance and a sampler of the user's choosing (see
        :ref:`own_lkl_sampler`).

        A ``ValueError`` is raised if a parameter is not recognized, if ``met`` or
        ``loga`` are outside of the ranges of the loaded isochrones, if any of
        ``alpha, beta, Av, DR, dm`` is negative, if ``Rv<=0``, or if any value is
        ``nan`` (for any of the models if ``batch=True``). The same values are
        accepted with ``batch=True`` and ``batch=False``.

        :param params: Dictionary containing the values for the fundamental parameters,
            a single value per parameter, or an array of values (one per model) if
            ``batch=True``.
        :type params: dict | None
        :param N_stars: Number of synthetic stars to generate. Only used if the
            object was not calibrated; after calling :py:meth:`calibrate` the
            synthetic clusters contain up to as many stars as the observed cluster
            (with valid photometry) and this argument is ignored.
        :type N_stars: int
        :param batch: If ``True``, generate one synthetic cluster for each set of
            parameters in ``params`` and return them in a list.
        :type batch: bool
        :param show_errors: If ``True``, possible errors will be printed to the screen.
            Useful for debugging. Only used if ``batch=False``.
        :type show_errors: bool

        :return: Returns a ``np.array`` containing a synthetic cluster
            with the data ``[mag, c1, (c2), mass, mass_b]``, where ``mag`` is
            the magnitude, ``c1`` is the color, ``c2`` is the optional second color,
            and ``mass, mass_b`` are the masses of the single and secondary components
            of the binary systems, respectively (if generated). If the system is a
            single star, then ``mass_b==np.nan``. If ``batch=True``, a list with one
            such array per model.
        :rtype: np.ndarray | list[np.ndarray]
        """
        if batch:
            return self._engine.generate(
                self._engine.params_array(params),
                self.DR_distribution,
                self._calib_values(N_stars),
            )

        _, synth_arr = self._get_model(params, N_stars, show_errors=show_errors)

        return synth_arr

    def cmd_overlap(  # noqa: DOC502, RUF100
        self,
        params: dict[str, float],
        h_mag: float = 0.3,
        h_col: float = 0.075,
    ) -> tuple[float, float]:
        """Goodness of fit of the synthetic cluster generated with the ``params``
        parameters, that can be compared across observed clusters with any number of
        stars.

        The synthetic cluster is generated internally with :py:meth:`generate`, i.e.
        ``params`` is defined in the same way (e.g.: the best fit or median
        parameters returned by :py:meth:`asteca.Sampler.summary`).

        The value is ``D = 1 - 2 * K_os / (K_oo + K_ss)``, where ``K_os`` is the mean
        value of a Gaussian kernel (with widths ``h_mag`` and ``h_col``) over all the
        pairs of observed and synthetic stars, and ``K_oo`` and ``K_ss`` are the mean
        values over the pairs of distinct observed stars and of distinct synthetic
        stars. ``D`` is the normalized squared difference between the observed and
        synthetic densities in the CMD (magnitude and first color), both smoothed
        with the kernel: ``0`` means the same distribution and ``1`` means no overlap.

        Since stars are never compared with themselves, a synthetic cluster drawn
        from the same distribution as the observed one gives ``D=0`` on average for
        any number of stars, and a given mismatch gives on average the same ``D`` for
        any number of stars. Individual values fluctuate around ``0`` and can be
        negative: ``D<0`` is a match within the noise (compare it with ``D_std``),
        not an error. Values are not clipped to ``0``, as this would bias their
        average upwards. ``D<=0.25`` is a reasonable general threshold for a good
        fit (see :ref:`cmd_overlap`). The kernel widths set the scale at which the
        CMDs are compared, so the same values must be used for all the clusters
        whose ``D`` values are compared.

        For large clusters (more than ~550 observed and synthetic stars) the sums
        over the pairs of stars are computed on a grid, which is much faster and
        gives the same ``D`` up to a small fraction of ``D_std``.

        The observed cluster is the one passed to :py:meth:`calibrate`.

        :param params: Dictionary with the values of the fundamental parameters used
            to generate the synthetic cluster, see :py:meth:`generate`
        :type params: dict[str, float]
        :param h_mag: Width of the Gaussian kernel in magnitude, defaults to ``0.3``
        :type h_mag: float
        :param h_col: Width of the Gaussian kernel in color, defaults to ``0.075``
        :type h_col: float

        :raises ValueError: If :py:meth:`calibrate` was not called, if the observed
            cluster has fewer than two stars with valid photometry, or in the same
            cases as :py:meth:`generate`

        :return: ``(D, D_std)``, where ``D_std`` is the approximate standard
            deviation of ``D`` for a synthetic cluster drawn from the same
            distribution as the observed one (it decreases with the number of stars).
            If the synthetic cluster has fewer than two stars, ``D=1`` and ``D_std=0``
        :rtype: tuple[float, float]
        """
        obs = self._check_calibrated()

        from .modules.synthetic import likelihood as lh

        # Same synthetic cluster as generate(params); the isochrone (with DR=0) is
        # used to locate the turn-off point
        # isochrone, synth_clust = self._get_model(params, isoch_flag=True)
        _, synth_clust = self._get_model(params, isoch_flag=False)

        D, D_std = lh.cmd_overlap(obs.mag, obs.colors[0], synth_clust, h_mag, h_col)

        # Turn-off check: the observed and synthetic stars inside a box around the
        # turn-off are compared with the cmd_overlap() statistic. The check is skipped
        # if the box contains few stars.

        # if synth_clust.size == 0:
        #     return D, D_std

        # # Turn-off point of the isochrone up to the maximum mass of the synthetic
        # # cluster, as in cluster_properties()
        # max_mass = synth_clust[self._engine.m_ini_idx].max()
        # if (isochrone[self._engine.m_ini_idx] < max_mass).sum() < 2:
        #     return D, D_std
        # to = sst.turn_off_points([isochrone[[0, 1, self._engine.m_ini_idx]]], [max_mass], 0)
        # to_col, to_mag = to["to_col"][0], to["to_mag"][0]

        # lh.cmd_overlap_turnoff_box(
        #     self.cluster_mag,
        #     self.cluster_colors[0],
        #     synth_clust[0],
        #     synth_clust[1],
        #     to_mag,
        #     to_col,
        #     to_box_mag,
        #     to_box_col,
        # )

        return D, D_std

    def cluster_properties(
        self,
        model: dict[str, float],
        model_std: dict[str, float],
        N_models: int = 200,
        color_idx: int = 0,
        isochs_young: Isochrones | None = None,
        turn_off: Sequence[float] | None = None,
        radec_c: tuple[float, float] | None = None,
        rho_amb: float | None = None,
        galaxy_potential: dict[str, float] | None = None,
        disruption: dict[str, float] | None = None,
        return_arrays: bool = False,
    ) -> dict:
        r"""Estimate the per-star and the cluster properties of the observed cluster,
        from synthetic clusters sampled around a given solution.

        ``N_models`` sets of parameters are sampled from a normal distribution with
        mean and STDDEV given by ``model`` and ``model_std``, truncated to the valid
        range of each parameter: the isochrones' ranges for ``met, loga``, ``>=0`` for
        ``alpha, beta, Av, DR, dm``, and ``>0`` for ``Rv``. A synthetic cluster is
        generated for each set, and all the properties below are estimated from the
        same synthetic clusters.

        Per-star properties (``stars`` key of the returned dictionary), one array per
        property with a value for each observed star:

        - Primary mass (``mass_single``), and secondary mass (``mass_binar``) if
          binary systems were generated
        - Mass assuming the star is a blue straggler (``mass_bss``), if
          ``isochs_young`` is given
        - Probabilities of being a single system, a binary system, or a blue
          straggler (``single_prob, binar_prob, bss_prob``; these add up to ``1``)

        Cluster properties (``cluster`` key of the returned dictionary):

        - ``binar_fr``: total binary fraction. The number of binary systems among the
          observed stars follows a Poisson-binomial distribution given by their
          ``binar_prob`` values, which is computed exactly
        - ``binar_fr_models``: total binary fraction of each sampled synthetic
          cluster, i.e.: the fraction of the observed stars matched to a binary
          system (weighted by ``1-bss_prob``). Unlike ``binar_fr``, its spread
          includes the uncertainty in the fundamental parameters
        - ``N_bss``: number of blue stragglers, from the Poisson-binomial distribution
          given by the ``bss_prob`` values
        - ``M_init, M_actual, M_obs, M_phot, M_evol, M_dyn``: initial, actual,
          observed, photometric, evolutionary and dynamical masses, where
          ``M_actual = M_obs + M_phot`` and ``M_init = M_actual + M_evol + M_dyn``

        Each cluster property is given as its ``(median, 16th, 84th)`` percentiles.
        With ``return_arrays=True`` the per-model values of ``binar_fr_models`` and of
        the masses are returned instead.

        The BSS region of each sampled synthetic cluster contains the stars brighter
        than its turn-off point and bluer than its turn-off color shifted by ``-0.05``
        (to detach the region from the isochrone). The turn-off point is located on
        the isochrone moved to the absorption ``Av + f_5 * DR``, where ``f_5`` is the
        5th percentile of the differential reddening fractions of the
        ``DR_distribution`` model, i.e.: at the blue edge of the bulk of the stars
        instead of the bluest edge (``DR=0``). This color limit is capped at the
        bluest edge (5th percentile) of the observed colors around the turn-off
        magnitude, so that the region never extends into the observed main sequence.
        If ``turn_off`` is given, the turn-off points are instead sampled from a
        normal distribution around it (one per synthetic cluster), and used as given
        (no color shift or cap). The turn-off points used (before the color shift)
        are stored in the ``turn_off_points`` attribute, a dictionary with the keys
        ``color_idx, to_col, to_mag, manual`` (one ``to_col, to_mag`` value per
        synthetic cluster).

        The initial, evolutionary and dynamical masses require the ambient density
        ``rho_amb``. If it is not given, it is estimated using the cluster's position
        ``radec_c`` and a model for the Galaxy's potential. If ``radec_c`` is not given
        either, the median ``(ra, dec)`` of the observed cluster stored by
        :py:meth:`calibrate` are used. If no coordinates are available these three
        masses are ``nan`` and a warning is shown.

        Stars with no valid photometry are assigned ``nan`` values and excluded from
        the cluster properties.

        :param model: Dictionary with the values for the fundamental parameters.
            Parameters not included take the default values of the object
        :type model: dict[str, float]
        :param model_std: Dictionary with the standard deviations for the fundamental
            parameters in the ``model`` argument. Parameters not included are fixed
        :type model_std: dict[str, float]
        :param N_models: Number of sampled models
        :type N_models: int
        :param color_idx: Index of the color used to locate the turn-off point (and
            the BSS region). If two colors are defined, set to ``0`` to use the first
            color, or ``1`` to use the second color. Defaults to ``0``.
        :type color_idx: int
        :param isochs_young: Isochrones object (same photometric system) that contains
            isochrones with ``log(age)<=7``. These are used to estimate the
            ``mass_bss`` values. If ``None``, these values are ``nan``
        :type isochs_young: Isochrones | None
        :param turn_off: Manual turn-off point that defines the BSS region, given as
            ``(color, magnitude, stddev_color, stddev_magnitude)`` in the color
            selected by ``color_idx``. If ``None``, the turn-off point of each
            sampled synthetic cluster is estimated from its isochrone
        :type turn_off: Sequence[float] | None
        :param radec_c: Center coordinates in ``(RA, DEC)``
        :type radec_c: tuple[float, float] | None
        :param rho_amb: Ambient density [:math:`M_{\odot}\,pc^{-3}`]. If ``None``, it
            is estimated using the cluster's position and a model for the Galaxy's
            potential
        :type rho_amb: float | None
        :param galaxy_potential: Parameters of the Galaxy's potential used to estimate
            the ambient density. Keys not given take their default values:
            ``M_B=2.5e10`` [:math:`M_{\odot}`] and ``r_B=0.5e3`` [pc] (bulge mass and
            characteristic radius), ``M_D=7.5e10`` [:math:`M_{\odot}`], ``a=5.4e3``
            [pc] and ``b=0.3e3`` [pc] (disc mass, scale radius and scaleheight), from
            `Haghi et al. 2015 <https://doi.org/10.1093/mnras/stv827>`__, Table 1;
            ``M_s=1.87e11`` [:math:`M_{\odot}`] and ``r_s=15.19e3`` [pc] (dark matter
            halo mass and scale radius), from `Sanderson et al. 2017
            <https://iopscience.iop.org/article/10.3847/1538-4357/aa5eb4>`__, Table 1
        :type galaxy_potential: dict[str, float] | None
        :param disruption: Parameters of the cluster's disruption. Keys not given take
            their default values: ``C_env=810e6`` [yr] (constant related to the
            disruption time, from `Lamers, Gieles & Zwart 2005
            <https://www.aanda.org/articles/aa/abs/2005/01/aa1476/aa1476.html>`__),
            ``gamma=0.62`` (constant related to the disruption time, from
            `Boutloukos & Lamers 2003
            <https://academic.oup.com/mnras/article/338/3/717/1458721>`__), and
            ``epsilon=0.08`` (eccentricity of the orbit, from `Angelo et al. (2023)
            <https://doi.org/10.1093/mnras/stad1038>`__)
        :type disruption: dict[str, float] | None
        :param return_arrays: If ``True``, return the per-model values of
            ``binar_fr_models`` and of the masses instead of their percentiles
        :type return_arrays: bool

        :raises ValueError:
            - If the `calibrate()` method was not previously called
            - If ``turn_off`` does not contain four values or has negative STDDEVs,
              or ``galaxy_potential`` or ``disruption`` contain unknown keys
            - If ``rho_amb`` is not a finite value larger than zero
            - If ``isochs_young`` contains no isochrones with ``log(age)<=7``
            - If all the synthetic arrays generated are empty

        :return: Dictionary with the keys ``stars`` (dictionary of per-star arrays:
            ``star_id`` if the :py:class:`Cluster` object contains it, ``mass_single,
            mass_single_std, mass_binar, mass_binar_std, mass_bss, mass_bss_std,
            single_prob, binar_prob, bss_prob``) and ``cluster`` (dictionary with the
            keys ``binar_fr, binar_fr_models, N_bss, M_init, M_actual, M_obs, M_phot,
            M_evol, M_dyn``)
        :rtype: dict
        """
        obs = self._check_calibrated()

        from .modules.synthetic import cluster_mass as cm
        from .modules.synthetic import precompute as pre
        from .modules.synthetic import stellar_stats as sst

        # Check the arguments before generating the synthetic clusters
        if turn_off is not None:
            if len(turn_off) != 4:
                raise ValueError(
                    "'turn_off' must contain four values: (color, magnitude, "
                    + "stddev_color, stddev_magnitude)"
                )
            if turn_off[2] < 0 or turn_off[3] < 0:
                raise ValueError("The STDDEVs in 'turn_off' must be >= 0")
        if rho_amb is not None and not (np.isfinite(rho_amb) and rho_amb > 0):
            raise ValueError(f"'rho_amb' must be a finite value > 0, got {rho_amb}")
        gal_pot = cm.update_pars(
            "galaxy_potential", cm.GALAXY_POTENTIAL, galaxy_potential
        )
        disrup = cm.update_pars("disruption", cm.DISRUPTION, disruption)
        if isochs_young is not None and isochs_young.amin > 7:
            raise ValueError("No isochrones with log(ages)<=7 found in 'isochs_young'")

        # Fill in the default values. New dictionaries are created so that the ones
        # passed by the user are not modified
        model = self.def_params | model
        model_std = {k: 0 for k in self.def_params} | model_std
        self._vp("\nGenerate synthetic models", 1)
        self._vp("par            : mean +/ STDDEV", 2)
        for k in ("met", "loga", "alpha", "beta", "Rv", "DR", "Av", "dm"):
            self._vp(f"{k:<15}: {round(model[k], 4)} +/- {round(model_std[k], 4)}", 2)

        # The isochrone with DR=0 is the bluest edge of the synthetic cluster. The
        # turn-off points are located on the isochrone moved to the absorption
        # 'Av + f_to * DR' instead (the blue edge of the bulk of the stars), where
        # 'f_to' is a percentile of the differential reddening fractions. Not used
        # if the turn-off point is given
        #
        # Percentile of the differential reddening fractions 'f' (A = Av + f * DR) at
        # which the turn-off points are located. The DR=0 isochrone is the bluest edge
        # of the synthetic cluster, which for large DR values lies far to the blue of
        # the bulk of the stars (e.g.: 'foreground' model, where
        # f ~ N(0.5, 1/6)) and the BSS region misses most BSSs
        TO_DR_PERCENTILE = 5
        f_to = 0.0
        if turn_off is None:
            f_to = float(
                np.percentile(
                    self._engine.rand_floats["dr_frac"][self.DR_distribution],
                    TO_DR_PERCENTILE,
                )
            )

        # Sample the models and generate their synthetic clusters and isochrones.
        # Models that result in empty synthetic clusters are discarded
        models, synthcls, isoch_tracks, max_masses = [], [], [], []
        mass_infos = []
        for smodel in pre.sample_params(
            self.rng, model, model_std, self._engine.limits, N_models
        ):
            mass_info = {}
            isoch_arr, synth_arr = self._get_model(
                smodel, isoch_flag=True, isoch_dr_frac=f_to, mass_info=mass_info
            )
            if not synth_arr.any():
                continue
            models.append(smodel)
            synthcls.append(synth_arr)
            mass_infos.append(mass_info)
            # Magnitude, color and initial mass of the isochrone, and maximum mass of
            # the synthetic cluster, used to locate the turn-off point
            isoch_tracks.append(isoch_arr[[0, color_idx + 1, self._engine.m_ini_idx]])
            max_masses.append(synth_arr[self._engine.m_ini_idx].max())

        N_models_non_empty = len(models)
        if N_models_non_empty == 0:
            raise ValueError("All sampled models produced empty synthetic clusters.")
        if N_models_non_empty < N_models:
            warnings.warn(
                f"Empty synthetic arrays generated: {N_models - N_models_non_empty}"
            )
        self._vp(f"N_models       : {N_models_non_empty}", 1)

        if turn_off is None:
            turn_off_points = sst.turn_off_points(isoch_tracks, max_masses, color_idx)
        else:
            turn_off_points = sst.turn_off_points_manual(
                self.rng, turn_off, N_models_non_empty, color_idx
            )
        self.turn_off_points = turn_off_points

        # Used to estimate the BSS masses
        extrinsic_pars = {
            "values": {k: model[k] for k in ("Av", "dm", "Rv")},
            "stddevs": {k: model_std[k] for k in ("Av", "dm", "Rv")},
            "limits": self._engine.limits,
        }

        stars, m2_obs, valid = sst.stellar_parameters(
            obs.mag,
            obs.e_mag,
            obs.colors,
            obs.e_colors,
            self._engine.m_ini_idx,
            synthcls,
            turn_off_points,
            isochs_young,
            extrinsic_pars,
            self.ext_law,
            self._engine.ext_coefs,
            self.rng,
        )
        # Add ID column if it exists
        if obs.star_id is not None:
            stars = {"star_id": obs.star_id, **stars}
        N_nan = int((~valid).sum())
        if N_nan > 0:
            self._vp("", 0)
            warnings.warn(
                f"\nN={N_nan} stars found with no valid photometric data. "
                + "These will be assigned 'nan' values\nfor masses and "
                + "probabilities"
            )

        # Binary fraction and number of BSSs, from the stars with valid photometry
        N_valid = int(valid.sum())
        binar_p, bss_p = stars["binar_prob"][valid], stars["bss_prob"][valid]
        binar_fr = sst.poisson_binomial_percentiles(binar_p) / N_valid
        N_bss = sst.poisson_binomial_percentiles(bss_p)
        # Binary fraction of each sampled model. The binary systems are weighted by
        # (1 - bss_prob), as the 'binar_prob' values
        binar_fr_models = (np.isfinite(m2_obs[:, valid]) * (1.0 - bss_p)).sum(
            axis=1
        ) / N_valid

        # Use the position of the observed cluster if no ambient density or
        # coordinates are given
        if (
            rho_amb is None
            and radec_c is None
            and obs.ra is not None
            and obs.dec is not None
        ):
            # The RA values are measured relative to their circular mean before
            # taking the median, so that clusters crossing RA=0 are handled
            ra = np.asarray(obs.ra, dtype=float)
            ra_rad = np.deg2rad(ra)
            ra_ref = np.rad2deg(
                np.arctan2(np.nanmean(np.sin(ra_rad)), np.nanmean(np.cos(ra_rad)))
            )
            ra_off = (ra - ra_ref + 180.0) % 360.0 - 180.0
            radec_c = (
                float((ra_ref + np.nanmedian(ra_off)) % 360.0),
                float(np.nanmedian(obs.dec)),
            )
            self._vp(
                "\nUsing ambient density estimated from observed cluster (ra, dec)", 1
            )
        masses = cm.cluster_masses(
            models,
            synthcls,
            mass_infos,
            self._engine.m_ini_idx,
            self.isochs.z_to_FeH,
            radec_c,
            rho_amb,
            gal_pot,
            disrup,
        )

        per_model = {"binar_fr_models": binar_fr_models, **masses}
        if not return_arrays:
            per_model = {
                k: tuple(float(_) for _ in np.percentile(v, (50, 16, 84)))
                for k, v in per_model.items()
            }
        cluster = {
            "binar_fr": tuple(float(_) for _ in binar_fr),
            "N_bss": tuple(int(_) for _ in N_bss),
            **per_model,
        }

        self._vp("\nCluster properties estimated", 1)

        return {"stars": stars, "cluster": cluster}

    def get_isochrone(
        self,
        fit_params: dict,
        color_idx: int = 0,
        full_track: bool = False,
        N_stars: int = 1000,
    ) -> np.ndarray:
        """Generate an isochrone for plotting.

        The isochrone is generated using the fundamental parameter values
        given in the ``fit_params`` dictionary.

        :param fit_params: Dictionary with the values for the fundamental parameters.
            Parameters not included take the default values defined in the
            ``def_params`` argument of the :py:class:`Synthetic` object
        :type fit_params: dict
        :param color_idx: Index of the color to plot. If ``0`` (default), plot the
            first color. If ``1`` plot the second color
        :type color_idx: int
        :param full_track: If ``True``, return the full isochrone (ie: do not apply
            the max mass filter)
        :type full_track: bool
        :param N_stars: Number of synthetic stars to generate
        :type N_stars: int

        :return: Array with the isochrone data to plot
        :rtype: np.ndarray
        """
        # Generate displaced isochrone and full synthetic CMD array
        isochrone, full_arr = self._get_model(fit_params, N_stars, isoch_flag=True)

        if full_track:
            # Return the full track without applying the max mass filter
            return np.array(isochrone)[[0, color_idx + 1]]

        if full_arr.size == 0:
            # No synthetic stars were generated (e.g.: all the isochrone is below
            # the maximum magnitude), the max mass filter can not be applied
            warnings.warn(
                "The synthetic cluster for these parameters is empty, returning "
                + "the full isochrone",
                stacklevel=2,
            )
            return np.array(isochrone)[[0, color_idx + 1]]

        # Extract magnitude and mass from the isochrone and the full synthetic cluster
        isoch_mag, isoch_mass = isochrone[0], isochrone[self._engine.m_ini_idx]
        mag_full, mass_full = full_arr[0], full_arr[self._engine.m_ini_idx]
        full_mag_min = mag_full.min()

        # Select the max mass for the isochrone that results in a minimum difference
        # between the minimum magnitude of the isochrone and the minimum magnitude of
        # the full synthetic cluster
        best_diff = np.inf
        best_max_mass_isoch = None
        for p in np.arange(95, 100.5, 0.5):
            max_mass_isoch = np.percentile(mass_full, p)
            msk = isoch_mass < max_mass_isoch
            if not np.any(msk):
                continue
            isoch_mag_min = isoch_mag[msk].min()
            diff = abs(isoch_mag_min - full_mag_min)
            if diff < best_diff:
                best_diff = diff
                best_max_mass_isoch = max_mass_isoch

        # Apply max mass filter to isochrone
        if best_max_mass_isoch is None:
            # No isochrone points below any of the mass percentiles: keep them all
            msk = np.full(isoch_mass.size, True)
        else:
            msk = isoch_mass < best_max_mass_isoch
        # Generate proper array for plotting
        isochrone = np.array([isoch_mag, isochrone[color_idx + 1]])[:, msk]

        return isochrone

    def _check_calibrated(self) -> ObservedData:
        """Data of the observed cluster stored by :py:meth:`calibrate`

        :raises ValueError: If :py:meth:`calibrate` was not called

        :return: Observed cluster data
        :rtype: ObservedData
        """
        if self._obs is None:
            raise ValueError(
                "This method requires running the 'calibrate()' method first\n"
                + "with an observed 'Cluster' object"
            )
        return self._obs

    def _calib_values(self, N_stars: int) -> tuple[float, list[np.ndarray], int]:
        """Maximum magnitude, error distribution and number of stars of the
        synthetic clusters: those of the observed cluster if :py:meth:`calibrate` was
        called, else no magnitude cut, no errors and ``N_stars`` stars

        :param N_stars: Number of synthetic stars if the object was not calibrated
        :type N_stars: int

        :return: Maximum magnitude, error distribution and number of stars
        :rtype: tuple[float, list[np.ndarray], int]
        """
        if self._obs is None:
            return np.inf, [], int(N_stars)
        return self._obs.calib_values()

    def _get_model(
        self,
        params: dict | None,
        N_stars: int = 100,
        isoch_flag: bool = False,
        show_errors: bool = False,
        isoch_dr_frac: float = 0.0,
        mass_info: dict | None = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Generate a synthetic cluster (and its isochrone) with
        :py:func:`generate.get_model`, from the data stored in the object

        :param params: Dictionary with the values for the fundamental parameters. If
            ``None``, the default values are used
        :type params: dict | None
        :param N_stars: Number of synthetic stars if the object was not calibrated
        :type N_stars: int
        :param isoch_flag: If ``True``, also return the isochrone (with ``DR=0``)
        :type isoch_flag: bool
        :param show_errors: If ``True``, print the reason when an empty synthetic
            cluster is returned
        :type show_errors: bool
        :param isoch_dr_frac: Fraction ``f`` of the differential reddening applied to
            the returned isochrone (absorption ``Av + f * DR``; ``0`` for ``DR=0``).
            It does not affect the synthetic cluster
        :type isoch_dr_frac: float
        :param mass_info: If a dictionary is given, it is filled with the data
            required to estimate the total mass of the simulated cluster, see
            :py:func:`generate.get_model`
        :type mass_info: dict | None

        :raises ValueError: If the parameters are not valid, or if more than one value
            is given for a parameter

        :return: Isochrone (empty if ``isoch_flag=False``) and synthetic cluster
        :rtype: tuple[np.ndarray, np.ndarray]
        """
        P = self._engine.params_array(params)
        if P.shape[0] != 1:
            raise ValueError(
                "A single value is required for each parameter (use 'batch=True' to "
                + "generate several synthetic clusters)"
            )
        return self._engine.get_model(
            P[0],
            self.DR_distribution,
            self._calib_values(N_stars),
            show_errors=show_errors,
            isoch_flag=isoch_flag,
            isoch_dr_frac=isoch_dr_frac,
            mass_info=mass_info,
        )

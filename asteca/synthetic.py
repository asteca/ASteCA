import warnings

import numpy as np

from .cluster import Cluster
from .isochrones import Isochrones


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
    :param def_params: Default values for all the required fundamental parameters. If
        ``None`` is given, the following values are used::

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
        ``Av + DR`` the maximum absorption
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

        if def_params is None:
            def_params = {
                "met": 0.0152,
                "loga": 8.0,
                "alpha": 0.09,
                "beta": 0.94,
                "Rv": 3.1,
                "DR": 0.0,
                "Av": 0.2,
                "dm": 9.0,
            }

        self.isochs = isochs
        self.def_params = def_params
        self.IMF_name = IMF_name
        self.max_mass = max_mass
        self.gamma = gamma
        self.ext_law = ext_law
        self.DR_distribution = DR_distribution
        self.seed = seed
        self.verbose = verbose

        # Set seed
        if self.seed is None:
            self.seed = np.random.randint(0, 100_000)
        self.rng = np.random.default_rng(self.seed)

        # Check gamma distribution
        gammas = ("D&K", "fisher_stepped", "fisher_peaked", "raghavan")
        if isinstance(self.gamma, str) and self.gamma not in gammas:
            raise ValueError(
                f"gamma '{self.gamma}' not recognized. Should be one of {gammas}"
            )

        # Check differential reddening function
        DR_funcs = ("foreground", "internal")
        if self.DR_distribution not in DR_funcs:
            raise ValueError(
                f"Differential reddening function '{self.DR_distribution}' not "
                + f"recognized. Should be one of {DR_funcs}"
            )

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

        from .modules import synth_cluster_priv as scp

        # Sample the selected IMF
        Nmets, Nages = self.isochs.theor_tracks.shape[:2]
        self._vp("(Sampling IMF)", 1)
        self.st_dist_mass = scp.sample_imf(
            self.rng, self.IMF_name, self.max_mass, Nmets, Nages
        )

        # Add binary systems
        self.theor_tracks = scp.add_binarity(
            self.rng,
            self.gamma,
            self.isochs.color,
            self.isochs.color2,
            self.isochs.theor_tracks,
            self.isochs.color_filters,
        )

        # Data used by the `generate()` method
        self.m_ini_idx = 2  # (0->mag, 1->color, 2->mass_ini)
        if self.isochs.color2_effl is not None:
            self.m_ini_idx = 3  # (0->mag, 1->color, 2->color2, 3->mass_ini)
        # Set tolerances for the ranges of the isochrones in `generate()`
        self.age_tol = 0.01
        self.met_tol = 0.001  # z
        if self.isochs.z_to_FeH is not None:
            self.met_tol = 0.01  # FeH

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
                    f"\nExtinction law '{self.ext_law}' does not apply to the second color.\n"
                    + "Use the 'CCMO' extinction law if you need to process a second color"
                )

            # The Rv parameter is not used by this law
            if self.def_params.get("Rv", 3.1) != 3.1:
                warnings.warn(
                    f"\nExtinction law '{self.ext_law}' does not use the 'Rv' parameter: "
                    + f"the default value Rv={self.def_params['Rv']} will have no effect."
                )

        # Get extinction coefficients for these filters
        self.ext_coefs = []  # It is important to pass an empty list because the
        # function `extinction()` checks its length later on
        if self.ext_law == "CCMO":
            if self.isochs.magnitude_effl is None or self.isochs.color_effl is None:
                raise ValueError(
                    f"Extinction law '{self.ext_law}' requires effective lambda\n"
                    + "values for the magnitude and first color."
                )
            self.ext_coefs = scp.ccmo_ext_coeffs(
                self.isochs.magnitude_effl,
                self.isochs.color_effl,
                self.isochs.color2_effl,
            )

        # Generate random floats used by `generate()`
        self.rand_floats = scp.randVals(
            self.rng, self.theor_tracks, self.st_dist_mass, self.DR_distribution
        )

        # Store for internal usage
        self.met_age_dict = self.isochs.met_age_dict

        # Pre-computed search indexes used by mass_interp()
        self.mass_idx = scp.mass_idx_lookup(
            self.theor_tracks, self.st_dist_mass, self.m_ini_idx
        )

        # Check that the ranges are respected
        for par in ("met", "loga"):
            if par in def_params:
                pmin, pmax = min(self.met_age_dict[par]), max(self.met_age_dict[par])
                if self.def_params[par] < pmin or self.def_params[par] > pmax:
                    warnings.warn(
                        f"Parameter {par}={self.def_params[par]} is out of range: [{pmin} - {pmax}]"
                    )

        required = {
            "met",
            "loga",
            "alpha",
            "beta",
            "Rv",
            "DR",
            "Av",
            "dm",
        }
        if not required.issubset(def_params.keys()):
            missing = required - set(def_params.keys())
            warnings.warn(
                f"Missing required keys in 'def_params': {missing}.\nThese must be all "
                "present in 'params' when calling 'synthcl.generate(params)'"
            )

        self._vp(f"Default params : {self.def_params}", 1)
        self._vp(f"Extinction law : {self.ext_law}", 1)
        self._vp(f"Diff reddening : {self.DR_distribution}", 1)
        self._vp(f"IMF            : {self.IMF_name}", 1)
        self._vp(f"Max init mass  : {self.max_mass}", 1)
        self._vp(f"Gamma dist     : {self.gamma}", 1)
        self._vp(f"Random seed    : {self.seed}", 1)
        self._vp("Synthetic clusters object generated")

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
        from .modules import likelihood_priv as lpriv
        from .modules import synth_cluster_priv as scp

        if cluster.mag is None:
            raise ValueError(
                "Attribute 'mag' is required to calibrate 'Cluster' object."
            )
        if cluster.e_mag is None:
            raise ValueError(
                "Attribute 'e_mag' is required to calibrate 'Cluster' object."
            )
        if cluster.color is None:
            raise ValueError(
                "Attribute 'color' is required to calibrate 'Cluster' object."
            )
        if cluster.e_color is None:
            raise ValueError(
                "Attribute 'e_color' is required to calibrate 'Cluster' object."
            )

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
        if cluster.mag.size * f_stars_mass > self.max_mass:
            warnings.warn(
                (
                    f"\nThe number of stars in this `Cluster` object ({cluster.mag.size}) is large.\n"
                    "This may lead to synthetic clusters with insufficient stars to\n"
                    "properly reproduce the observations.\n"
                    "Consider increasing the 'max_mass' parameter in the Synthetic cluster \n"
                    f"object (current value: {self.max_mass}). As a rule of thumb, set 'max_mass' \n"
                    f"to at least {f_stars_mass} times the number of observed stars."
                ),
                UserWarning,
            )

        # Calibrate parameters using the observed cluster
        self.max_mag_syn_obs = float(np.nanmax(cluster.mag))
        self.N_stars_obs = cluster.N_stars - cluster.N_nan_phot
        self.err_dist_obs = scp.error_distribution(
            cluster.mag,
            cluster.e_mag,
            cluster.e_color,
            cluster.e_color2,
            self.rand_floats["norm"][1],
        )

        # Used by the `cluster_properties()` method
        self.star_id = cluster.star_id
        self.cluster_mag = cluster.mag
        self.cluster_colors = [cluster.color, cluster.color2]
        # Used by get_stellar_masses()
        self.cluster_e_mag = cluster.e_mag
        self.cluster_e_colors = [cluster.e_color, cluster.e_color2]

        # Used by the `cluster_properties()` method
        self.cluster_ra = cluster.ra
        self.cluster_dec = cluster.dec

        # Data used by the likelihood (asteca.Sampler)
        colors = [cluster.color]
        if cluster.color2 is not None:
            colors += [cluster.color2]
        self._lkl = lpriv.Tremmel(cluster.mag, colors)

        self._vp("\nCalibrated observed cluster")
        self._vp(f"N_stars_obs    : {self.N_stars_obs}", 1)
        self._vp(f"Max magnitude  : {self.max_mag_syn_obs:.2f}", 1)
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

        A ``ValueError`` is raised if ``met`` or ``loga`` are outside of the ranges of
        the loaded isochrones, if any of ``alpha, beta, Av, DR, dm`` is negative, or if
        ``Rv<=0`` (for any of the models if ``batch=True``).

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
            from .modules import synth_batch as sb

            P = self._batch_params(params)
            if not sb.available(self.ext_law):
                return [
                    self.generate(dict(zip(sb.FULL_PARAMS, p)), N_stars=N_stars)
                    for p in P
                ]

            out = sb._generate_batch(*self._batch_args(P, N_stars))
            # An empty synthetic cluster is returned as an empty 1D array, as with
            # batch=False
            return [a if a.shape[1] else np.array([]) for a in out]

        from .modules import synth_model as sm

        if params is None:
            params = {}

        try:
            # Extract from calibrated observed cluster
            max_mag_syn = self.max_mag_syn_obs
            err_dist_synth = self.err_dist_obs
            N_synth_stars = self.N_stars_obs
        except AttributeError:
            # If no calibration was done, use the provided 'N_stars' value (or the
            # default value if not given)
            max_mag_syn = np.inf
            err_dist_synth = []
            N_synth_stars = int(N_stars)

        _, synth_arr = sm.get_model(
            params,
            self.met_age_dict,
            self.met_tol,
            self.age_tol,
            self.def_params,
            self.m_ini_idx,
            self.theor_tracks,
            self.ext_law,
            self.ext_coefs,
            self.rand_floats,
            self.st_dist_mass,
            self.mass_idx,
            max_mag_syn,
            err_dist_synth,
            N_synth_stars,
            show_errors,
        )

        return synth_arr

    def cmd_overlap(  # noqa: DOC502, RUF100
        self,
        params: dict[str, float],
        h_mag: float = 0.3,
        h_col: float = 0.075,
        to_box_mag: tuple[float, float] = (-1.0, 0.5),
        to_box_col: tuple[float, float] = (-0.2, 0.8),
    ) -> tuple[float, float]:
        """Goodness of fit of the synthetic cluster generated with the ``params``
        parameters, that can be compared across observed clusters with any number of
        stars. The method also checks that the turn-off region of the observed cluster
        is reproduced by the synthetic cluster, and shows a warning if it is not.

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

        Turn-off check: the turn-off point of the isochrone of the synthetic cluster
        is located (as in :py:meth:`cluster_properties`) and a box is defined around it,
        from ``to_box_mag[0]`` to ``to_box_mag[1]`` mag in magnitude and from
        ``to_box_col[0]`` to ``to_box_col[1]`` mag in color, relative to the
        turn-off (by default from 1 mag brighter to 0.5 mag fainter, including the
        subgiant and lower red giant branches). The distributions of the observed and
        synthetic stars inside this box are compared with the same ``D`` statistic,
        using narrower kernels (``0.1`` mag in magnitude and ``0.025`` mag in color)
        to resolve the turn-off. A warning is shown if this ``D`` is larger than
        four times its standard deviation, i.e. if the turn-off region is clearly
        not reproduced. Since the faint main sequence contains most of the stars, it
        dominates the likelihood (and ``D``), and the fit can sacrifice the turn-off
        region, which constrains the age and metallicity. This happens e.g. when the
        observed luminosity function is flatter than the one of the synthetic
        cluster (dynamically evolved clusters that lost low-mass stars, or an
        incomplete faint end), see :ref:`cmd_overlap`. The check requires at least
        ``30`` observed and ``30`` synthetic stars in the box, i.e. it is not
        performed for most young clusters, whose turn-off is sparsely populated.

        :param params: Dictionary with the values of the fundamental parameters used
            to generate the synthetic cluster, see :py:meth:`generate`
        :type params: dict[str, float]
        :param h_mag: Width of the Gaussian kernel in magnitude, defaults to ``0.3``
        :type h_mag: float
        :param h_col: Width of the Gaussian kernel in color, defaults to ``0.075``
        :type h_col: float
        :param to_box_mag: Magnitude limits of the turn-off box, relative to the
            magnitude of the turn-off point, defaults to ``(-1.0, 0.5)``
        :type to_box_mag: tuple[float, float]
        :param to_box_col: Color limits of the turn-off box, relative to the color of
            the turn-off point, defaults to ``(-0.2, 0.8)``
        :type to_box_col: tuple[float, float]

        :raises ValueError: If :py:meth:`calibrate` was not called, if the observed
            cluster has fewer than two stars with valid photometry, or in the same
            cases as :py:meth:`generate`

        :return: ``(D, D_std)``, where ``D_std`` is the approximate standard
            deviation of ``D`` for a synthetic cluster drawn from the same
            distribution as the observed one (it decreases with the number of stars).
            If the synthetic cluster has fewer than two stars, ``D=1`` and ``D_std=0``
        :rtype: tuple[float, float]
        """
        self._check_calibrated()

        from .modules import likelihood_priv as lpriv
        from .modules import stellar_stats_funcs as ssf
        from .modules import synth_model as sm

        # Same synthetic cluster as generate(params); the isochrone (with DR=0) is
        # used to locate the turn-off point
        isochrone, synth_clust = sm.get_model(
            params,
            self.met_age_dict,
            self.met_tol,
            self.age_tol,
            self.def_params,
            self.m_ini_idx,
            self.theor_tracks,
            self.ext_law,
            self.ext_coefs,
            self.rand_floats,
            self.st_dist_mass,
            self.mass_idx,
            self.max_mag_syn_obs,
            self.err_dist_obs,
            self.N_stars_obs,
            isoch_flag=True,
        )

        D, D_std = lpriv.cmd_overlap(
            self.cluster_mag, self.cluster_colors[0], synth_clust, h_mag, h_col
        )
        if synth_clust.size == 0:
            return D, D_std

        # Turn-off point of the isochrone up to the maximum mass of the synthetic
        # cluster, as in cluster_properties()
        msk = isochrone[self.m_ini_idx] < synth_clust[self.m_ini_idx].max()
        if msk.sum() < 2:
            return D, D_std
        to_col, to_mag, _ = ssf.to_point_find(
            isochrone[1][msk], isochrone[0][msk], self.isochs.model, 0
        )

        # Turn-off check: the observed and synthetic stars inside a box around the
        # turn-off are compared with the cmd_overlap() statistic. The check is skipped
        # if the box contains few stars.
        lpriv.turnoff_box(
            self.cluster_mag,
            self.cluster_colors[0],
            synth_clust[0],
            synth_clust[1],
            to_mag,
            to_col,
            to_box_mag,
            to_box_col,
        )

        return D, D_std

    def _log_lkl_batch(  # noqa: DOC502, RUF100
        self,
        params: dict | None,
        temperature: float,
        parallel: bool,
        n_threads: int | None,
    ) -> np.ndarray:
        """Log-likelihoods of a batch of synthetic clusters, used by
        :py:class:`asteca.Sampler`.

        Same as ``self._lkl.log_lkl(self.generate(p), temperature)`` for each set of
        parameters ``p`` in ``params`` (defined as in :py:meth:`generate` with
        ``batch=True``), with identical results. With numba, the compiled functions in
        :py:mod:`synth_batch` generate each synthetic cluster and evaluate its
        likelihood in the same call (the synthetic clusters never go back to
        Python, which is where the gain comes from), optionally distributing the
        models among ``n_threads`` threads if ``parallel=True`` (if ``None``,
        numba's current setting is used). Without numba, ``parallel`` and
        ``n_threads`` are ignored. The ``temperature`` is assumed positive (it is
        validated by :py:class:`asteca.Sampler`).

        :param params: Dictionary with an array of values (or a single value) for each
            parameter, as in :py:meth:`generate` with ``batch=True``. If ``None``, the
            default parameters are used
        :type params: dict | None
        :param temperature: Likelihood temperature, assumed positive (validated by
            :py:class:`asteca.Sampler`)
        :type temperature: float
        :param parallel: Evaluate the models in parallel threads (requires numba)
        :type parallel: bool
        :param n_threads: Number of threads used if ``parallel=True``. If ``None``,
            numba's current setting is used
        :type n_threads: int | None

        :raises ValueError: In the same cases as :py:meth:`generate` for any of the
            models, if :py:meth:`calibrate` was not called, or if ``n_threads`` is
            used and it is not an integer between ``1`` and ``NUMBA_NUM_THREADS``

        :return: Array with the log-likelihood of each model, ``0`` is the best value
        :rtype: np.ndarray
        """
        from .modules import synth_batch as sb

        self._check_calibrated()

        lk = self._lkl
        P = self._batch_params(params)
        if not sb.available(self.ext_law):
            return np.array(
                [
                    lk.log_lkl(self.generate(dict(zip(sb.FULL_PARAMS, p))), temperature)
                    for p in P
                ]
            )

        args: tuple = (
            *self._batch_args(P),
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

        # 'n_threads' is only validated here, where it is used
        sb.check_n_threads(n_threads)
        # Set the number of threads for this call only
        prev_threads = numba.get_num_threads()
        numba.set_num_threads(int(n_threads))
        try:
            return sb._lnl_batch_parallel(*args)
        finally:
            numba.set_num_threads(prev_threads)

    def _check_calibrated(self) -> None:
        """Raise a ``ValueError`` if :py:meth:`calibrate` was not called"""
        if not (
            hasattr(self, "max_mag_syn_obs")
            and hasattr(self, "err_dist_obs")
            and hasattr(self, "N_stars_obs")
            and hasattr(self, "_lkl")
        ):
            raise ValueError(
                "This method requires running the 'calibrate()' method first\n"
                + "with an observed 'Cluster' object"
            )

    def _batch_params(self, params: dict | None) -> np.ndarray:
        """Full parameter array for the batch methods, see synth_batch.full_params()"""
        from .modules import synth_batch as sb

        return sb.full_params(
            {} if params is None else params,
            self.def_params,
            self.met_age_dict,
            self.met_tol,
            self.age_tol,
        )

    def _batch_args(self, P: np.ndarray, N_stars: int = 100) -> tuple:
        """Arguments of the compiled batch functions in synth_batch"""
        from .modules import synth_batch as sb

        # Data that does not change after the object is created, obtained once
        if getattr(self, "_batch_data", None) is None:
            self._batch_data = sb.batch_data(self)
        bd = self._batch_data

        # Same as in generate()
        try:
            max_mag_syn = self.max_mag_syn_obs
            err_dist_synth = self.err_dist_obs
            N_synth_stars = self.N_stars_obs
        except AttributeError:
            max_mag_syn = np.inf
            err_dist_synth = []
            N_synth_stars = int(N_stars)
        if len(err_dist_synth):
            err = np.array(err_dist_synth, dtype=float)
        else:
            err = np.empty((0, 0))

        return (
            P,
            bd["met_grid"],
            bd["age_grid"],
            self.theor_tracks,
            self.m_ini_idx,
            bd["st_list"],
            bd["midx_list"],
            bd["N_a"],
            bd["dr_frac"],
            bd["ext_kind"],
            bd["ccm"],
            bd["unif1"],
            float(max_mag_syn),
            err,
            int(N_synth_stars),
        )

    def cluster_properties(
        self,
        model: dict[str, float],
        model_std: dict[str, float],
        N_models: int = 200,
        color_idx: int = 0,
        isochs_young: Isochrones | None = None,
        bss_box: tuple[float, float, float, float] = (0.5, -0.1, -0.5, -0.05),
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

        The BSS region is defined by the turn-off point of each sampled synthetic
        cluster plus the offsets in ``bss_box``. Its color limits are capped at the
        bluest edge (5th percentile) of the observed colors around the turn-off
        magnitude, so that the region never extends into the observed main sequence.

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
        :param bss_box: Offsets applied to the turn-off point to define the two points
            of the BSS region: ``(mag_offset_1, col_offset_1, mag_offset_2,
            col_offset_2)``
        :type bss_box: tuple[float, float, float, float]
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
            - If ``bss_box`` does not contain four values, or ``galaxy_potential`` or
              ``disruption`` contain unknown keys
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
        self._check_calibrated()

        from .modules import cluster_mass_funcs as cmf
        from .modules import stellar_stats_funcs as ssf

        # Check the arguments before generating the synthetic clusters
        if len(bss_box) != 4:
            raise ValueError(
                "'bss_box' must contain four values: (mag_offset_1, col_offset_1, "
                + "mag_offset_2, col_offset_2)"
            )
        gal_pot = cmf.update_pars(
            "galaxy_potential", cmf.GALAXY_POTENTIAL, galaxy_potential
        )
        disrup = cmf.update_pars("disruption", cmf.DISRUPTION, disruption)
        if isochs_young is not None and isochs_young.amin > 7:
            raise ValueError("No isochrones with log(ages)<=7 found in 'isochs_young'")

        models, synthcls, turn_off_points, extrinsic_pars = self._sample_models(
            model, model_std, N_models, color_idx
        )

        stars, m2_obs, valid = self._stellar_parameters(
            synthcls, turn_off_points, extrinsic_pars, isochs_young, bss_box
        )

        # Binary fraction and number of BSSs, from the stars with valid photometry
        N_valid = int(valid.sum())
        binar_p, bss_p = stars["binar_prob"][valid], stars["bss_prob"][valid]
        binar_fr = ssf.poisson_binomial_percentiles(binar_p) / N_valid
        N_bss = ssf.poisson_binomial_percentiles(bss_p)
        # Binary fraction of each sampled model. The binary systems are weighted by
        # (1 - bss_prob), as the 'binar_prob' values
        binar_fr_models = (np.isfinite(m2_obs[:, valid]) * (1.0 - bss_p)).sum(
            axis=1
        ) / N_valid

        masses = self._cluster_masses(
            models, synthcls, radec_c, rho_amb, gal_pot, disrup
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

    def _sample_models(
        self,
        model: dict[str, float],
        model_std: dict[str, float],
        N_models: int,
        color_idx: int,
    ) -> tuple[list[dict], list[np.ndarray], dict, dict]:
        """Sample models around the given solution and generate their synthetic
        clusters, used by :py:meth:`cluster_properties`

        :param model: Dictionary with the values for the fundamental parameters
        :type model: dict[str, float]
        :param model_std: Dictionary with the standard deviations for the fundamental
            parameters in the ``model`` argument
        :type model_std: dict[str, float]
        :param N_models: Number of sampled models
        :type N_models: int
        :param color_idx: Index of the color used to locate the turn-off point
        :type color_idx: int

        :raises ValueError: If all the synthetic arrays generated are empty

        :return: The sampled models and their synthetic clusters (empty ones are
            discarded), the turn-off points of the synthetic clusters, and the values
            and STDDEVs of the extrinsic parameters
        :rtype: tuple[list[dict], list[np.ndarray], dict, dict]
        """
        from .modules import stellar_stats_funcs as ssf
        from .modules import synth_cluster_priv as scp
        from .modules import synth_model as sm

        # Update dictionaries with default values, if required. New dictionaries are
        # created so that the ones passed by the user are not modified
        model = self.def_params | model
        model_std = {k: 0 for k in self.def_params} | model_std

        pars_k = ("met", "loga", "alpha", "beta", "Rv", "DR", "Av", "dm")
        self._vp("\nGenerate synthetic models", 1)
        self._vp("par            : mean +/ STDDEV", 2)
        for k in pars_k:
            self._vp(
                f"{k:<15}: {round(model[k], 4)} +/- {round(model_std[k], 4)}",
                2,
            )

        # Generate the 'N_models' models via sampling a Gaussian centered on
        # 'model', with standard deviation given by 'model_std', truncated to the
        # valid range of each parameter: the isochrones' ranges for (met, loga),
        # and non-negative values for the rest
        bounds = {
            "met": (self.met_age_dict["met"][0], self.met_age_dict["met"][-1]),
            "loga": (self.met_age_dict["loga"][0], self.met_age_dict["loga"][-1]),
        }
        sampled_pars = {}
        for k, f_val in model.items():
            low, high = bounds.get(k, (scp.PARAMS_MIN.get(k, -np.inf), np.inf))
            sampled_pars[k] = scp.truncnorm_sample(
                self.rng, f_val, model_std[k], low, high, N_models
            )
        sampled_models = [
            dict(zip(sampled_pars, t)) for t in zip(*sampled_pars.values())
        ]

        # Generate the synthetic arrays using the above models
        models_no_empty, synthcls = [], []
        turn_off_points = {
            "color_idx": color_idx,
            "to_col": [],
            "to_mag": [],
        }
        # Initial minimum isochrone index for identifying the turn-off point
        fallback_min_idx = 10
        to_idx_lst = [fallback_min_idx]
        for smodel in sampled_models:
            # Get both the isochrone and the synthetic cluster
            isoch_arr, synth_arr = sm.get_model(
                smodel,
                self.met_age_dict,
                self.met_tol,
                self.age_tol,
                self.def_params,
                self.m_ini_idx,
                self.theor_tracks,
                self.ext_law,
                self.ext_coefs,
                self.rand_floats,
                self.st_dist_mass,
                self.mass_idx,
                self.max_mag_syn_obs,
                self.err_dist_obs,
                self.N_stars_obs,
                isoch_flag=True,
            )

            # Do not store models that result in empty arrays
            if not synth_arr.any():
                continue

            models_no_empty.append(smodel)
            synthcls.append(synth_arr)

            # Estimate turn-off point(s)
            # Extract maximum mass
            max_mass_isoch = synth_arr[self.m_ini_idx].max()
            # Apply max mass filter to isochrone
            msk = isoch_arr[self.m_ini_idx] < max_mass_isoch
            # Extract isochrone up to the max mass
            mag, color = np.array([isoch_arr[0], isoch_arr[color_idx + 1]])[:, msk]
            # Find turn-off point for the isochrone
            to_col, to_mag, to_idx = ssf.to_point_find(
                color, mag, self.isochs.model, fallback_min_idx
            )

            # If the turn-off point is below the minimum allowed, use the median of the
            # previous (non fallback) turn-off points
            if (to_idx < fallback_min_idx) and turn_off_points["to_col"]:
                to_col = np.median(turn_off_points["to_col"])
                to_mag = np.median(turn_off_points["to_mag"])

            to_idx_lst.append(to_idx)
            # Use the 10% percentile as the minimum fallback value for the turn-off index
            fallback_min_idx = int(np.percentile(to_idx_lst, 10))

            turn_off_points["to_col"].append(to_col)
            turn_off_points["to_mag"].append(to_mag)

        N_models_non_empty = len(models_no_empty)
        if N_models_non_empty == 0:
            raise ValueError("All sampled models produced empty synthetic clusters.")
        if N_models_non_empty < N_models:
            warnings.warn(
                f"Empty synthetic arrays generated: {N_models - N_models_non_empty}"
            )

        # Used to estimate the BSS masses
        extrinsic_pars = {
            "values": {"Av": model["Av"], "dm": model["dm"], "Rv": model["Rv"]},
            "stddevs": {
                "Av": model_std["Av"],
                "dm": model_std["dm"],
                "Rv": model_std["Rv"],
            },
        }

        self._vp(f"N_models       : {N_models_non_empty}", 1)

        return models_no_empty, synthcls, turn_off_points, extrinsic_pars

    def _stellar_parameters(
        self,
        synthcls: list[np.ndarray],
        turn_off_points: dict,
        extrinsic_pars: dict,
        isochs_young: Isochrones | None,
        bss_box: tuple[float, float, float, float],
    ) -> tuple[dict, np.ndarray, np.ndarray]:
        """Per-star masses and probabilities, used by :py:meth:`cluster_properties`

        :param synthcls: Synthetic clusters of the sampled models
        :type synthcls: list[np.ndarray]
        :param turn_off_points: Turn-off points of the synthetic clusters
        :type turn_off_points: dict
        :param extrinsic_pars: Values and STDDEVs of the extrinsic parameters
        :type extrinsic_pars: dict
        :param isochs_young: Isochrones with ``log(age)<=7``, or ``None``
        :type isochs_young: Isochrones | None
        :param bss_box: Offsets that define the BSS region
        :type bss_box: tuple[float, float, float, float]

        :return: Dictionary of per-star arrays, the secondary masses assigned to the
            observed stars by each synthetic cluster (``nan`` for single systems),
            and the mask of the stars with valid photometry
        :rtype: tuple[dict, np.ndarray, np.ndarray]
        """
        from .modules import stellar_stats_funcs as ssf

        # Stellar masses
        m1_med, m1_std, m2_med, m2_std, nan_msk, m2_obs = ssf.get_stellar_masses(
            self.cluster_mag,
            self.cluster_e_mag,
            self.cluster_colors,
            self.cluster_e_colors,
            self.m_ini_idx,
            synthcls,
        )

        # Binary probabilities
        binar_probs = ssf.get_binar_probabilities(m2_obs)

        # BSS probabilities
        mag_offset_1, col_offset_1, mag_offset_2, col_offset_2 = bss_box
        bss_probs = ssf.get_bss_probabilities(
            self.cluster_mag,
            self.cluster_colors,
            turn_off_points,
            len(synthcls),
            mag_offset_1,
            mag_offset_2,
            col_offset_1,
            col_offset_2,
        )

        binar_probs = (1 - bss_probs) * binar_probs
        single_probs = 1 - binar_probs - bss_probs

        N_cluster = len(self.cluster_mag)
        bss_mass_median, bss_mass_stddev = (
            np.full(N_cluster, np.nan),
            np.full(N_cluster, np.nan),
        )
        if isochs_young is not None:
            bss_mass_median, bss_mass_stddev = ssf.get_bss_masses(
                self.cluster_mag,
                self.cluster_colors[0],
                extrinsic_pars,
                self.ext_law,
                self.ext_coefs,
                self.rng,
                isochs_young.met_age_dict,
                isochs_young.theor_tracks,
            )
        else:
            warnings.warn(
                "No young isochrones were provided, BSSs masses can not be estimated"
            )

        # Store as dictionary
        data = {
            "mass_single": m1_med,
            "mass_single_std": m1_std,
            "mass_binar": m2_med,
            "mass_binar_std": m2_std,
            "mass_bss": bss_mass_median,
            "mass_bss_std": bss_mass_stddev,
            "single_prob": single_probs,
            "binar_prob": binar_probs,
            "bss_prob": bss_probs,
        }

        # Assign all nans to stars with a photometric nan in any dimension
        for value in data.values():
            value[nan_msk] = np.nan

        # Add ID column if it exists
        if self.star_id is not None:
            data = {"star_id": self.star_id, **data}

        if nan_msk.sum() > 0:
            self._vp("", 0)
            warnings.warn(
                f"\nN={nan_msk.sum()} stars found with no valid photometric data. "
                + "These will be assigned 'nan' values\nfor masses and "
                + "probabilities"
            )

        return data, m2_obs, ~nan_msk

    def _cluster_masses(
        self,
        models: list[dict],
        synthcls: list[np.ndarray],
        radec_c: tuple[float, float] | None,
        rho_amb: float | None,
        gal_pot: dict[str, float],
        disrup: dict[str, float],
    ) -> dict[str, np.ndarray]:
        """Masses of the cluster for each sampled model, used by
        :py:meth:`cluster_properties`

        :param models: Sampled models
        :type models: list[dict]
        :param synthcls: Synthetic clusters of the sampled models
        :type synthcls: list[np.ndarray]
        :param radec_c: Center coordinates in ``(RA, DEC)``, or ``None``
        :type radec_c: tuple[float, float] | None
        :param rho_amb: Ambient density, or ``None``
        :type rho_amb: float | None
        :param gal_pot: Parameters of the Galaxy's potential
        :type gal_pot: dict[str, float]
        :param disrup: Parameters of the cluster's disruption
        :type disrup: dict[str, float]

        :return: Dictionary with an array of values (one per model) for each mass
        :rtype: dict[str, np.ndarray]
        """
        from .modules import cluster_mass_funcs as cmf

        # The number of stars in the synthetic clusters is not constant so we estimate
        # its median
        N_stars_isoch = int(np.median([np.shape(_)[-1] for _ in synthcls]))
        # Compare the number of observed vs generated synthetic stars
        if N_stars_isoch < self.N_stars_obs:
            warnings.warn(
                f"\nThe number of synthetic stars ({N_stars_isoch}) is smaller than "
                + f"the number of\nobserved stars ({self.N_stars_obs}). Increase "
                + "the 'max_mass' argument in `asteca.Synthetic()`\nfor a more "
                + "accurate mass estimation"
            )

        # Resolve radec_c if not explicitly provided
        if (
            rho_amb is None
            and radec_c is None
            and self.cluster_ra is not None
            and self.cluster_dec is not None
        ):
            radec_c = (
                float(np.median(self.cluster_ra)),
                float(np.median(self.cluster_dec)),
            )
            self._vp(
                "\nUsing ambient density estimated from observed cluster (ra, dec)", 1
            )

        # Ambient density for each model
        rho_amb_arr = None
        if rho_amb is not None:
            rho_amb_arr = np.full(len(models), rho_amb)
        elif radec_c is not None:
            rho_amb_arr = cmf.ambient_density(
                models,
                radec_c,
                gal_pot["M_B"],
                gal_pot["r_B"],
                gal_pot["M_D"],
                gal_pot["a"],
                gal_pot["b"],
                gal_pot["r_s"],
                gal_pot["M_s"],
            )
        else:
            warnings.warn(
                "\nNo 'rho_amb', 'radec_c', or (ra, dec) values in the calibrated "
                + "Cluster object:\nthe initial, evolutionary and dynamical masses "
                + "can not be estimated"
            )

        C_env, gamma, epsilon = disrup["C_env"], disrup["gamma"], disrup["epsilon"]
        masses_all = []
        for i, model in enumerate(models):
            # Estimate the actual mass, ie: the sum of the observed and photometric
            # masses
            M_obs, M_phot, M_a = cmf.get_M_actual(
                self.rng,
                self.m_ini_idx,
                self.st_dist_mass,
                synthcls[i],
            )

            if rho_amb_arr is None:
                masses_all.append([M_obs, M_phot, np.nan, np.nan])
                continue

            # Extract met and loga
            z_met, loga = model["met"], model["loga"]
            # Convert back to Z if values are stored as FeH. This is required for the
            # 'mu_ev' estimation in the stellar_evol_mass_loss() function
            if self.isochs.z_to_FeH is not None:
                z_met = self.isochs.z_to_FeH * 10**z_met

            # Dissolution parameter
            t0 = cmf.dissolution_param(C_env, epsilon, gamma, rho_amb_arr[i])

            # Fraction of the initial mass that remains after the mass loss by stellar
            # evolution
            mu_ev = cmf.stellar_evol_mass_loss(z_met, loga)

            # Initial mass
            M_i = cmf.minit_LGB05(loga, M_a, gamma, t0, mu_ev)

            # Obtain evolutionary and dynamical masses
            M_evol = max(0, M_i * (1 - mu_ev))
            M_dyn = max(0, M_i - M_evol - M_a)

            masses_all.append([M_obs, M_phot, M_evol, M_dyn])

        masses_all = np.array(masses_all, dtype=float).T

        M_actual = masses_all[0] + masses_all[1]
        M_init = M_actual + masses_all[2] + masses_all[3]

        return {
            "M_init": M_init,
            "M_actual": M_actual,
            "M_obs": masses_all[0],
            "M_phot": masses_all[1],
            "M_evol": masses_all[2],
            "M_dyn": masses_all[3],
        }

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
        # # Generate displaced isochrone
        # fit_params_copy = dict(fit_params)

        from .modules import synth_model as sm

        try:
            # Extract from calibrated observed cluster
            max_mag_syn = self.max_mag_syn_obs
            err_dist_synth = self.err_dist_obs
            N_synth_stars = self.N_stars_obs
        except AttributeError:
            max_mag_syn = np.inf
            err_dist_synth = []
            N_synth_stars = int(N_stars)

        # Generate displaced isochrone and full synthetic CMD array
        isochrone, full_arr = sm.get_model(
            fit_params,
            self.met_age_dict,
            self.met_tol,
            self.age_tol,
            self.def_params,
            self.m_ini_idx,
            self.theor_tracks,
            self.ext_law,
            self.ext_coefs,
            self.rand_floats,
            self.st_dist_mass,
            self.mass_idx,
            max_mag_syn,
            err_dist_synth,
            N_synth_stars,
            isoch_flag=True,
        )

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
        isoch_mag, isoch_mass = isochrone[0], isochrone[self.m_ini_idx]
        mag_full, mass_full = full_arr[0], full_arr[self.m_ini_idx]
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

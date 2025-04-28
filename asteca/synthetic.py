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
    :param def_params: Default values for all the required fundamental parameters
    :type def_params: dict
    :param ext_law: Extinction law, one of ``CCMO, GAIADR3``.
        If ``GAIADR3`` is selected, the magnitude and first
        color defined in the :py:class:`Isochrones <asteca.isochrones.Isochrones>` and
        :py:class:`Cluster <asteca.cluster.Cluster>` objects are assumed to be Gaia's
        (E)DR3 **G** and **(BP-RP)** respectively. The second color (if defined) will
        always be affected by the ``CCMO`` model, defaults to ``CCMO``
    :type ext_law: str
    :param DR_distribution: Distribution function for the differential reddening,
        one of ``uniform, normal``; defaults to ``uniform``
    :type DR_distribution: str
    :param IMF_name: Name of the initial mass function used to populate the isochrones,
        one of ``salpeter_1955, kroupa_2001, chabrier_2014``;
        defaults to ``chabrier_2014``
    :type IMF_name: str
    :param max_mass: Maximum total initial mass. Should be large enough to allow
        generating as many synthetic stars as observed stars, defaults to ``100_000``
    :type max_mass: int
    :param gamma: Distribution function for the mass ratio of the binary systems,
        float or one of ``D&K, fisher_stepped, fisher_peaked, raghavan``;
        defaults to ``D&K``
    :type gamma: float | str
    :param seed: Random seed. If ``None`` a random integer will be generated and used,
        defaults to ``None``
    :type seed: int | None
    :param verbose: Verbose level. A value of ``0`` hides all output, defaults to ``1``
    :type verbose: int

    :raises ValueError: If any of the attributes is not recognized as a valid option
    """

    def __init__(
        self,
        isochs: Isochrones,
        def_params: dict = {
            "met": 0.0152,
            "loga": 8.0,
            "alpha": 0.09,
            "beta": 0.94,
            "Rv": 3.1,
            "DR": 0.0,
            "Av": 0.2,
            "dm": 9.0,
        },
        ext_law: str = "CCMO",
        DR_distribution: str = "uniform",
        IMF_name: str = "chabrier_2014",
        max_mass: int = 10_000,
        gamma: float | str = "D&K",
        seed: int | None = None,
        verbose: int = 1,
    ) -> None:
        self.isochs = isochs
        self.def_params = def_params
        self.ext_law = ext_law
        self.DR_distribution = DR_distribution
        self.IMF_name = IMF_name
        self.max_mass = max_mass
        self.gamma = gamma
        self.seed = seed
        self.verbose = verbose

        # Set seed
        self.rng = np.random.default_rng(self.seed)

        # Check gamma distribution
        gammas = ("D&K", "fisher_stepped", "fisher_peaked", "raghavan")
        if isinstance(self.gamma, str):
            if self.gamma not in gammas:
                raise ValueError(
                    f"gamma '{self.gamma}' not recognized. Should be one of {gammas}"
                )

        # Check differential reddening function
        DR_funcs = ("uniform", "normal")
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
        self.st_dist_mass, self.st_dist_mass_ordered = scp.sample_imf(
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
        if self.ext_law == "GAIADR3":
            if (
                self.isochs.magnitude_effl is not None
                or self.isochs.color_effl is not None
            ):
                warnings.warn(
                    f"\nExtinction law '{self.ext_law}' does not require effective lambda"
                    + " values for the\nmagnitude and first color (assumed to be 'G' "
                    + "and 'BP-RP', respectively)."
                )

        # Generate random floats used by `generate()`
        self.rand_floats = scp.randVals(self.rng, self.theor_tracks, self.st_dist_mass)

        # Store for internal usage
        self.met_age_dict = self.isochs.met_age_dict

        # Check that the ranges are respected
        for par in ("met", "loga"):
            pmin, pmax = min(self.met_age_dict[par]), max(self.met_age_dict[par])
            if self.def_params[par] < pmin or self.def_params[par] > pmax:
                warnings.warn(
                    f"Parameter {par}={self.def_params[par]} is out of range: [{pmin} - {pmax}]"
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

        See the :ref:`synthetic_module` section for more details.

        :param cluster: :py:class:`Cluster` object with the
            processed data from your observed cluster
        :type cluster: Cluster

        :raises ValueError:
            - If required are 'Cluster' attributes missing.
            - If the number of colors defined in the :py:class:`Cluster` and
              :py:class:`Synthetic` objects do not match
        """
        from .modules import synth_cluster_priv as scp

        for attr in ("mag", "color", "ra", "dec"):
            if getattr(cluster, attr) is None:
                raise ValueError(
                    f"Attribute '{attr}' is required to calibrate 'Cluster' object."
                )
        if cluster.e_mag is None:
            raise ValueError(
                "Attribute 'e_mag' is required to calibrate 'Cluster' object."
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
        if len(cluster.mag) > cluster.N_clust_max:
            raise ValueError(
                f"The number of stars in this `Cluster` object ({len(cluster.mag)})\n"
                + f"is larger than the 'N_clust_max={cluster.N_clust_max}' parameter.\n"
                + "Either define a `Cluster` object with fewer members or increase "
                + "the 'N_clust_max' value."
            )

        # Calibrate parameters using the observed cluster
        self.max_mag_syn_obs = max(cluster.mag)
        self.N_stars_obs = len(cluster.mag)
        self.err_dist_obs = scp.error_distribution(
            cluster.mag,
            cluster.e_mag,
            cluster.e_color,
            cluster.e_color2,
            self.rand_floats["norm"][1],
        )

        self._vp("\nCalibrated observed cluster")
        self._vp(f"N_stars_obs    : {self.N_stars_obs}", 1)
        self._vp(f"Max magnitude  : {self.max_mag_syn_obs:.2f}", 1)
        self._vp("Error distribution loaded", 1)

        # Used by the `stellar_masses()` method
        self.cluster_mag = cluster.mag
        self.cluster_colors = [cluster.color, cluster.color2]

        # Used by the `cluster_masses()` method
        self.cluster_ra = cluster.ra
        self.cluster_dec = cluster.dec

    def generate(self, params: dict, N_stars: int = 100) -> np.ndarray:
        """Generate a synthetic cluster.

        The synthetic cluster is generated according to the parameters given in
        the ``params`` dictionary and the already calibrated
        :py:class:`Synthetic` object.

        :param params: Dictionary containing the values for the fundamental parameters.
            The dictionary must include values for all the parameters, e.g.:
            ``params = {met: 0.0152, loga: 8.1, alpha: 0.1, beta: 1, Av: 0.2, DR: 0., Rv: 3.1, dm: 9.7}``
        :type params: dict
        :param N_stars: Number of synthetic stars to generate; defaults to ``100``
        :type N_stars: int

        :return: Returns a ``np.array`` containing a synthetic cluster
            with the data ``[mag, c1, (c2), mass, mass_b]``, where ``mag`` is
            the magnitude, ``c1`` is the color, ``c2`` is the optional second color,
            and ``mass, mass_b`` are the masses of the single and secondary components
            of the binary systems, respectively (if generated). If the system is a
            single star, then ``mass_b==np.nan``.
        :rtype: np.ndarray
        """
        from .modules import synth_cluster_priv as scp

        try:
            # Extract from calibrated observed cluster
            max_mag_syn = self.max_mag_syn_obs
            N_synth_stars = self.N_stars_obs
            err_dist_synth = self.err_dist_obs
        except AttributeError:
            max_mag_syn = np.inf
            N_synth_stars = int(N_stars)
            err_dist_synth = []

        # Return proper values for fixed parameters and parameters required
        # for the (z, log(age)) isochrone averaging.
        met, loga, alpha, beta, av, dr, rv, dm, ml, mh, al, ah = scp.properModel(
            self.met_age_dict, self.def_params, params
        )

        # If (z, a) are both fixed, use the single processed isochrone
        if ml == al == mh == ah == 0:
            # The np.array() is important to avoid overwriting 'theor_tracks'
            isochrone = np.array(self.theor_tracks[0][0])
        else:
            # Generate a weighted average isochrone from the (z, log(age)) values in
            # the 'model'.
            isochrone = scp.zaWAverage(
                self.theor_tracks,
                self.met_age_dict,
                self.m_ini_idx,
                met,
                loga,
                ml,
                mh,
                al,
                ah,
            )

        binar_flag = True
        if alpha == 0.0 and beta == 0.0:
            binar_flag = False

            # # TODO: this was not tested thoroughly (April 2025)
            # # Remove binary photometry
            # isochrone = isochrone[: self.m_ini_idx + 2]
            # # TODO: this was not tested thoroughly

        # Move theoretical isochrone using the distance modulus
        isoch_moved = scp.move_isochrone(isochrone, binar_flag, self.m_ini_idx, dm)

        # Apply extinction correction
        isoch_extin = scp.extinction(
            self.ext_law,
            self.ext_coefs,
            self.rand_floats["norm"][0],
            self.rand_floats["unif"][0],
            self.DR_distribution,
            self.m_ini_idx,
            binar_flag,
            av,
            dr,
            rv,
            isoch_moved,
        )

        # Remove isochrone stars beyond the maximum magnitude
        isoch_cut = scp.cut_max_mag(isoch_extin, max_mag_syn)
        if not isoch_cut.any():
            return np.array([])

        # This is an internal trick to return the array at this point. It is used
        # to generate an isochrone to plot
        if N_stars == -1:
            return isoch_cut

        # Interpolate IMF's sampled masses into the isochrone.
        isoch_mass = scp.mass_interp(
            isoch_cut,
            self.m_ini_idx,
            self.st_dist_mass[ml][al],
            N_synth_stars,
            binar_flag,
        )
        if not isoch_mass.any():
            return np.array([])

        # import matplotlib.pyplot as plt
        # # plt.title('2000, steps 1')
        # plt.title('2000, steps 2')
        # plt.scatter(isoch_mass[1], isoch_mass[0], alpha=.25)
        # plt.scatter(isoch_mass[5], isoch_mass[4], alpha=.25)
        # plt.gca().invert_yaxis()
        # plt.show()

        # Assignment of binarity.
        isoch_binar = scp.binarity(
            alpha,
            beta,
            binar_flag,
            self.m_ini_idx,
            self.rand_floats["unif"][1],
            isoch_mass,
        )

        # Assign errors according to errors distribution.
        synth_clust = scp.add_errors(isoch_binar, err_dist_synth)

        return synth_clust

    def get_models(
        self,
        model: dict[str, float],
        model_std: dict[str, float],
        N_models: int = 200,
    ) -> None:
        """Generate random sampled models from the selected solution. Use these models
        to generate full synthetic clusters.

        :param model: Dictionary with the values for the fundamental parameters
        :type model: dict[str, float]
        :param model_std: Dictionary with the standard deviations for the fundamental
            parameters in the ``model`` argument
        :type model_std: dict[str, float]
        :param N_models: Number of sampled models, defaults to ``200``
        :type N_models: int

        :raises ValueError: If any of the (met, age) parameters are out of range or
            if all the synthetic arrays generated are empty
        """
        from .modules import mass_binary as mb

        # Update dictionaries with default values, if required
        model = self.def_params | model
        for k in self.def_params.keys():
            if k not in model_std.keys():
                model_std[k] = 0

        pars_k = ("met", "loga", "alpha", "beta", "Rv", "DR", "Av", "dm")
        self._vp("\nGenerate synthetic models", 1)
        self._vp("par            : mean +/ STDDEV", 2)
        for k in pars_k:
            self._vp(
                f"{k:<15}: {round(model[k], 4)} +/- {round(model_std[k], 4)}",
                2,
            )

        # Check isochrones ranges
        for par in ("met", "loga"):
            try:
                pmin, pmax = min(self.met_age_dict[par]), max(self.met_age_dict[par])
                if model[par] < pmin or model[par] > pmax:
                    raise ValueError(
                        f"Parameter '{par}' out of range: [{pmin} - {pmax}]"
                    )
            except KeyError:
                pass

        # Generate random model parameter values
        sampled_models = mb.ranModels(model, model_std, N_models, self.rng)

        # Generate the synthetic arrays using the above models
        sampled_models_no_empty, sampled_synthcls = [], []
        for i, smodel in enumerate(sampled_models):
            isoch = self.generate(smodel)
            # Do not store models that result in empty arrays
            if not isoch.any():
                continue
            sampled_models_no_empty.append(smodel)
            sampled_synthcls.append(isoch)

        N_models_non_empty = len(sampled_models_no_empty)
        if N_models_non_empty == 0:
            raise ValueError("All sampled models produced empty synthetic clusters.")
        if len(sampled_models_no_empty) < N_models:
            warnings.warn(
                f"Empty synthetic arrays generated: {N_models - N_models_non_empty}"
            )

        # cluster_masses
        self.sampled_models = sampled_models_no_empty
        # stellar_masses, cluster_masses
        self.sampled_synthcls = sampled_synthcls

        self._vp(f"N_models       : {len(sampled_models_no_empty)}", 1)
        self._vp("Attributes stored in Synthetic object", 1)

    def stellar_masses(self) -> dict:
        """Estimate individual masses for the observed stars, along with their binary
        probabilities (if binarity was estimated).

        :return: Data frame containing per-star primary and secondary masses along with
            their uncertainties, and their probability of being a binary system
        :rtype: dict

        :raises ValueError: If the `get_models()` method was not previously called
        """
        if (
            hasattr(self, "cluster_mag") is False
            or hasattr(self, "cluster_colors") is False
        ):
            raise ValueError(
                "This method requires running the 'calibrate()' method first\n"
                + "with an observed 'Cluster' object"
            )

        if hasattr(self, "sampled_synthcls") is False:
            raise ValueError(
                "This method requires running the 'get_models()' method first"
            )

        from .modules import mass_binary as mb

        # Observed photometry
        cluster_colors = [self.cluster_colors[0]]
        if self.cluster_colors[1] is not None:
            cluster_colors.append(self.cluster_colors[1])
        obs_phot = np.array([self.cluster_mag] + [_ for _ in cluster_colors])
        # Replace nans in mag and colors to avoid crashing KDTree()
        nan_msk = np.full(obs_phot.shape[1], False)
        for ophot in obs_phot:
            nan_msk = nan_msk | np.isnan(ophot)
        obs_phot[:, nan_msk] = -10.0
        obs_phot = obs_phot.T

        #
        close_stars_idxs = []
        for isoch in self.sampled_synthcls:
            close_stars_idxs.append(mb.get_close_idxs(self.m_ini_idx, obs_phot, isoch))

        #
        m12_masses = []
        for i, isoch in enumerate(self.sampled_synthcls):
            m1_obs, m2_obs = mb.get_m1m2(self.m_ini_idx, isoch, close_stars_idxs[i])
            m12_masses.append([m1_obs, m2_obs])
        m12_masses = np.array(m12_masses)

        # Primary masses (median + stddev)
        m1_med = np.median(m12_masses[:, 0, :], 0)
        m1_std = np.std(m12_masses[:, 0, :], 0)
        # Secondary masses  (median + stddev). Hide 'All-nan slice' warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m2_med = np.nanmedian(m12_masses[:, 1, :], 0)
            # m2 can not be larger than m1
            m2_med = np.min([m1_med, m2_med], 0)
            m2_std = np.nanstd(m12_masses[:, 1, :], 0)

        # Binary probability per observed star. Check how many times the secondary mass
        # of an observed star was assigned a value 'not np.nan', i.e.: was identified
        # as a binary system. Dividing this value by the number of synthetic models
        # used, results in the per observed star probability of being a binary system.
        binar_prob = (~np.isnan(m12_masses[:, 1, :])).sum(0) / m12_masses.shape[0]

        # Store as dictionary
        data = {
            "m1": m1_med,
            "m1_std": m1_std,
            "m2": m2_med,
            "m2_std": m2_std,
            "binar_prob": binar_prob,
        }

        # Assign all nans to stars with a photometric nan in any dimension
        for value in data.values():
            value[nan_msk] = np.nan

        if nan_msk.sum() > 0:
            warnings.warn(
                f"\nN={nan_msk.sum()} stars found with no valid photometric data. "
                + "These will be assigned 'nan' values\nfor masses and "
                + "binarity probability"
            )

        self._vp("\nStellar masses and binary probabilities estimated", 1)

        return data

    def binary_fraction(
        self,
        binar_prob: np.ndarray,
        Nsamples: int = 10_000,
    ) -> tuple[float, float]:
        """Estimate the total binary fraction for the observed cluster, given
        the per-star probability of being a binary system.

        :param binar_prob: Per-star binary system probabilities
        :type binar_prob: np.ndarray
        :param Nsamples: Number of samples generated to produce the final values;
            defaults to ``10_000``.
        :type Nsamples: int

        :return: Median and STDDEV values for the total binary fraction
        :rtype: tuple[float, float]

        :raises ValueError: If the `get_models()` method was not previously called
        """
        if binar_prob.min() < 0.0 or binar_prob.max() > 1.0:
            raise ValueError("The values in 'binar_prob' must be in the range [0, 1]")

        # Calculates the fraction of stars that are binaries for Nsamples random samples
        N_stars = len(binar_prob)
        # Observe systems (Nsamples) and store how many are single/binaries
        b_fr_vals = self.rng.random((Nsamples, N_stars)) < binar_prob
        # Average the binary fraction for each observation
        b_fr_sum = np.sum(b_fr_vals, axis=1) / N_stars
        # Return the median and STDDEV values
        bfr_med, bfr_std = np.median(b_fr_sum), np.std(b_fr_sum)
        self._vp(f"Binary fraction: {bfr_med:.3f} +/- {bfr_std:.3f}", 1)

        return float(bfr_med), float(bfr_std)

    def cluster_masses(
        self,
        rho_amb: float | None = None,
        M_B: float = 2.5e10,
        r_B: float = 0.5e3,
        M_D: float = 7.5e10,
        a: float = 5.4e3,
        b: float = 0.3e3,
        M_s: float = 1.87e11,
        r_s: float = 15.19e3,
        C_env: float = 810e6,
        gamma: float = 0.62,
        epsilon: float = 0.08,
    ) -> dict:
        """Estimate the different total masses for the observed cluster.

        The returned dictionary contains distributions for
        ``M_init, M_actual, M_obs, M_phot, M_evol, M_dyn``, where:

        - ``M_init`` : Initial mass
        - ``M_actual`` : Actual mass
        - ``M_obs`` : Observed mass
        - ``M_phot``: Photometric mass, ie: mass not observed that is located below the
          maximum observed magnitude
        - ``M_evol``: Mass lost via stellar evolution
        - ``M_dyn`` : Mass lost via dynamical effects

        The actual and initial masses can also be obtained as:

        - ``M_actual = M_obs + M_phot``
        - ``M_init = M_actual + M_evol + M_dyn``

        :param rho_amb: Ambient density. If ``None``, it is estimated using the
            cluster's position and a model for the Galaxy's potential; defaults to
            ``None``.
        :type rho_amb: float | None
        :param M_B: Bulge mass (in solar masses); defaults to ``2.5e10`` (from
            `Haghi et al. 2015 <https://doi.org/10.1093/mnras/stv827>`__, Table 1)
        :type M_B: float
        :param r_B: Characteristic radius of the bulge (in pc); defaults to ``0.5e3``
            (from `Haghi et al. 2015 <https://doi.org/10.1093/mnras/stv827>`__,
            Table 1)
        :type r_B: float
        :param M_D: Disc mass (in solar masses); defaults to ``7.5e10`` (from
            `Haghi et al. 2015 <https://doi.org/10.1093/mnras/stv827>`__, Table 1)
        :type M_D: float
        :param a: Disc scale radius (in pc); defaults to ``5.4e3`` (from
            `Haghi et al. 2015 <https://doi.org/10.1093/mnras/stv827>`__, Table 1)
        :type a: float
        :param b: Disc scaleheight (in pc); defaults to ``0.3e3`` (from
            `Haghi et al. 2015 <https://doi.org/10.1093/mnras/stv827>`__, Table 1)
        :type b: float
        :param M_s: Dark matter halo mass (in solar masses); defaults to ``1.87e11``
            (from `Sanderson et al. 2017
            <https://iopscience.iop.org/article/10.3847/1538-4357/aa5eb4>`__, Table 1)
        :type M_s: float
        :param r_s: Dark matter halo scale radius (in pc); defaults to ``15.19e3`` (from
            `Sanderson et al. 2017
            <https://iopscience.iop.org/article/10.3847/1538-4357/aa5eb4>`__, Table 1)
        :type r_s: float
        :param C_env: Constant related to the disruption time (in Myr); defaults to
            ``810e6`` (from `Lamers, Gieles & Zwart 2005
            <https://www.aanda.org/articles/aa/abs/2005/01/aa1476/aa1476.html>`__)
        :type C_env: float
        :param gamma: Constant related to the disruption time (no units); defaults to
            ``0.62`` (from `Lamers, Gieles & Zwart 2005
            <https://www.aanda.org/articles/aa/abs/2005/01/aa1476/aa1476.html>`__)
        :type gamma: float
        :param epsilon: Eccentricity of the orbit; defaults to ``0.08`` (from
            `Angelo et al. (2023) <https://doi.org/10.1093/mnras/stad1038>`__)
        :type epsilon: float

        :raises ValueError:
            If no synthetic models were generated via the get_models() method.

        :return: Dictionary with the mass distributions for the initial, actual,
            observed, photometric, evolutionary, and dynamical masses:
            ``M_init, M_actual, M_obs, M_phot, M_evol, M_dyn``
        :rtype: dict
        """
        if hasattr(self, "sampled_synthcls") is False:
            raise ValueError(
                "This method requires running the get_models() method first"
            )

        from .modules import mass_binary as mb

        # The number of stars in the synthetic clusters is not constant so we estimate
        # its median
        N_stars_isoch = int(np.median([np.shape(_)[-1] for _ in self.sampled_synthcls]))
        # Compare the number of observed vs generated synthetic stars
        if N_stars_isoch < self.N_stars_obs:
            warnings.warn(
                f"\nNumber of synthetic stars ({N_stars_isoch}) is smaller than "
                + f"observed stars ({self.N_stars_obs}). Increase the 'max_mass' "
                + "argument for a more accurate mass estimation"
            )

        # Estimate the center coordinates from the cluster's median values
        radec_c = (np.median(self.cluster_ra), np.median(self.cluster_dec))  # pyright: ignore

        if rho_amb is not None:
            # Input float to array
            rho_amb_arr = np.ones(len(self.sampled_models)) * rho_amb
        else:
            # Obtain galactic vertical distance and distance to center
            Z, R_GC, R_xy = mb.galactic_coords(self.sampled_models, radec_c)
            rho_amb_arr = mb.ambient_density(
                M_B, r_B, M_D, a, b, r_s, M_s, Z, R_GC, R_xy
            )

        masses_all = []
        for i, model in enumerate(self.sampled_models):
            # Extract met and loga
            z_met, loga = model["met"], model["loga"]

            # Convert back to Z if values are stored as FeH. This is required for the
            # 'mu_ev' estimation in the stellar_evol_mass_loss() function
            if self.isochs.z_to_FeH is not None:
                z_met = self.isochs.z_to_FeH * 10**z_met

            # Estimate the actual mass, ie: the sum of the observed and photometric
            # masses
            sampled_synthcl = self.sampled_synthcls[i]
            M_obs, M_phot, M_a = mb.get_M_actual(
                self.rng,
                self.m_ini_idx,
                self.st_dist_mass,
                self.st_dist_mass_ordered,
                sampled_synthcl,
            )

            # Dissolution parameter
            t0 = mb.dissolution_param(C_env, epsilon, gamma, rho_amb_arr[i])

            # Fraction of the initial mass that is lost by stellar evolution
            mu_ev = mb.stellar_evol_mass_loss(z_met, loga)

            # Initial mass
            M_i = mb.minit_LGB05(loga, M_a, gamma, t0, mu_ev)

            # Obtain evolutionary and dynamical masses
            M_evol = max(0, M_i * (1 - mu_ev))
            M_dyn = max(0, M_i - M_evol - M_a)

            masses_all.append([M_obs, M_phot, M_evol, M_dyn])

        masses_all = np.array(masses_all).T

        M_actual = masses_all[0] + masses_all[1]
        M_init = M_actual + masses_all[2] + masses_all[3]

        self._vp("\nMass values estimated", 1)

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
    ) -> np.ndarray:
        """Generate an isochrone for plotting.

        The isochrone is generated using the fundamental parameter values
        given in the ``fit_params`` dictionary.

        :param fit_params: Dictionary with the values for the fundamental parameters
            that were **not** included in the ``fix_params`` dictionary when the
            :py:class:`Synthetic` object was calibrated (:py:meth:`calibrate` method).
        :type fit_params: dict
        :param color_idx: Index of the color to plot. If ``0`` (default), plot the
            first color. If ``1`` plot the second color. Defaults to ``0``
        :type color_idx: int

        :raises ValueError: If either parameter (met, age) is outside of allowed range

        :return: Array with the isochrone data to plot
        :rtype: np.ndarray
        """
        # Generate displaced isochrone
        fit_params_copy = dict(fit_params)

        # Check isochrones ranges
        for par in ("met", "loga"):
            try:
                pmin, pmax = min(self.met_age_dict[par]), max(self.met_age_dict[par])
                if fit_params_copy[par] < pmin or fit_params_copy[par] > pmax:
                    raise ValueError(
                        f"Parameter '{par}' out of range: [{pmin} - {pmax}]"
                    )
            except KeyError:
                pass

        # Generate physical synthetic cluster to extract the max mass
        max_mass = self.generate(fit_params_copy)[self.m_ini_idx].max()

        # Generate displaced isochrone, The 'N_stars=-1' indicates to return the
        # array right after the magnitude cut
        fit_params_copy["DR"] = 0.0
        isochrone = self.generate(fit_params_copy, N_stars=-1)

        # Apply max mass filter to isochrone
        msk = isochrone[self.m_ini_idx] < max_mass

        # Generate proper array for plotting
        isochrone = np.array([isochrone[0], isochrone[color_idx + 1]])[:, msk]

        return isochrone

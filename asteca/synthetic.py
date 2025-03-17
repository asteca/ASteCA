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

    See the :ref:`synth_clusters` section for more details.

    :param isochs: :py:class:`Isochrones <asteca.isochrones.Isochrones>` object with
        the loaded files for the theoretical isochrones
    :type isochs: Isochrones
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
        ext_law: str = "CCMO",
        DR_distribution: str = "uniform",
        IMF_name: str = "chabrier_2014",
        max_mass: int = 20_000,
        gamma: float | str = "D&K",
        seed: int | None = None,
        verbose: int = 1,
    ) -> None:
        self.isochs = isochs
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

        self._vp("\nInstantiating synthetic...")

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

        # Generate random floats used by `synth_clusters.synthcl_generate()`
        self.rand_floats = scp.randVals(self.rng, self.theor_tracks, self.st_dist_mass)

        # Store for internal usage
        self.met_age_dict = self.isochs.met_age_dict

        self._vp(f"IMF            : {self.IMF_name}", 1)
        self._vp(f"Max init mass  : {self.max_mass}", 1)
        self._vp(f"Gamma dist     : {self.gamma}", 1)
        self._vp(f"Extinction law : {self.ext_law}", 1)
        self._vp(f"Diff reddening : {self.DR_distribution}", 1)
        self._vp(f"Random seed    : {self.seed}", 1)
        self._vp("Synthetic clusters object generated")

    def _vp(self, mssg: str, level: int = 0) -> None:
        """Verbose print method"""
        if self.verbose > level:
            print(mssg)

    def calibrate(self, cluster: Cluster, fix_params: dict = {}):
        """Calibrate a :py:class:`Synthetic` object based on a
        :py:class:`Cluster <asteca.cluster.Cluster>` object and a dictionary of fixed
        fundamental parameters (``fix_params``).

        Use the data obtained from your observed cluster stored in the
        :py:class:`Cluster <asteca.cluster.Cluster>` object, to calibrate a
        :py:class:`Synthetic` object. Additionally, a dictionary of fixed fundamental
        parameters (metallicity, age, distance, extinction, etc.) can be passed.

        See the :ref:`synth_clusters` section for more details.

        :param cluster: :py:class:`Cluster <asteca.cluster.Cluster>` object with the
            processed data from your observed cluster
        :type cluster: Cluster
        :param fix_params: Dictionary with the values for the fixed parameters (if any),
            defaults to ``{}``
        :type fix_params: dict

        :raises ValueError:
            -If the number of colors defined in the
            :py:class:`Cluster <asteca.cluster.Cluster>` and
            :py:class:`Synthetic <asteca.synthetic.Synthetic>` objects do not match
            -If the metallicity or age parameters are not fixed to a
            single value but there ranges are.
        """
        from .modules import synth_cluster_priv as scp

        # Check that the number of colors match
        if self.isochs.color2_effl is not None and cluster.color2 is None:
            raise ValueError(
                "Two colors were defined in 'synthetic' but a single color\n"
                + "was defined in 'cluster'."
            )
        if self.isochs.color2_effl is None and cluster.color2 is not None:
            raise ValueError(
                "Two colors were defined in 'cluster' but a single color\n"
                + "was defined in 'synthetic'."
            )
        if len(cluster.mag_v) > cluster.N_clust_max:
            raise ValueError(
                f"The number of stars in this `Cluster` object ({len(cluster.mag_v)})\n"
                + f"is larger than the N_clust_max={cluster.N_clust_max} parameter.\n"
                + "Either define a `Cluster` object with fewer members or increase "
                + "the N_clust_max value."
            )

        # Data used by the `generate()` method
        self.m_ini_idx = 2  # (0->mag, 1->color, 2->mass_ini)
        if self.isochs.color2_effl is not None:
            self.m_ini_idx = 3  # (0->mag, 1->color, 2->color2, 3->mass_ini)

        self.max_mag_syn = max(cluster.mag_v)
        self.N_obs_stars = len(cluster.mag_v)
        self.err_dist = scp.error_distribution(
            cluster.mag_v,
            cluster.e_mag_v,
            cluster.e_colors_v,
            self.rand_floats["norm"][1],
        )

        # Used by the `get_models()` method and its result by the `stellar_masses()`
        # and `binary_fraction()` methods
        self.mag_v = cluster.mag_v
        self.colors_v = cluster.colors_v

        self.fix_params = fix_params

        self.binar_flag = True
        if "alpha" in list(fix_params.keys()) and "beta" in list(fix_params.keys()):
            if fix_params["alpha"] == 0.0 and fix_params["beta"] == 0.0:
                self.binar_flag = False

        # Check that the ranges are respected
        for par in ("met", "loga"):
            if par not in self.fix_params.keys():
                N_par = len(self.met_age_dict[par])
                if N_par == 1:
                    raise ValueError(
                        f"Parameter '{par}' is not fixed and its range is limited to "
                        + f"a single value: {self.met_age_dict[par]}"
                    )
            else:
                pmin, pmax = min(self.met_age_dict[par]), max(self.met_age_dict[par])
                if self.fix_params[par] < pmin or self.fix_params[par] > pmax:
                    raise ValueError(
                        f"Parameter {par}={self.fix_params[par]} out of range: [{pmin} - {pmax}]"
                    )

    def generate(
        self, fit_params: dict, plot_flag: bool = False, full_arr_flag: bool = False
    ) -> np.ndarray:
        """Generate a synthetic cluster.

        The synthetic cluster is generated according to the parameters given in
        the ``fit_params`` dictionary and the already calibrated
        :py:class:`Synthetic` object.

        :param fit_params: Dictionary with the values for the fundamental parameters
            that were **not** included in the ``fix_params`` dictionary when the
            :py:class:`Synthetic` object was calibrated
            (:py:meth:`calibrate` method).
        :type fit_params: dict
        :param plot_flag: If ``True`` returns the isochrone after the maximum magnitude
            cut  is applied. Used mainly for plotting, defaults to ``False``
        :type plot_flag: bool
        :param full_arr_flag: If ``True`` returns the full array for the synthetic
            cluster, including the binary data (if any). Used mainly for plotting,
            defaults to ``False``
        :type full_arr_flag: bool

        :return: By default it returns a ``np.array`` containing a synthetic cluster
            with the shape ``[mag, c1, (c2)]``, where ``mag`` is the magnitude
            dimension, and``c1`` and ``c2`` (last one is optional) are the color
            dimension(s). This changes depending on the flags above.
        :rtype: np.ndarray
        """
        from .modules import synth_cluster_priv as scp

        # Return proper values for fixed parameters and parameters required
        # for the (z, log(age)) isochrone averaging.
        met, loga, alpha, beta, av, dr, rv, dm, ml, mh, al, ah = scp.properModel(
            self.met_age_dict, self.fix_params, fit_params
        )

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

        # Move theoretical isochrone using the distance modulus
        isoch_moved = scp.move_isochrone(isochrone, self.m_ini_idx, dm)

        # Apply extinction correction
        isoch_extin = scp.extinction(
            self.ext_law,
            self.ext_coefs,
            self.rand_floats["norm"][0],
            self.rand_floats["unif"][0],
            self.DR_distribution,
            self.m_ini_idx,
            self.binar_flag,
            av,
            dr,
            rv,
            isoch_moved,
        )

        # Remove isochrone stars beyond the maximum magnitude
        isoch_cut = scp.cut_max_mag(isoch_extin, self.max_mag_syn)
        if not isoch_cut.any():
            return np.array([])
        if plot_flag:
            return isoch_cut

        # Interpolate IMF's sampled masses into the isochrone.
        isoch_mass = scp.mass_interp(
            isoch_cut, self.m_ini_idx, self.st_dist_mass[ml][al], self.N_obs_stars
        )
        if not isoch_mass.any():
            return np.array([])

        # Assignment of binarity.
        isoch_binar = scp.binarity(
            alpha,
            beta,
            self.binar_flag,
            self.m_ini_idx,
            self.rand_floats["unif"][1],
            isoch_mass,
        )

        # Assign errors according to errors distribution.
        synth_clust = scp.add_errors(isoch_binar, self.err_dist)

        if full_arr_flag:
            return synth_clust
        return synth_clust[: self.m_ini_idx]

    def get_models(
        self,
        model: dict[str, float],
        model_std: dict[str, float],
        N_models: int = 200,
    ) -> None:
        """Generate random sampled models from the selected solution. Use these models
        to generate full synthetic clusters.

        :param model: Dictionary with the values for the fundamental parameters that
            were **not** included in the ``fix_params`` dictionary when the
            :py:class:`Synthetic` object was calibrated
            (:py:meth:`synthetic.Synthetic.calibrate` method)
        :type model: dict[str, float]
        :param model_std: Dictionary with the standard deviations for the fundamental
            parameters in the ``model`` argument
        :type model_std: dict[str, float]
        :param N_models: Number of sampled models, defaults to ``200``
        :type N_models: int

        :raises ValueError: If any of the (met, age) parameters are out of range
        """
        from .modules import mass_binary as mb

        self._vp("\nGenerate synthetic models...", 1)
        self._vp(f"N_models       : {N_models}", 1)
        self._vp(
            "Model          :"
            + ", ".join(f"{k}: {round(v, 3)}" for k, v in model.items()),
            2,
        )
        self._vp(
            "Model STDDEV   :"
            + ", ".join(f"{k}: {round(v, 3)}" for k, v in model_std.items()),
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

        # Observed photometry
        obs_phot = np.array([self.mag_v] + [_ for _ in self.colors_v])
        # Replace nans in mag and colors to avoid crashing KDTree()
        nan_msk = np.full(obs_phot.shape[1], False)
        for ophot in obs_phot:
            nan_msk = nan_msk | np.isnan(ophot)
        obs_phot[:, nan_msk] = -10.0
        obs_phot = obs_phot.T

        sampled_models = mb.ranModels(model, model_std, N_models, self.rng)

        sampled_synthcls, close_stars_idxs = [], []
        remove_model_index = []
        for i, smodel in enumerate(sampled_models):
            isoch = self.generate(smodel, full_arr_flag=True)
            if not isoch.any():
                remove_model_index.append(i)
                continue
            sampled_synthcls.append(isoch)
            idxs = mb.get_close_idxs(self.m_ini_idx, obs_phot, isoch)
            close_stars_idxs.append(idxs)

        # Remove models associated to empty isochrones
        if len(remove_model_index) > 0:
            sampled_models = list(
                np.delete(np.array(sampled_models), np.array(remove_model_index))
            )
        self.sampled_models = sampled_models
        self.sampled_synthcls = sampled_synthcls
        self.close_stars_idxs = close_stars_idxs
        self.obs_nan_msk = nan_msk

        self._vp("Attributes stored in Synthetic object", 1)

    def stellar_masses(self) -> dict:
        """Estimate individual masses for the observed stars, along with their binary
        probabilities (if binarity was estimated).

        :return: Data frame containing per-star primary and secondary masses along with
            their uncertainties, and their probability of being a binary system
        :rtype: dict

        """
        from .modules import mass_binary as mb

        m12_masses = []
        for i, isoch in enumerate(self.sampled_synthcls):
            m1_obs, m2_obs = mb.get_m1m2(
                self.m_ini_idx, isoch, self.close_stars_idxs[i]
            )
            m12_masses.append([m1_obs, m2_obs])
        m12_masses = np.array(m12_masses)

        # Primary masses (median + stddev)
        m1_med = np.median(m12_masses[:, 0, :], 0)
        m1_std = np.std(m12_masses[:, 0, :], 0)
        # Secondary masses  (median + stddev). Hide 'All-nan slice' warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # TODO: this was changed in Python>3.10
            m2_med = np.nanmedian(m12_masses[:, 1, :], 0)
            # m2 can not be larger than m1
            m2_med = np.min([m1_med, m2_med], 0)
            m2_std = np.nanstd(m12_masses[:, 1, :], 0)

        # Binary probability per star
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
        nan_msk = self.obs_nan_msk
        for value in data.values():
            value[nan_msk] = np.nan

        if nan_msk.sum() > 0:
            warnings.warn(
                f"\nN={nan_msk.sum()} stars found with no valid photometric data. "
                + "These will be assigned 'nan' values\nfor masses and "
                + "binarity probability"
            )

        self._vp("\nBinary fraction estimated", 1)

        return data

    def binary_fraction(
        self,
    ) -> np.ndarray:
        """Estimate the total binary fraction for the observed cluster.

        :return: Distribution of total binary fraction values for the cluster
        :rtype: np.ndarray
        """
        from .modules import mass_binary as mb

        b_fr_all = []
        for i, isoch in enumerate(self.sampled_synthcls):
            b_fr = mb.get_bpr(isoch, self.close_stars_idxs[i])
            b_fr_all.append(b_fr)

        return np.array(b_fr_all)

    def cluster_masses(
        self,
        radec_c: tuple[float, float] | None = None,
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

        :param radec_c: Right ascension and declination center coordinates for the
            cluster. If ``None`` then the ``rho_amb`` must be given; defaults to
            ``None``.
        :type radec_c: tuple[float, float] | None
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
            -if no synthetic models were generated.
            - if both ``radec_c`` and ``rho_amb`` are ``None``.

        :return: Dictionary with the mass distributions for the initial, actual,
            observed, photometric, evolutionary, and dynamical masses:
            ``M_init, M_actual, M_obs, M_phot, M_evol, M_dyn``
        :rtype: dict
        """
        from .modules import mass_binary as mb

        if len(self.sampled_models) == 0:
            raise ValueError(
                "No models were generated. Run the `get_models()` method before "
                + "estimating the cluster masses."
            )

        # Number of observed stars
        N_obs = len(self.mag_v)
        # The number of stars in a synthetic isochrones is not constant so we estimate
        # its median
        N_stars_isoch = int(np.median([np.shape(_)[-1] for _ in self.sampled_synthcls]))
        # Compare the number of observed vs generated synthetic stars
        if N_stars_isoch < N_obs:
            warnings.warn(
                f"\nNumber of synthetic stars ({N_stars_isoch}) is smaller than "
                + f"observed stars ({N_obs}). Increase the 'max_mass' argument for "
                + "a more accurate mass estimation"
            )

        if rho_amb is not None:
            # Input float to array
            rho_amb_arr = np.ones(len(self.sampled_models)) * rho_amb
            if radec_c is not None:
                warnings.warn(
                    "Both 'radec_c' and 'rho_amb' were given. Using 'rho_amb' only."
                )
        else:
            if radec_c is not None:
                # Obtain galactic vertical distance and distance to center
                Z, R_GC, R_xy = mb.galactic_coords(
                    self.sampled_models, self.fix_params, radec_c
                )
                rho_amb_arr = mb.ambient_density(
                    M_B, r_B, M_D, a, b, r_s, M_s, Z, R_GC, R_xy
                )
            else:
                raise ValueError(
                    "Either the 'radec_c' or 'rho_amb' arguments must be given."
                )

        masses_all = []
        for i, model in enumerate(self.sampled_models):
            # Extract met and loga
            model_comb = self.fix_params | model
            z_met, loga = model_comb["met"], model_comb["loga"]
            if self.isochs.z_to_FeH is not None:
                z_met = self.isochs.z_to_FeH * 10**z_met

            # Estimate the actual mass, ie: the sum of the observed and photometric
            # masses
            sampled_synthcl = self.sampled_synthcls[i]
            # TODO: don't pass the entire object, just the necessary values
            M_obs, M_phot = mb.get_M_actual(
                self.rng,
                self.m_ini_idx,
                self.st_dist_mass,
                self.st_dist_mass_ordered,
                sampled_synthcl,
            )
            M_a = M_obs + M_phot

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

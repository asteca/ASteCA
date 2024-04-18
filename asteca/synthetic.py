import warnings
from dataclasses import dataclass
from typing import Optional
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from .isochrones import isochrones
from .modules import synth_cluster_priv as scp
from .modules import mass_binary as mb


@dataclass
class synthetic:
    r"""Define a ``synthetic`` object.

    Use the isochrones loaded in the :py:mod:`asteca.isochrones` object to generate a
    :py:mod:`asteca.synthetic` object. This object is used to generate synthetic clusters
    given a :py:mod:`asteca.cluster` object  and a set of input fundamental parameters
    (metallicity, age, distance, extinction, etc.).

    See the :ref:`synth_clusters` section for more details.

    Parameters
    ----------
    isochs : :class:`isochrones`
         :py:mod:`asteca.isochrones` object with the loaded files for the theoretical isochrones.
    IMF_name : str, {"salpeter_1955", "kroupa_2001", "chabrier_2014"}, default="chabrier_2014"
        Name of the initial mass function used to populate the isochrones.
    max_mass : int, default=100_000
        Maximum total initial mass. Should be large enough to allow generating as many
        synthetic stars as observed stars.
    gamma : str, float, {"D&K", "fisher_stepped", "fisher_peaked", "raghavan"}, default="D&K"
        Distribution function for the mass ratio of the binary systems.
    DR_distribution : str, {"uniform", "normal"}, default="uniform"
        Distribution function for the differential reddening.
    seed: int, optional, default=None
        Random seed. If ``None`` a random integer will be generated and used.

    """
    isochs: isochrones
    IMF_name: str = "chabrier_2014"
    max_mass: int = 100_000
    gamma: float | str = "D&K"
    DR_distribution: str = "uniform"
    seed: Optional[int] = None

    def __post_init__(self):

        if self.seed is None:
            self.seed = np.random.randint(100000)
        print(f"Setting random seed to {self.seed}")

        # Check gamma distribution
        gammas = (
            "D&K",
            "fisher_stepped",
            "fisher_peaked",
            "raghavan",
        )
        if isinstance(self.gamma, str):
            if self.gamma not in gammas:
                raise ValueError(
                    f"gamma '{self.gamma}' not recognized. Should be one of {gammas}"
                )

        # Check differential reddening function
        DR_funcs = (
            "uniform",
            "normal",
        )
        if self.DR_distribution not in DR_funcs:
            raise ValueError(
                f"Differential reddening function '{self.DR_distribution}' not recognized. "
                + f"Should be one of {DR_funcs}"
            )

        # Check IMF function
        imfs = (
            "salpeter_1955",
            "kroupa_2001",
            "chabrier_2014",
        )
        if self.IMF_name not in imfs:
            raise ValueError(
                f"IMF '{self.IMF_name}' not recognized. Should be one of {imfs}"
            )

        # Sample the selected IMF
        # Nmets, Nages = self.theor_tracks.shape[:2]
        Nmets, Nages = self.isochs.theor_tracks.shape[:2]
        self.st_dist_mass = scp.sample_imf(self, Nmets, Nages)

        # Add binary systems
        self.theor_tracks = scp.add_binarity(self)

        # Get extinction coefficients for these filters
        self.ext_coefs = scp.extinction_coeffs(self)

        # Generate random floats used by `synth_clusters.synthcl_generate()`
        self.rand_floats = scp.randVals(self)

        print("Synthetic clusters object generated\n")

    def calibrate(self, cluster, fix_params: dict = {}, z_to_FeH: float | None = None):
        r"""Calibrate a :py:mod:`asteca.synthetic` object based on a
        :py:mod:`asteca.cluster` object and a dictionary of fixed fundamental parameters
        (``fix_params``).

        Use the data obtained from your observed cluster stored in the
        :py:mod:`asteca.cluster` object, to calibrate a :py:mod:`asteca.synthetic`
        object. Additionally, a dictionary of fixed fundamental parameters 
        (metallicity, age, distance, extinction, etc.) can be passed.

        See the :ref:`synth_clusters` section for more details.

        Parameters
        ----------
        cluster : :class:`cluster`
             :py:mod:`asteca.cluster` object with the processed data from your observed
             cluster.
        fix_params : dict, optional, default={}
            Dictionary with the values for the fixed parameters (if any).
        z_to_FeH : float, optional, default=None
            If ``None``, the default ``z`` values (defined when loading the isochrones
            via the :py:mod:`asteca.isochrones` object) will be used to generate the
            synthetic clusters. If ``float``, it must represent the solar metallicity
            for these isochrones. The metallicity values will then be converted to
            ``[FeH]`` values, to be used by the :meth:`synthetic.generate()` method.

        """
        # Used by the mass and binary probability estimation
        self.mag_p = cluster.mag_p
        self.colors_p = cluster.colors_p

        # Data used by the `generate()` method
        self.max_mag_syn = max(cluster.mag_p)
        self.N_obs_stars = len(cluster.mag_p)
        self.m_ini_idx = len(cluster.colors_p) + 1
        self.err_dist = scp.error_distribution(
            self, cluster.mag_p, cluster.e_mag_p, cluster.e_colors_p
        )

        self.fix_params = fix_params

        # Convert z to FeH if requested
        self.met_age_dict = self.isochs.met_age_dict
        if z_to_FeH is not None:
            self._func_z_to_FeH(z_to_FeH)

        # Check if binary systems are to be produced
        self.binar_flag = True
        if "alpha" in list(fix_params.keys()) and "beta" in list(fix_params.keys()):
            if fix_params["alpha"] == 0.0 and fix_params["beta"] == 0.0:
                self.binar_flag = False

        # # Remove low masses if required
        # if dm_min is not None:
        #     self._rm_low_masses(dm_min)

    def _func_z_to_FeH(self, z_to_FeH):
        """ """
        feh = np.log10(self.met_age_dict["met"] / z_to_FeH)
        N_old = len(feh)
        round_n = 4
        while True:
            feh_r = np.round(feh, round_n)
            N_new = len(set(feh_r))
            # If no duplicated values exist after rounding
            if N_old == N_new:
                break
            round_n += 1
        # Replace old values
        self.met_age_dict["met"] = feh_r

    def min_max(self) -> tuple[float]:
        r"""Return the minimum and maximum values for the metallicity and age defined
        in the theoretical isochrones.

        Returns
        -------
        tuple[float]
            Tuple of (minimum_metallicity, maximum_metallicity, minimum_age, maximum_age)

        """
        zmin = self.met_age_dict["met"].min()
        zmax = self.met_age_dict["met"].max()
        amin = self.met_age_dict["a"].min()
        amax = self.met_age_dict["a"].max()
        return zmin, zmax, amin, amax

    def generate(self, fit_params: dict) -> np.ndarray:
        r"""Generate a synthetic cluster.

        The synthetic cluster is generated according to the parameters given in
        the ``fit_params`` dictionary and the already calibrated
        :py:mod:`asteca.synthetic` object.

        Parameters
        ----------
        fit_params : dict
            Dictionary with the values for the fundamental parameters that were **not**
            included in the ``fix_params`` dictionary when the
            :py:mod:`asteca.synthetic` object was calibrated 
            (:meth:`synthetic.calibrate()` method).

        Returns
        -------
        array[mag, c1, (c2)]
            Return a ``np.array`` containing a synthetic cluster with the shape
            ``[mag, c1, (c2)]``, where ``mag`` is the magnitude dimension, and
            ``c1`` and ``c2`` (last one is optional) are the color dimension(s).

        """

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
            self.ext_coefs,
            self.rand_floats["norm"][0],
            self.rand_floats["unif"][0],
            self.DR_distribution,
            self.m_ini_idx,
            av,
            dr,
            rv,
            isoch_moved,
        )

        # Remove isochrone stars beyond the maximum magnitude
        isoch_cut = scp.cut_max_mag(isoch_extin, self.max_mag_syn)
        if not isoch_cut.any():
            return np.array([])

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
            self.gamma,
            self.m_ini_idx,
            self.rand_floats["unif"][1],
            isoch_mass,
        )

        # Assign errors according to errors distribution.
        synth_clust = scp.add_errors(
            isoch_binar, self.err_dist, self.rand_floats["norm"][1]
        )

        return synth_clust[: self.m_ini_idx]

    def synthplot(self, fit_params, ax=None, isochplot=False):
        r"""Generate a color-magnitude plot for a synthetic cluster.

        The synthetic cluster is generated using the fundamental parameter values
        given in the ``fit_params`` dictionary.

        Parameters
        ----------
        fit_params : dict
            Dictionary with the values for the fundamental parameters that were **not**
            included in the ``fix_params`` dictionary when the
            :py:mod:`asteca.synthetic` object was calibrated
            (:meth:`synthetic.calibrate()` method).
        ax : matplotlib.axis, optional, default=None
            Matplotlib axis where to draw the plot.
        isochplot : bool, default=False
            If ``True``, the accompanying isochrone will be plotted.

        Returns
        -------
        matplotlib.axis
            Matplotlib axis object

        """
        if ax is None:
            f, ax = plt.subplots()

        # Generate synthetic cluster.
        synth_clust = scp.generate(self, fit_params)
        if self.binar_flag is True:
            binar_idx = ~np.isnan(synth_clust[-1])
        else:
            binar_idx = np.full(synth_clust.shape[1], False)

        x_synth, y_synth = synth_clust[1], synth_clust[0]
        # Single synthetic systems
        ax.scatter(
            x_synth[~binar_idx],
            y_synth[~binar_idx],
            marker="^",
            c="#519ddb",
            alpha=0.5,
            label=f"Synthetic (single), N={len(x_synth[~binar_idx])}",
        )
        # Binary synthetic systems
        ax.scatter(
            x_synth[binar_idx],
            y_synth[binar_idx],
            marker="v",
            c="#F34C4C",
            alpha=0.5,
            label=f"Synthetic (binary), N={len(x_synth[binar_idx])}",
        )

        plt.ylabel(list(self.isochs.magnitude.keys())[0])
        c1, c2 = list(self.isochs.color.keys())
        plt.xlabel(f"{c1}-{c2}")
        ax.set_ylim(max(self.mag_p) + .5, min(self.mag_p) - .5)
        ax.legend()

        if isochplot is False:
            return ax

        # Generate displaced isochrone
        fit_params_copy = dict(fit_params)
        fit_params_copy["DR"] = 0.0
        isochrone = scp.generate(self, fit_params_copy, True)
        # Remove stars beyond the color limits
        xmin, xmax = x_synth[~binar_idx].min(), x_synth[~binar_idx].max()
        msk = (isochrone[1] >= xmin) & (isochrone[1] <= xmax)
        isochrone = isochrone[:, msk]
        ax.plot(isochrone[1], isochrone[0], c="k")

        return ax

    def masses_binary_probs(self, model, model_std):
        r"""Estimate individual masses for the observed stars, along with their binary
        probabilities (if binarity was estimated).

        Parameters
        ----------
        model : dict
            Dictionary with the values for the fundamental parameters that were **not**
            included in the ``fix_params`` dictionary when the
            :py:mod:`asteca.synthetic` object was calibrated 
            (:meth:`synthetic.calibrate()` method).
        model_std : dict
            Dictionary with the standard deviations for the fundamental parameters in
            the ``model`` argument.

        Returns
        -------
        pandas.DataFrame
            Data frame containing per-star primary and secondary masses along with
            their uncertainties, and their probability of being a binary system.
        numpy.array
            Distribution of total binary fraction values for the cluster.

        """
        # Generate random models from the selected solution
        models = mb.ranModels(model, model_std, self.seed)

        # Observed photometry
        obs_phot = np.array([self.mag_p] + [_ for _ in self.colors_p])
        # Identify nans in mag and colors and re-write them as -10. This ensures
        # that the photometric distances to these stars is large enough to not be
        # matched.
        nan_msk = np.full(obs_phot.shape[1], False)
        for col in obs_phot:
            nan_msk = nan_msk | np.isnan(col)
        obs_phot[:, nan_msk] = -10.0
        obs_phot = obs_phot.T

        m12_masses, b_fr_all,  = [], []
        for model in models:
            isoch = scp.generate(self, model)
            if not isoch.any():
                continue
            m1_obs, m2_obs, b_fr = mb.get_m1m2_bpr(self, isoch, obs_phot)
            m12_masses.append([m1_obs, m2_obs])
            b_fr_all.append(b_fr)
        m12_masses = np.array(m12_masses)

        # Primary masses (median + stddev)
        m1_med = np.median(m12_masses[:, 0, :], 0)
        m1_std = np.std(m12_masses[:, 0, :], 0)
        # Secondary masses  (median + stddev). Hide 'All-nan slice' warnings
        with warnings.catch_warnings(action="ignore"):
            m2_med = np.nanmedian(m12_masses[:, 1, :], 0)
            # m2 can not be larger than m1
            m2_med = np.min([m1_med, m2_med], 0)
            m2_std = np.nanstd(m12_masses[:, 1, :], 0)

        # Binary probability per star
        binar_prob = (~np.isnan(m12_masses[:, 1, :])).sum(0) / m12_masses.shape[0]

        # Store as pandas.DataFrame
        df = pd.DataFrame(
            {
                "m1": m1_med,
                "m1_std": m1_std,
                "m2": m2_med,
                "m2_std": m2_std,
                "binar_prob": binar_prob,
            }
        )

        return df, np.array(b_fr_all)

    def _get_masses(self, fit_params, model_std, ra_c, dec_c):
        """
        Estimate the different total masses for the observed cluster
        """
        print("Estimating total initial and actual masses")
        # Generate random models from the selected solution
        models = mb.ranModels(fit_params, model_std, self.seed)

        masses_all = []
        for i, model in enumerate(models):
            isoch = scp.generate(self, model)
            if not isoch.any():
                continue
            masses_all.append(mb.get_masses(self, model, ra_c, dec_c, isoch, i))
        masses_all = np.array(masses_all).T

        # Check the number of generated synthetic stars
        N_obs = len(self.mag_p)
        if N_obs > isoch.shape[1]:
            warnings.warn(
                "Number of synthetic stars is smaller than observed stars. Increase "
                + "the 'max_mass' argument for a more accurate mass estimation"
            )

        return {"M_init": masses_all[0], "M_actual": masses_all[1]}

    # def _rm_low_masses(self, dm_min):
    #     """
    #     dm_min: float | None = None

    #     dm_min : float, optional, default=None
    #         Value for the minimum distance modulus. Used to constrain the lower masses
    #         in the theoretical isochrones to make the generating process more
    #         efficient.
    #     """
    #     min_masses = []
    #     for met_arr in self.theor_tracks:
    #         met_lst = []
    #         for age_arr in met_arr:
    #             mag, mass = age_arr[0], age_arr[self.m_ini_idx]
    #             i = np.argmin(abs(self.max_mag_syn - (mag + dm_min)))
    #             met_lst.append(mass[i])
    #         min_masses.append(met_lst)

    #     st_dist_mass_lmass = []
    #     for i, met_arr in enumerate(self.st_dist_mass):
    #         met_lst = []
    #         for j, mass_sample in enumerate(met_arr):
    #             min_mass = min_masses[i][j]
    #             msk = mass_sample > min_mass
    #             sampled_IMF = mass_sample[msk]
    #             met_lst.append(sampled_IMF)
    #         st_dist_mass_lmass.append(met_lst)

    #     # Make copy of original array, used for mass estimation in cluster() class
    #     self.st_dist_mass_full = self.st_dist_mass.copy()
    #     # Update this parameter with the new array
    #     self.st_dist_mass = st_dist_mass_lmass

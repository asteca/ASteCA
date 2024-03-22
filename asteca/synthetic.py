from dataclasses import dataclass
import numpy as np
from .isochrones import isochrones
from .modules import synth_cluster_priv as scp


@dataclass
class synthetic:
    r"""Define a ``synthetic`` object.

    Use the data loaded in the :class:`isochrones` object to generate a
    :class:`synthetic` object. This object is used to generate synthetic clusters
    given a set of input fundamental parameters (metallicity, age, distance,
    extinction, etc.).

    Parameters
    ----------
    isochs : :class:`isochrones`
         :doc:`isochrones` object with the loaded files for the theoretical isochrones.
    alpha : None, float
        First parameter of the 
    gamma : str, float
        xxxxx
    DR_dist : str
        xxxxx
    DR_percentage : float
        xxxxx
    IMF_name : str
        xxxxx
    max_mass : int
        xxxxx

    """
    isochs: isochrones
    alpha: None | float = .0
    gamma: float | str = "D&K"
    DR_dist: str = "uniform"
    DR_percentage: float = 1.0
    IMF_name: str = "kroupa_2001"
    max_mass: int = 100_000

    def __post_init__(self):
        # Check gamma distribution
        gammas = (
            "D&K",
            "fisher_stepped",
            "fisher_peaked",
            "raghavan",
        )
        if type(self.gamma) is str:
            if self.gamma not in gammas:
                raise ValueError(
                    f"gamma '{self.gamma}' not recognized. Should be one of {gammas}"
                )

        # Check differential reddening function
        DR_funcs = (
            "uniform",
            "normal",
        )
        if self.DR_dist not in DR_funcs:
            raise ValueError(
                f"Differential reddening function '{self.DR_dist}' not recognized. "
                + f"Should be one of {DR_funcs}"
            )

        # Check IMF function
        imfs = (
            "salpeter_1955",
            "kroupa_1993",
            "kroupa_2001",
            "chabrier_2001_log",
            "chabrier_2001_exp",
            "popescu_2009",
        )
        if self.IMF_name not in imfs:
            raise ValueError(
                f"IMF '{self.IMF_name}' not recognized. Should be one of {imfs}"
            )

        # Add binary systems
        self.theor_tracks = scp.add_binarity(self)

        # Get extinction coefficients for these filters
        self.ext_coefs = scp.extinction_coeffs(self)

        # Sample the selected IMF
        Nmets, Nages = self.theor_tracks.shape[:2]
        self.st_dist_mass = scp.sample_imf(self, Nmets, Nages)

        # Generate random floats used by `synth_clusters.synthcl_generate()`
        self.rand_floats = scp.randVals(self)
        print("Finished generating parameters for synthetic clusters")

    def adjust_mass_sampling(self, my_cluster, dm_min: float) -> tuple[list, list]:
        """ """
        min_masses = []
        for met_arr in self.theor_tracks:
            met_lst = []
            for age_arr in met_arr:
                mag, mass = age_arr[0], age_arr[my_cluster.cluster_dict['m_ini_idx']]
                i = np.argmin(abs(my_cluster.cluster_dict['max_mag_syn'] - (mag + dm_min)))
                met_lst.append(mass[i])
            min_masses.append(met_lst)

        st_dist_mass_lmass = []
        for i, met_arr in enumerate(self.st_dist_mass):
            met_lst = []
            for j, mass_sample in enumerate(met_arr):
                min_mass = min_masses[i][j]
                msk = mass_sample > min_mass
                sampled_IMF = mass_sample[msk]
                met_lst.append(sampled_IMF)
            st_dist_mass_lmass.append(met_lst)

        # Make copy of original array, used for mass estimation in cluster() class
        self.st_dist_mass_full = self.st_dist_mass.copy()
        # Update this parameter with the new array
        self.st_dist_mass = st_dist_mass_lmass

    def generate(self, my_cluster, model_fit: dict) -> np.ndarray:
        r"""Return a synthetic cluster.

        The synthetic cluster is generated according to the parameters given in
        the ``model_fixed`` and ``model_fit`` dictionaries.

        The synthetic cluster returned has the shape:

        synth_clust = [mag, c1, (c2)]

        where c1 and c2 colors defined.

        """

        # Return proper values for fixed parameters and parameters required
        # for the (z, log(age)) isochrone averaging.
        met, loga, beta, av, dr, rv, dm, ml, mh, al, ah = scp.properModel(
            self.isochs.met_age_dict, my_cluster.model_fixed, model_fit
        )

        # Generate a weighted average isochrone from the (z, log(age)) values in
        # the 'model'.
        isochrone = scp.zaWAverage(
            self.theor_tracks,
            self.isochs.met_age_dict,
            my_cluster.cluster_dict['m_ini_idx'],
            met,
            loga,
            ml,
            mh,
            al,
            ah,
        )

        # Move theoretical isochrone using the distance modulus
        isoch_moved = scp.move_isochrone(
            isochrone, my_cluster.cluster_dict['m_ini_idx'], dm)

        # Apply extinction correction
        isoch_extin = scp.extinction(
            self.ext_coefs,
            self.rand_floats['norm'][0],
            self.rand_floats['unif'][0],
            self.DR_dist,
            self.DR_percentage,
            my_cluster.cluster_dict['m_ini_idx'],
            av,
            dr,
            rv,
            isoch_moved,
        )

        # Remove isochrone stars beyond the maximum magnitude
        isoch_cut = scp.cut_max_mag(
            isoch_extin, my_cluster.cluster_dict["max_mag_syn"])
        if not isoch_cut.any():
            return np.array([])

        # Interpolate IMF's sampled masses into the isochrone.
        isoch_mass = scp.mass_interp(
            isoch_cut, my_cluster.cluster_dict['m_ini_idx'], self.st_dist_mass[ml][al],
            my_cluster.cluster_dict["N_obs_stars"])
        if not isoch_mass.any():
            return np.array([])

        # Assignment of binarity.
        isoch_binar = scp.binarity(
            self.alpha, beta, my_cluster.cluster_dict['m_ini_idx'],
            self.rand_floats['unif'][1], isoch_mass
        )

        # Assign errors according to errors distribution.
        synth_clust = scp.add_errors(
            isoch_binar, my_cluster.cluster_dict["err_lst"],
            self.rand_floats['norm'][1]
        )

        return synth_clust[:my_cluster.cluster_dict['m_ini_idx']]

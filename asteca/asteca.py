import os
import numpy as np
from dataclasses import dataclass
from typing import Optional
#
from .cluster import cluster_load, get_masses_binar, mplot
from .isochrones import isochs_load
from .synth_cluster import add_binarity, extinction_coeffs, sample_imf, randVals, \
    synthcl_generate, remove_low_masses
from .likelihood import get_lkl


@dataclass
class isochrones:
    r"""Define an `isochrones` object.

    The object contains the loaded theoretical isochrones used by the
    `synth_clusters` class to generate synthetic clusters.

    Parameters
    ----------
    isochs_path : str
        Path to the folder that contains the files for the theoretical isochrones.
        The name of the folder must be one of the supported isochrone services:
        [PARSEC][1], [MIST][2], or [BASTI][3]. Examples of valid paths:
        'isochrones/PARSEC/', 'Mist/', 'basti'.
        See [REFERENCE][TODO] on how to properly store the isochrone files.
    magnitude : dict
        Dictionary containing the magnitude's filter name (as defined in the files
        of the theoretical isochrones) as the key, and its effective lambda
        (in Angstrom) as the value. Example for Gaia's 'G' magnitude:
        `{"Gmag": 6390.21}`
    color : dict
        Dictionary containing the color used in the cluster's analysis. The correct
        format is: `{"filter1": 1111.11, "filter2": 2222.22}`, where `filter1`
        and `filter2` are the names of the filters that are combined to generate
        the color. The order is important because the color will be generated as:
        `filter1-filter2``. The values `1111.11` and ``2222.22` are the effective
        lambdas (in Angstrom) for each filter. The color does not need to be defined
        in the same photometric system as the magnitude.
        Example for Gaia's 'BP-RP' color: `{"G_BPmag": 5182.58, "G_RPmag": 7825.08}`
    color2 : dict, optional
        Optional second color to use in the analysis. Same format as that used by the
        `color` parameter.
    column_names : dict
        Column names for the initial mass, metallicity, and age for the photometric
        system's isochrones files. Example:
        `{"mass_col": "Mini", "met_col": "Zini", "age_col": "logAge"}`.
        This dictionary is defined internally in `ASteCA` and should only be given
        by the user if the isochrone service changes its format and the `isochrones`
        class fails to load the files.
    N_interp : int
        Number of interpolation points used to ensure that all isochrones are the
        same shape.

    References
    ----------
    .. [1] http://stev.oapd.inaf.it/cgi-bin/cmd_3.7
    .. [2] https://waps.cfa.harvard.edu/MIST/
    .. [3] http://basti-iac.oa-abruzzo.inaf.it/isocs.html
    .. [TODO] url_here

    """
    isochs_path: str
    magnitude: dict
    color: dict
    color2: Optional[dict] = None
    column_names: Optional[dict] = None
    N_interp: int = 2500

    def __post_init__(self):
        # Extract model from isochrones' path
        if self.isochs_path.endswith("/"):
            self.model = self.isochs_path[:-1].split("/")[-1].upper()
        else:
            self.model = self.isochs_path.split("/")[-1].upper()

        # Check model input
        models = ("PARSEC", "MIST", "BASTI")
        if self.model not in models:
            raise ValueError(
                f"Model '{self.model}' not recognized. Should be one of {models}"
            )

        # Check path to isochrones
        if os.path.isdir(self.isochs_path) is False:
            raise ValueError(
                f"Path '{self.isochs_path}' must point to the folder that "
                + f"contains the '{self.model}' isochrones"
            )

        # Extract magnitude, color(s), and lambdas
        self.mag_color_lambdas = self.magnitude | self.color
        self.mag_filter_name = list(self.magnitude.keys())[0]
        self.color_filter_name = [list(self.color.keys())]
        # Add second color if any
        if self.color2 is not None:
            self.color_filter_name.append(list(self.color2.keys()))
            self.mag_color_lambdas = self.mag_color_lambdas | self.color2

        # Load isochrone files
        self.theor_tracks, self.color_filters, self.met_age_dict = isochs_load(self)

    def min_max(self) -> tuple[float]:
        """Return the minimum and maximum values for the metallicity and age defined
        in the theoretical isochrones.
        """
        zmin = self.met_age_dict['z'].min()
        zmax = self.met_age_dict['z'].max()
        amin = self.met_age_dict['a'].min()
        amax = self.met_age_dict['a'].max()
        return zmin, zmax, amin, amax


@dataclass
class synth_clusters:
    r"""Define a `synth_clusters` object.

    Use the data loaded in the `isochrones` object to generate a `synth_clusters`
    object, used to generate synthetic clusters.

    Parameters
    ----------
    isochs : isochrones
        xxxxx
    alpha : xxx
        xxxxx
    gamma : xxx
        xxxxx
    DR_dist : xxx
        xxxxx
    DR_percentage : xxx
        xxxxx
    IMF_name : xxx
        xxxxx
    max_mass : xxx
        xxxxx

    """
    isochs: isochrones
    alpha: float = 0.45
    gamma: float | str = "D&K"
    DR_dist: str = "uniform"
    DR_percentage: float = 1
    IMF_name: str = "kroupa_2001"
    max_mass: float = 1e5

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

        # Store this parameter as part of the `synth_clusters` class
        self.met_age_dict = self.isochs.met_age_dict

        # Add binary systems
        self.theor_tracks = add_binarity(self)

        # Get extinction coefficients for these filters
        self.ext_coefs = extinction_coeffs(self)

        # Sample the selected IMF
        Nmets, Nages = self.theor_tracks.shape[:2]
        self.st_dist_mass = sample_imf(self, Nmets, Nages)

        # Generate random floats used by `synth_clusters.synthcl_generate()`
        self.rand_floats = randVals(self.theor_tracks, self.st_dist_mass)

    def adjust_mass_sampling(self, my_cluster, dm_min):
        """ """
        self.st_dist_mass_full, self.st_dist_mass = remove_low_masses(
            self, my_cluster, dm_min)

    def generate(self, model_fit, my_cluster, full_synth_arr=False):
        """ """
        return synthcl_generate(
            self, model_fit, my_cluster.model_fixed, my_cluster.cluster_dict,
            full_synth_arr
        )


@dataclass
class cluster:
    source_id: str
    magnitude: str
    e_mag: str
    color: str
    e_color: str
    model_fixed: dict
    cluster_df: Optional[None] = None
    ra: Optional[str] = None
    dec: Optional[str] = None
    bin_method: str = "knuth"
    N_mag: int = 15
    N_col: int = 10
    N_models: int = 100
    color2: str | None = None
    e_color2: str | None = None

    def __post_init__(self):
        bin_methods = ("knuth", "fixed", "bayes_blocks")
        if self.bin_method not in bin_methods:
            raise ValueError(
                f"Binning '{self.bin_method}' not recognized. "
                + f"Should be one of {bin_methods}"
            )

        self.cluster_dict = cluster_load(self)

    def masses_binar(self, synthcl, model_fit, model_std):
        """ """
        return get_masses_binar(self, synthcl, model_fit, model_std)

    def clust_plot(self, synthcl, model_fit):
        mplot(self, synthcl, model_fit)


@dataclass
class likelihood:
    lkl_name: str = "plr"

    def __post_init__(self):
        likelihoods = ("plr", "bins_distance")
        if self.lkl_name not in likelihoods:
            raise ValueError(
                f"'{self.lkl_name}' not recognized. Should be one of {likelihoods}"
            )

    def get(self, my_cluster, synth):
        return get_lkl(self.lkl_name, my_cluster.cluster_dict, synth)

    def get_optm(self, my_cluster):
        cl_dict = my_cluster.cluster_dict
        obs_cluster = np.array([cl_dict["mag"], *cl_dict["colors"]])
        return get_lkl(self.lkl_name, cl_dict, obs_cluster)

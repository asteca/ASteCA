import os
import numpy as np
from dataclasses import dataclass
from typing import Tuple, Dict

#
from .cluster import cluster_load, get_masses_binar, mplot
from .isochrones import isochs_load, isochs_minmax
from .synth_cluster import synthcl_load, synthcl_generate
from .likelihood import get_lkl


@dataclass
class cluster:
    id_column: str
    mag_column: str
    e_mag_column: str
    col_column: Tuple
    e_col_column: Tuple
    cluster_df: None
    ra: float
    dec: float
    model_fixed: dict
    bin_method: str = "knuth"
    N_mag: int = 15
    N_col: int = 10

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
        return get_masses_binar(self, synthcl, model_fit,model_std)

    def clust_plot(self, synthcl, model_fit):
        mplot(self, synthcl, model_fit)


@dataclass
class isochrones:
    magnitude: Dict
    color: Dict
    # mag_color_lambdas: Dict
    isochs_path: str
    color2: Dict | None = None
    metallicity_col: str | None = None
    age_col: str | None = None
    initial_mass_col: str | None = None
    N_interp: int = 2500

    def __post_init__(self):
        # Extract model from isochrones' path
        if self.isochs_path.endswith("/"):
            self.model = self.isochs_path[:-1].split("/")[-1].upper()
        else:
            self.model = self.isochs_path.split("/")[-1].upper()

        models = ("PARSEC", "MIST", "BASTI")
        if self.model not in models:
            raise ValueError(
                f"Model '{self.model}' not recognized. Should be one of {models}"
            )

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

        self.isochs_dict = isochs_load(self)

    def min_max(self):
        """ """
        return isochs_minmax(self.isochs_dict)


@dataclass
class synth_clusters:
    isochs: None
    gamma: str = "D&K"
    DR_dist: str = "uniform"
    DR_percentage: float = 1
    IMF_name: str = "kroupa_2001"
    max_mass: float = 8000
    completeness_f: Tuple = None
    binarity: bool = True

    def __post_init__(self):
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

        # def load(self, isochs_dict):
        self.synthcl_dict = synthcl_load(self)

    def generate(
        self, model_fit, my_cluster, M_total_flag=False
    ):
        return synthcl_generate(
            self, model_fit, my_cluster.model_fixed, my_cluster.cluster_dict, M_total_flag
        )


@dataclass
class likelihood:
    lkl_name: str = "plr"

    def __post_init__(self):
        likelihoods = ("plr",)
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

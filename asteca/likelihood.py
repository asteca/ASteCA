from dataclasses import dataclass
import numpy as np
from .modules import likelihood_priv


@dataclass
class likelihood:
    r""" """
    lkl_name: str = "plr"

    def __post_init__(self):
        likelihoods = ("plr", "bins_distance")
        if self.lkl_name not in likelihoods:
            raise ValueError(
                f"'{self.lkl_name}' not recognized. Should be one of {likelihoods}"
            )

    def get(self, my_cluster, synth_clust=None):
        """ """
        if synth_clust is None:
            synth_clust = np.array([my_cluster.cluster_dict["mag"], *my_cluster.cluster_dict["colors"]])

        if self.lkl_name == "plr":
            return likelihood_priv.tremmel(my_cluster.cluster_dict, synth_clust)
        if self.lkl_name == "visual":
            return likelihood_priv.visual(my_cluster.cluster_dict, synth_clust)
        if self.lkl_name == "mean_dist":
            return likelihood_priv.mean_dist(my_cluster.cluster_dict, synth_clust)
        if self.lkl_name == "bins_distance":
            return likelihood_priv.bins_distance(my_cluster.cluster_dict, synth_clust)

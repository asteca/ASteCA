from dataclasses import dataclass
import numpy as np
from .cluster import cluster
from .modules import likelihood_priv as lpriv


@dataclass
class likelihood:
    r"""Define a ``likelihood`` object.

    This object is used to assess how similar your observed cluster is, stored in
    a :py:class:`cluster` object,  compared to a given synthetic cluster, generated by
    the :py:meth:`synthetic.generate()` method.

    Parameters
    ----------
    my_cluster : :class:`cluster`
         :py:mod:`asteca.cluster` object with the loaded data for the observed cluster.
    lkl_name : str, {"plr"}, default="plr"
        Currently only the Poisson likelihood ratio defined in
        `Tremmel et al. (2013) <https://ui.adsabs.harvard.edu/abs/2013ApJ...766...19T/abstract)>`_
        is accepted.
    bin_method: str, {"knuth", "fixed", "bayes_blocks", "manual"}, default="knuth"
        Bin method used to split the color-magnitude diagram into cells
        (`Hess diagram <https://en.wikipedia.org/wiki/Hess_diagram>`_). If ``manual``
        is selected, a list containing an array of edge values for the magnitude,
        followed by one or two arrays (depending on the number of colors defined) for
        the color(s), also with edge values.

    """
    my_cluster: cluster
    lkl_name: str = "plr"
    bin_method: str = "knuth"

    def __post_init__(self):
        likelihoods = ("plr", "bins_distance")
        if self.lkl_name not in likelihoods:
            raise ValueError(
                f"'{self.lkl_name}' not recognized. Should be one of {likelihoods}"
            )

        bin_methods = ("knuth", "fixed", "bayes_blocks", "manual")
        if self.bin_method not in bin_methods:
            raise ValueError(
                f"Binning '{self.bin_method}' not recognized. "
                + f"Should be one of {bin_methods}"
            )

        # Obtain data used by the ``likelihood.get()`` method
        lpriv.lkl_data(self)

        # Evaluate cluster against itself
        self.max_lkl = self.get(
            np.array([self.my_cluster.mag_p, *self.my_cluster.colors_p])
        )

        print("Likelihood object generated\n")

    def get(self, synth_clust):
        r"""Evaluate the selected likelihood function.

        Parameters
        ----------
        synth_clust : array
            ``np.array`` containing the synthetic cluster. The shape of this array must
            be: ``[magnitude, color1, (color2)]``, where ``magnitude`` and ``color``
            are arrays with the magnitude and color photometric data (``color2`` is the
            optional second color defined).

        Returns
        -------
        float
            Likelihood value.

        """
        if self.lkl_name == "plr":
            return lpriv.tremmel(self, synth_clust)
        if self.lkl_name == "visual":
            return lpriv.visual(self, synth_clust)
        if self.lkl_name == "mean_dist":
            return lpriv.mean_dist(self, synth_clust)
        if self.lkl_name == "bins_distance":
            return lpriv.bins_distance(self, synth_clust)
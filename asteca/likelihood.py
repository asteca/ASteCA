import numpy as np

from .cluster import Cluster
from .modules import likelihood_priv as lpriv


class Likelihood:
    """Define a :py:class:`Likelihood` object.

    This object is used to assess how similar your observed cluster is, stored in
    a :py:class:`Cluster <asteca.cluster.Cluster>` object,  compared to a given
    synthetic cluster, generated by the
    :py:meth:`Synthetic.generate() <asteca.synthetic.Synthetic.generate>` method.

    :param my_cluster: :py:class:`Cluster <asteca.cluster.Cluster>` object with the
        loaded data for the observed cluster
    :type my_cluster: Cluster
    :param lkl_name: Currently only the Poisson likelihood ratio (``plr``) defined in
        `Tremmel et al. (2013)
        <https://ui.adsabs.harvard.edu/abs/2013ApJ...766...19T/abstract>`__
        is accepted, defaults to ``plr``
    :type lkl_name: str
    :param bin_method: Bin method used to split the color-magnitude diagram into cells
        (`Hess diagram <https://en.wikipedia.org/wiki/Hess_diagram>`__); one of:
        ``knuth, blocks, scott, freedman or fixed``.The method ``fixed`` uses (15, 10)
        bins in magnitude and color(s) respectively; defaults to ``knuth``
    :type bin_method: str

    :raises ValueError: If any of the attributes is not recognized as a valid option
    """

    def __init__(
        self, my_cluster: Cluster, lkl_name: str = "plr", bin_method: str = "knuth"
    ) -> None:
        self.my_cluster = my_cluster
        self.lkl_name = lkl_name
        self.bin_method = bin_method

        likelihoods = ("plr", "bins_distance", "chisq")
        if self.lkl_name not in likelihoods:
            raise ValueError(
                f"'{self.lkl_name}' not recognized. Should be one of {likelihoods}"
            )

        bin_methods = ("knuth", "blocks", "scott", "freedman", "fixed")
        if self.bin_method not in bin_methods:
            raise ValueError(
                f"Binning '{self.bin_method}' not recognized. "
                + f"Should be one of {bin_methods}"
            )

        # Obtain data used by the ``likelihood.get()`` method
        self.ranges, self.Nbins, self.cl_z_idx, self.cl_histo_f_z = lpriv.lkl_data(
            bin_method, my_cluster.mag_p, my_cluster.colors_p
        )

        self.max_lkl = 1
        if self.lkl_name == "plr":
            # Evaluate cluster against itself to obtain the maximum likelihood.
            # Since the initial max_lkl=1, subtracting 1 inverts it back to the
            # original likelihood value
            self.max_lkl = 1 - self.get(
                np.array([self.my_cluster.mag_p, *self.my_cluster.colors_p])
            )

        print("\nLikelihood object generated")

    def get(self, synth_clust: np.ndarray) -> float:
        """Evaluate the selected likelihood function.

        :param synth_clust:  Numpy array containing the synthetic cluster. The shape of
            this array must be: ``[magnitude, color1, (color2)]``, where ``magnitude``
            and ``color`` are arrays with the magnitude and color photometric data
            (``color2`` is the optional second color defined)
        :type synth_clust: np.ndarray

        :raise ValueError: If the likelihood function is not recognized

        :return: Likelihood value
        :rtype: float
        """
        if self.lkl_name == "plr":
            return lpriv.tremmel(
                self.ranges,
                self.Nbins,
                self.cl_z_idx,
                self.cl_histo_f_z,
                self.max_lkl,
                synth_clust,
            )
        # if self.lkl_name == "visual":
        #     return lpriv.visual(self, synth_clust)
        # if self.lkl_name == "mean_dist":
        #     return lpriv.mean_dist(self, synth_clust)
        elif self.lkl_name == "bins_distance":
            return lpriv.bins_distance(
                self.my_cluster.mag_p, self.my_cluster.colors_p, synth_clust
            )
        elif self.lkl_name == "chisq":
            return lpriv.chi_square(
                self.ranges, self.Nbins, self.cl_z_idx, self.cl_histo_f_z, synth_clust
            )
        else:
            raise ValueError(f"Likelihood '{self.lkl_name}' not recognized")

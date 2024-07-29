import numpy as np
from .cluster import Cluster
from .modules import cluster_priv as cp
from .modules.fastmp import fastMP
from .modules.bayesian_da import bayesian_mp


class Membership:
    """Define a :py:class:`Membership` object.

    This object is used as a container for membership probabilities methods. Currently
    two methods are included:

    :py:meth:`bayesian` : The algorithm was described in detail in
    the `article <https://doi.org/10.1051/0004-6361/201424946>`__
    were we originally introduced ``ASteCA``. The method requires ``(RA, DEC)``
    data and will use any extra data dimensions stored in the
    :py:class:`Cluster <asteca.cluster.Cluster>` object, i.e.: photometry,
    proper motions, and parallax. A minimum of two data dimensions are required.

    :py:meth:`fastmp` : The algorithm was described in detail in
    the `article <https://academic.oup.com/mnras/article/526/3/4107/7276628>`__
    were we introduced the `Unified Cluster Catalogue (UCC) <https://ucc.ar/>`__.
    The method requires proper motions, and parallax data dimensions stored in the
    :py:class:`Cluster <asteca.cluster.Cluster>` object. Photometry data is not
    employed.

    :param my_field: :py:class:`Cluster <asteca.cluster.Cluster>` object with the
        loaded data for the observed field
    :type my_field: Cluster
    :param seed: Random seed. If ``None`` a random integer will be generated and used,
        defaults to ``None``
    :type seed: int | None

    :raises ValueError: If there are missing required attributes in the
        :py:class:`Cluster <asteca.cluster.Cluster>` object
    """

    def __init__(
        self,
        my_field: Cluster,
        seed: int | None = None,
    ) -> None:
        self.my_field = my_field
        self.seed = seed

        noradeg_flag = False
        try:
            self.my_field.ra
        except AttributeError:
            noradeg_flag = True
        try:
            self.my_field.dec
        except AttributeError:
            noradeg_flag = True
        if noradeg_flag:
            raise ValueError(
                "The 'Membership' class requires (RA, DEC) data\n"
                + "to be present in the 'cluster' object"
            )
        if hasattr(self.my_field, "N_cluster") is False:
            raise ValueError(
                "The 'Membership' class requires the 'N_cluster' attribute\n"
                + "to be present in the 'cluster' object"
            )

        # Set seed
        self.rng = np.random.default_rng(seed)

        print(f"\nN_cluster      : {self.my_field.N_cluster}")
        print(f"Random seed    : {self.seed}")
        print("Membership object generated")

    def bayesian(self, N_runs: int = 1000, eq_to_gal: bool = False) -> np.ndarray:
        """Assign membership probabilities.

        Estimate the probability of being a true cluster member for all observed
        stars, using a Bayesian algorithm. The ``radec_c`` and ``radius`` attributes
        are required to be present in the
        :py:class:`Cluster <asteca.cluster.Cluster>` object.

        :param N_runs: Maximum number of runs, defaults to ``1000``
        :type N_runs: int
        :param eq_to_gal: Convert ``(RA, DEC)`` to ``(lon, lat)``. Useful for clusters
            with large ``DEC`` values to reduce the frame's distortion,
            defaults to ``False``
        :type eq_to_gal: bool

        :raises AttributeError: If either the ``radec_c`` or ``radius`` attributes
            are missing from the :py:class:`Cluster <asteca.cluster.Cluster>` object

        :return: Membership probabilities for all stars in the frame
        :rtype: np.ndarray
        """
        for attrib in ("radec_c", "radius"):
            if hasattr(self.my_field, attrib) is False:
                raise AttributeError(
                    f"Attribute '{attrib}' is required to be present "
                    + "in the 'cluster' object"
                )
        if N_runs < 10:
            raise AttributeError("Parameter 'N_runs' should be > 10")

        xc, yc = self.my_field.ra_v, self.my_field.dec_v
        center = self.my_field.radec_c
        cent_str = "radec_c "
        # Convert (RA, DEC) to (lon, lat)
        if eq_to_gal is True:
            xc, yc = cp.radec2lonlat(xc, yc)
            center = cp.radec2lonlat(*center)
            cent_str = "lonlat_c"

        print("\nRunning Bayesian DA...")
        print("{}       : ({:.4f}, {:.4f})".format(cent_str, *center))
        print(f"radius         : {self.my_field.radius} [deg]")
        print(f"N_cluster      : {self.my_field.N_cluster}")
        print(f"N_runs         : {N_runs}")

        # Generate input data array
        X = [xc, yc]
        e_X = []
        if hasattr(self.my_field, "mag_p"):
            X.append(self.my_field.mag_p)
            e_X.append(self.my_field.e_mag_p)
        if hasattr(self.my_field, "colors_p"):
            X.append(self.my_field.colors_p[0])
            e_X.append(self.my_field.e_colors_p[0])
            try:
                X.append(self.my_field.colors_p[1])
                e_X.append(self.my_field.e_colors_p[1])
            except IndexError:
                pass
        if hasattr(self.my_field, "plx_v"):
            X.append(self.my_field.plx_v)
            e_X.append(self.my_field.e_plx_v)
        if hasattr(self.my_field, "pmra_v"):
            X.append(self.my_field.pmra_v)
            e_X.append(self.my_field.e_pmra_v)
        if hasattr(self.my_field, "pmde_v"):
            X.append(self.my_field.pmde_v)
            e_X.append(self.my_field.e_pmde_v)
        X = np.array(X)
        e_X = np.array(e_X)

        probs = bayesian_mp(
            X, e_X, center, self.my_field.radius, self.my_field.N_cluster, N_runs
        )
        return probs

    def fastmp(
        self,
        fixed_centers: bool = False,
        N_runs: int = 1000,
        eq_to_gal: bool = True,
    ) -> np.ndarray:
        """Assign membership probabilities.

        Estimate the probability of being a true cluster member for all observed
        stars using the fastMP algorithm. The following data dimensions are required:
        ``(pmRA, pmDE, plx)``; photometry is not employed. Center estimates in
        ``(RA, DEC)``, as well as ``(pmRA, pmDE)`` and ``plx`` are required.

        :param fixed_centers: If ``True`` the center values (radec_c, pms_c, plx_c)
            stored in the :py:class:`Cluster <asteca.cluster.Cluster>` object will be
            kept fixed throughout the process, defaults to ``False``
        :type fixed_centers: bool
        :param N_runs: Maximum number of resamples, defaults to ``1000``
        :type N_runs: int
        :param eq_to_gal: Convert ``(RA, DEC)`` to ``(lon, lat)``. Useful for clusters
            with large ``DEC`` values to reduce the frame's distortion,
            defaults to ``True``
        :type eq_to_gal: bool

        :raises ValueError: If the :py:class:`Cluster <asteca.cluster.Cluster>` object
            is missing a required attribute:
            ``(ra, dec, pmra, pmde, plx, e_pmra, e_pmde, e_plx, radec_c, pms_c, plx_c)``

        :return: Membership probabilities for all stars in the frame
        :rtype: np.ndarray
        """
        for k in (
            "ra",
            "dec",
            "pmra",
            "pmde",
            "plx",
            "e_pmra",
            "e_pmde",
            "e_plx",
            "radec_c",
            "pms_c",
            "plx_c",
        ):
            if hasattr(self.my_field, k) is False:
                raise ValueError(f"'{k}' must be present as a 'cluster' attribute")
        if N_runs < 10:
            raise AttributeError("Parameter 'N_runs' should be > 10")

        xv, yv = self.my_field.ra_v, self.my_field.dec_v
        xy_center = self.my_field.radec_c
        cent_str = "radec_c "
        # Convert (RA, DEC) to (lon, lat)
        if eq_to_gal is True:
            xv, yv = cp.radec2lonlat(xv, yv)
            xy_center = cp.radec2lonlat(*xy_center)
            cent_str = "lonlat_c"

        print("\nRunning fastMP...")
        print("{}       : ({:.4f}, {:.4f})".format(cent_str, *xy_center))
        print("pms_c          : ({:.3f}, {:.3f})".format(*self.my_field.pms_c))
        print(f"plx_c          : {self.my_field.plx_c}")
        print(f"fixed_centers  : {fixed_centers}")
        print(f"N_cluster      : {self.my_field.N_cluster}")
        # centers_ex_flag = False
        # if centers_ex is not None:
        #     centers_ex_flag = True
        # print(f"centers_ex     : {centers_ex_flag}")
        print(f"N_resample     : {N_runs}")

        # Generate input data array for fastMP
        X = np.array(
            [
                xv,
                yv,
                self.my_field.pmra_v,
                self.my_field.pmde_v,
                self.my_field.plx_v,
                self.my_field.e_pmra_v,
                self.my_field.e_pmde_v,
                self.my_field.e_plx_v,
            ]
        )
        probs = fastMP(
            X,
            list(xy_center),
            list(self.my_field.pms_c),
            self.my_field.plx_c,
            fixed_centers,
            self.my_field.N_cluster,
            self.my_field.N_clust_min,
            self.my_field.N_clust_max,
            self.rng,
            N_runs,
        )

        return probs

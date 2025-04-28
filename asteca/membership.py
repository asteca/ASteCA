import numpy as np

from .cluster import Cluster


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
    :param verbose: Verbose level. A value of ``0`` hides all output, defaults to ``1``
    :type verbose: int

    :raises ValueError: If there are missing required attributes in the
        :py:class:`Cluster <asteca.cluster.Cluster>` object
    """

    def __init__(
        self,
        my_field: Cluster,
        seed: int | None = None,
        verbose: int = 1,
    ) -> None:
        self.my_field = my_field
        self.seed = seed
        self.verbose = verbose

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
        self.rng = np.random.default_rng(self.seed)

        self._vp("\nMembership object generated")
        self._vp(f"N_cluster      : {self.my_field.N_cluster}", 1)
        self._vp(f"Random seed    : {self.seed}", 1)

    def _vp(self, mssg: str, level: int = 0) -> None:
        """Verbose print method"""
        if self.verbose > level:
            print(mssg)

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
        :raises ValueError: If ``eq_to_gal`` is True and (RA, DEC) or center data are
            missing from the :py:class:`Cluster <asteca.cluster.Cluster>` object

        :return: Membership probabilities for all stars in the frame
        :rtype: np.ndarray
        """
        from .modules import cluster_priv as cp
        from .modules.bayesian_da import bayesian_mp

        txt = "Attribute '{}' is required to be present in the 'cluster' object"
        if hasattr(self.my_field, "radec_c") is False:
            raise AttributeError(txt.format("radec_c"))
        if hasattr(self.my_field, "radius") is False:
            raise AttributeError(txt.format("radius"))
        if N_runs < 10:
            raise AttributeError("Parameter 'N_runs' should be > 10")

        xc, yc = self.my_field.ra, self.my_field.dec
        center = self.my_field.radec_c
        cent_str = "radec_c "
        # Convert (RA, DEC) to (lon, lat)
        if eq_to_gal is True:
            if xc is None or yc is None or center is None:
                raise ValueError("(RA, DEC) data and center values are required")
            xc, yc = cp.radec2lonlat(xc, yc)
            center = cp.radec2lonlat(*center)
            cent_str = "lonlat_c"

        self._vp("\nRunning Bayesian DA")
        self._vp("{}       : ({:.4f}, {:.4f})".format(cent_str, *center), 1)
        self._vp(f"radius         : {self.my_field.radius:.4f} [deg]", 1)
        self._vp(f"N_cluster      : {self.my_field.N_cluster}", 1)
        self._vp(f"N_runs         : {N_runs}", 1)

        # Generate input data array
        X = [self.my_field.ra, self.my_field.dec]
        e_X = []
        if self.my_field.mag is not None:
            X.append(self.my_field.mag)
            e_X.append(self.my_field.e_mag)
        if self.my_field.color is not None:
            X.append(self.my_field.color)
            e_X.append(self.my_field.e_color)
        if self.my_field.color2 is not None:
            X.append(self.my_field.color2)
            e_X.append(self.my_field.e_color2)
        if self.my_field.plx is not None:
            X.append(self.my_field.plx)
            e_X.append(self.my_field.e_plx)
        if self.my_field.pmra is not None:
            X.append(self.my_field.pmra)
            e_X.append(self.my_field.e_pmra)
        if self.my_field.pmde is not None:
            X.append(self.my_field.pmde)
            e_X.append(self.my_field.e_pmde)
        X = np.array(X)
        e_X = np.array(e_X)

        out_mssg, probs = bayesian_mp(
            self.my_field.N_cluster,
            X,
            e_X,
            np.array(center),
            self.my_field.radius,
            N_runs,
            self.rng,
        )
        self._vp(out_mssg, 1)

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

        :raises AttributeError: If the :py:class:`Cluster <asteca.cluster.Cluster>`
            object is missing a required attribute:
            ``(ra, dec, pmra, pmde, plx, e_pmra, e_pmde, e_plx, radec_c, pms_c, plx_c)``
        :raises ValueError: If the ``N_runs`` parameter is less than 10

        :return: Membership probabilities for all stars in the frame
        :rtype: np.ndarray
        """
        from .modules import cluster_priv as cp
        from .modules.fastmp import fastMP

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
                raise AttributeError(f"'{k}' must be present as a 'cluster' attribute")
        if N_runs < 10:
            raise ValueError("Parameter 'N_runs' should be > 10")

        xv, yv = self.my_field.ra, self.my_field.dec
        xy_center = self.my_field.radec_c
        cent_str = "radec_c "
        # Convert (RA, DEC) to (lon, lat)
        if eq_to_gal is True:
            if xv is None or yv is None or xy_center is None:
                raise ValueError("(RA, DEC) data and center values are required")
            xv, yv = cp.radec2lonlat(xv, yv)
            xc, yc = cp.radec2lonlat(*xy_center)
            xy_center = (float(xc), float(yc))
            cent_str = "lonlat_c"

        self._vp("\nRunning fastMP")
        self._vp("{}       : ({:.4f}, {:.4f})".format(cent_str, *xy_center), 1)
        self._vp("pms_c          : ({:.3f}, {:.3f})".format(*self.my_field.pms_c), 1)
        self._vp(f"plx_c          : {self.my_field.plx_c:.3f}", 1)
        self._vp(f"fixed_centers  : {fixed_centers}", 1)
        self._vp(f"N_cluster      : {self.my_field.N_cluster}", 1)
        # centers_ex_flag = False
        # if centers_ex is not None:
        #     centers_ex_flag = True
        # self._vp(f"centers_ex     : {centers_ex_flag}", 1)
        self._vp(f"N_resample     : {N_runs}", 1)

        # Generate input data array for fastMP
        X = np.array(
            [
                xv,
                yv,
                self.my_field.pmra,
                self.my_field.pmde,
                self.my_field.plx,
                self.my_field.e_pmra,
                self.my_field.e_pmde,
                self.my_field.e_plx,
            ]
        )
        out_mssg, probs = fastMP(
            self.rng,
            xy_center,
            self.my_field.pms_c,
            self.my_field.plx_c,
            self.my_field.N_cluster,
            self.my_field.N_clust_min,
            self.my_field.N_clust_max,
            fixed_centers,
            X,
            N_runs,
        )

        self._vp(out_mssg, 1)

        return probs

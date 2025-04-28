import warnings

import numpy as np


class Cluster:
    """Define a :py:class:`Cluster` object.

    This object contains the basic data required to load a group of observed stars
    that could represent a cluster or an entire field.

    :param ra: Array that contains the right ascension (RA), defaults to ``None``
    :type ra: np.ndarray | None
    :param dec: Array that contains the declination (DEC), defaults to ``None``
    :type dec: np.ndarray | None
    :param magnitude: Array that contains the magnitude, defaults to ``None``
    :type magnitude: np.ndarray | None
    :param e_mag: Array that contains the magnitude's uncertainty, defaults to ``None``
    :type e_mag: np.ndarray | None
    :param color: Array that contains the color, defaults to ``None``
    :type color: np.ndarray | None
    :param e_color: Array that contains the color's uncertainty, defaults to ``None``
    :type e_color: np.ndarray | None
    :param color2: Array that contains the second color, defaults to ``None``
    :type color2: np.ndarray | None
    :param e_color2: Array that contains the second color's uncertainty, defaults to ``None``
    :type e_color2: np.ndarray | None
    :param plx: Array that contains the parallax, defaults to ``None``
    :type plx: np.ndarray | None
    :param e_plx: Array that contains the parallax uncertainty, defaults to ``None``
    :type e_plx: np.ndarray | None
    :param pmra: Array that contains the RA proper motion, defaults to ``None``
    :type pmra: np.ndarray | None
    :param e_pmra: Array that contains the RA proper motion's uncertainty, defaults to ``None``
    :type e_pmra: np.ndarray | None
    :param pmde: Array that contains the DEC proper motion, defaults to ``None``
    :type pmde: np.ndarray | None
    :param e_pmde: Array that contains the DEC proper motion's uncertainty, defaults to ``None``
    :type e_pmde: np.ndarray | None
    :param N_clust_min: Minimum number of cluster members, defaults to ``25``
    :type N_clust_min: int
    :param N_clust_max: Maximum number of cluster members, defaults to ``5000``
    :type N_clust_max: int
    :param verbose: Verbose level. A value of ``0`` hides all output, defaults to ``1``
    :type verbose: int

    """

    def __init__(
        self,
        ra: np.ndarray | None = None,
        dec: np.ndarray | None = None,
        magnitude: np.ndarray | None = None,
        e_mag: np.ndarray | None = None,
        color: np.ndarray | None = None,
        e_color: np.ndarray | None = None,
        color2: np.ndarray | None = None,
        e_color2: np.ndarray | None = None,
        plx: np.ndarray | None = None,
        e_plx: np.ndarray | None = None,
        pmra: np.ndarray | None = None,
        e_pmra: np.ndarray | None = None,
        pmde: np.ndarray | None = None,
        e_pmde: np.ndarray | None = None,
        N_clust_min: int = 25,
        N_clust_max: int = 5000,
        verbose: int = 1,
    ) -> None:
        self.ra = ra
        self.dec = dec
        self.magnitude = magnitude
        self.e_mag = e_mag
        self.color = color
        self.e_color = e_color
        self.color2 = color2
        self.e_color2 = e_color2
        self.plx = plx
        self.pmra = pmra
        self.pmde = pmde
        self.e_plx = e_plx
        self.e_pmra = e_pmra
        self.e_pmde = e_pmde
        self.N_clust_min = N_clust_min
        self.N_clust_max = N_clust_max
        self.verbose = verbose

        self._vp("\nInstantiating cluster")
        self._load_column_data()
        self._vp(f"N_stars        : {self.N_stars}", 1)
        self._vp(f"N_clust_min    : {self.N_clust_min}", 1)
        self._vp(f"N_clust_max    : {self.N_clust_max}", 1)
        self._vp("Cluster object generated")

    def _vp(self, mssg: str, level: int = 0) -> None:
        """Verbose print method"""
        if self.verbose > level:
            print(mssg)

    def _load_column_data(self):
        cols_read = []
        N_stars_lst = []

        if self.ra is not None:
            self.ra = np.asarray(self.ra, dtype=float)
            cols_read.append("RA")
            N_stars_lst.append(len(self.ra))

        if self.dec is not None:
            self.dec = np.asarray(self.dec, dtype=float)
            cols_read.append("DEC")
            N_stars_lst.append(len(self.dec))

        if self.magnitude is not None:
            if self.e_mag is None:
                raise ValueError("Magnitude uncertainty is required")
            self.mag = np.asarray(self.magnitude, dtype=float)
            self.e_mag = np.asarray(self.e_mag, dtype=float)
            cols_read.append("Magnitude")
            cols_read.append("e_mag")
            N_stars_lst.append(len(self.mag))
            N_stars_lst.append(len(self.e_mag))

        if self.color is not None:
            if self.e_color is None:
                raise ValueError("Color uncertainty is required")
            self.color = np.asarray(self.color, dtype=float)
            self.e_color = np.asarray(self.e_color, dtype=float)
            cols_read.append("Color")
            cols_read.append("e_color")
            N_stars_lst.append(len(self.color))
            N_stars_lst.append(len(self.e_color))

        if self.color is None and self.color2 is not None:
            raise ValueError("Argument 'color' must be defined if 'color2' is defined")

        if self.color2 is not None:
            if self.e_color2 is None:
                raise ValueError("Color2 uncertainty is required")
            self.color2 = np.asarray(self.color2, dtype=float)
            self.e_color2 = np.asarray(self.e_color2, dtype=float)
            cols_read.append("Color2")
            cols_read.append("e_color2")
            N_stars_lst.append(len(self.color2))
            N_stars_lst.append(len(self.e_color2))

        if self.plx is not None:
            if self.e_plx is None:
                raise ValueError("Parallax uncertainty is required")
            self.plx = np.asarray(self.plx, dtype=float)
            self.e_plx = np.asarray(self.e_plx, dtype=float)
            cols_read.append("Plx")
            N_stars_lst.append(len(self.plx))
            N_stars_lst.append(len(self.e_plx))

        if self.pmra is not None:
            if self.e_pmra is None:
                raise ValueError("pmRA uncertainty is required")
            self.pmra = np.asarray(self.pmra, dtype=float)
            self.e_pmra = np.asarray(self.e_pmra, dtype=float)
            cols_read.append("pmRA")
            N_stars_lst.append(len(self.pmra))
            N_stars_lst.append(len(self.e_pmra))

        if self.pmde is not None:
            if self.e_pmde is None:
                raise ValueError("pmDE uncertainty is required")
            self.pmde = np.asarray(self.pmde, dtype=float)
            self.e_pmde = np.asarray(self.e_pmde, dtype=float)
            cols_read.append("pmDE")
            N_stars_lst.append(len(self.pmde))
            N_stars_lst.append(len(self.e_pmde))

        if len(cols_read) == 0:
            raise ValueError("No column names defined for cluster")
        self._vp(f"Columns read   : {', '.join(cols_read)}", 1)
        if not np.all(N_stars_lst):
            raise ValueError("Data arrays have different lengths")
        self.N_stars = N_stars_lst[0]
        if self.N_stars == 0:
            raise ValueError("Arrays are empty")

    def get_center(
        self,
        algo: str = "knn_5d",
        data_2d: str = "radec",
        radec_c: tuple[float, float] | None = None,
        pms_c: tuple[float, float] | None = None,
        plx_c: float | None = None,
    ) -> None:
        """Estimate center coordinates for the cluster

        Use the available data (ra, dec, pmra, pmde, plx) to estimate a cluster's
        center coordinates as the point(s) of maximum density. Algorithms:

        - ``knn_5d``: Estimates the 5-dimensional center values (in
          ``ra, dec, pmra, pmde, plx``) as the median position of the k
          (``k=N_clust_min``) stars with the largest nearest-neighbor density to a
          5D center value that can be either given (fully or partially) or estimated.
        - ``kde_2d``: Estimates the 2-dimensional center values (either in ``(ra, dec)``
          or in ``(pmra, pmde)``; determined by the ``data_2d`` argument) using a
          Kernel Density Estimator (KDE).

        :param algo: Algorithm used to estimate center values, one of
            (``knn_5d``, ``kde_2d``), defaults to ``knn_5d``
        :type algo: str
        :param data_2d: String indicating the data to be used to estimate
            the center value, either: ``radec`` or ``pms``, defaults to ``radec``
        :type data_2d: str
        :param radec_c: Estimated value for the (RA, DEC) center, defaults to ``None``
        :type radec_c: tuple[float, float] | None
        :param pms_c: Estimated value for the (pmRA, pmDE) center, defaults to ``None``
        :type pms_c: tuple[float, float] | None
        :param plx_c: Estimated value for the plx center, defaults to ``None``
        :type plx_c: float | None

        :raises ValueError: If required data is missing from the
            :py:class:`Cluster <asteca.cluster.Cluster>` object
        """
        from .modules import cluster_priv as cp

        if algo == "knn_5d":
            if (
                self.ra is None
                or self.dec is None
                or self.pmra is None
                or self.pmde is None
                or self.plx is None
            ):
                raise ValueError(
                    "Algorithm 'knn_5d' requires (ra, dec, pmra, pmde, plx) data"
                )

            # To galactic coordinates (not optional, always use instead of equatorial)
            glon, glat = cp.radec2lonlat(self.ra, self.dec)
            lonlat_c = None
            if radec_c is not None:
                lon, lat = cp.radec2lonlat(radec_c[0], radec_c[1])
                lonlat_c = (lon, lat)

            # Remove nans
            X = np.array([glon, glat, self.pmra, self.pmde, self.plx])
            # Reject nan values and extract clean data
            _, X_no_nan = cp.reject_nans(X)
            lon, lat, pmRA, pmDE, plx = X_no_nan

            x_c, y_c, pmra_c, pmde_c, plx_c = cp.get_knn_5D_center(
                lon,
                lat,
                pmRA,
                pmDE,
                plx,
                xy_c=lonlat_c,
                vpd_c=pms_c,
                plx_c=plx_c,
                N_clust_min=self.N_clust_min,
                N_clust_max=self.N_clust_max,
            )
            ra_c, dec_c = cp.lonlat2radec(x_c, y_c)

            self.radec_c = (ra_c, dec_c)
            self.pms_c = (pmra_c, pmde_c)
            self.plx_c = plx_c

            mssg = f"radec_c        : ({ra_c:.4f}, {dec_c:.4f})\n"
            mssg += f"pms_c          : ({pmra_c:.3f}, {pmde_c:.3f})\n"
            mssg += f"plx_c          : {plx_c:.3f}"

        elif algo == "kde_2d":
            if data_2d == "radec":
                if any([_ is None for _ in (self.ra, self.dec)]):
                    raise ValueError("Data for (ra, dec) is required")
                c_str = "radec_c"
                x, y = self.ra, self.dec
            elif data_2d == "pms":
                if any([_ is None for _ in (self.pmra, self.pmde)]):
                    raise ValueError("Data for (pmra, pmde) is required")
                c_str = "pms_c"
                x, y = self.pmra, self.pmde
            else:
                raise ValueError(f"Algorithm '{data_2d}' argument not recognized")

            # Remove rows containing any NaNs
            array_2d = np.array([x, y])
            x, y = array_2d[:, ~np.isnan(array_2d).any(axis=0)]

            x_c, y_c = cp.get_2D_center(x, y)

            if c_str == "radec_c":
                self.radec_c = (x_c, y_c)
            if c_str == "pms_c":
                self.pms_c = (x_c, y_c)

            mssg = "{}        : ({:.4f}, {:.4f})".format(c_str, x_c, y_c)

        else:
            raise ValueError(f"Selected method '{algo}' not recognized")

        self._vp("\nCenter coordinates found")
        self._vp(mssg, 1)

    def _get_radius(self, algo: str = "king") -> None:
        """Estimate the cluster radius

        - ``field_dens``: xxx
        - ``king``: xxx

        :param algo: Algorithm used to estimate the radius, one of
            (``field_dens, king``); defaults to ``field_dens``
        :type algo: str

        :raises ValueError: If required attributes are  missing from the
            :py:class:`Cluster <asteca.cluster.Cluster>` object
        """

        if algo == "field_dens":
            for k in ("ra", "dec", "radec_c"):
                if hasattr(self, k) is False:
                    raise ValueError(f"'{k}' must be present as a 'cluster' attribute")
        elif algo == "king":
            for k in ("ra", "dec", "radec_c"):
                if hasattr(self, k) is False:
                    raise ValueError(f"'{k}' must be present as a 'cluster' attribute")
        else:
            raise ValueError(f"Unrecognized method '{algo}")

        from .modules import cluster_priv as cp

        xv, yv = self.ra, self.dec
        xy_center = self.radec_c

        if algo == "field_dens":
            radius = cp.fdens_radius(xv, yv, list(xy_center))
            self._vp("\nRadius estimated")
            self._vp(f"Radius         : {radius:.3f}", 1)
            self.radius = radius
        elif algo == "king":
            amplitude, radius_core, radius, field_dens = cp.king_radius(
                xv, yv, xy_center
            )
            self._vp("\nKing parameters density estimated")
            self._vp(f"Field density  : {field_dens:.3f}", 1)
            self._vp(f"Amplitude      : {amplitude:.3f}", 1)
            self._vp(f"Radius (core)  : {radius_core:.3f}", 1)
            self._vp(f"Radius (tidal) : {radius:.3f}", 1)
            self.field_density = field_dens
            self.amplitude = amplitude
            self.radius_core = radius_core
            self.radius = radius

    def get_nmembers(self, algo: str = "ripley", eq_to_gal: bool = True) -> None:
        """Estimate the number of members for the cluster. Algorithms:

        - ``ripley``: Originally introduced with the ``fastMP`` membership method
          in `Perren et al.
          (2023) <https://academic.oup.com/mnras/article/526/3/4107/7276628>`__.
          Requires ``(ra, dec, pmra, pmde, plx)`` and their center estimates.
        - ``density``: Simple algorithm that counts the number of stars within the
          cluster region (center+radius) and subtracts the expected number of field
          stars within that region. Requires the ``(ra, dec)`` center and the radius of
          the cluster to be defined.

        :param algo: Algorithm used to estimate center values, one of
            (``ripley, density``); defaults to ``ripley``
        :type algo: str
        :param eq_to_gal: Convert ``(ra, dec)`` to ``(lon, lat)``. Useful for clusters
            with large ``dec`` values to reduce the frame's distortion,
            defaults to ``False``
        :type eq_to_gal: bool

        :raises ValueError: If ``algo`` argument is not recognized
        :raises AttributeError: If required attributes are  missing from the
            :py:class:`Cluster <asteca.cluster.Cluster>` object
        """
        from .modules import cluster_priv as cp
        from .modules import nmembers as nm

        if algo not in ("ripley", "density"):
            raise ValueError(f"Selected method '{algo}' not recognized")

        # Both methods require these
        if self.ra is None or self.dec is None:
            raise AttributeError(
                "The arguments (ra, dec) must be present as a 'cluster' attribute"
            )

        # (x, y) coordinates
        xv, yv = self.ra, self.dec
        # Convert (RA, DEC) to (lon, lat)
        if eq_to_gal is True:
            xv, yv = cp.radec2lonlat(xv, yv)

        if algo == "ripley":
            if (
                self.pmra is None
                or self.pmde is None
                or self.plx is None
                or hasattr(self, "pms_c") is False
                or hasattr(self, "plx_c") is False
            ):
                raise AttributeError(
                    "The arguments (pmra, pmdec, plx, pms_c, plx_c) must be present as a 'cluster' attribute"
                )
            N_cluster = nm.ripley_nmembs(
                xv,
                yv,
                self.pmra,
                self.pmde,
                self.plx,
                self.pms_c,
                self.plx_c,
            )
        elif algo == "density":
            if hasattr(self, "radec_c") is False or hasattr(self, "radius") is False:
                raise AttributeError(
                    "The arguments (radec_c, radius) must be present as a 'cluster' attribute"
                )

            ra_c, dec_c = self.radec_c
            xy_center = (ra_c, dec_c)
            if eq_to_gal is True:
                lon, lat = cp.radec2lonlat(*xy_center)
                xy_center = (lon, lat)

            N_cluster = nm.density_nmembs(xv, yv, xy_center, self.radius)

        if N_cluster < self.N_clust_min:
            warnings.warn(
                "The estimated number of cluster members is " + f"<{self.N_clust_min}"
            )
            N_cluster = self.N_clust_min
        if N_cluster > self.N_clust_max:
            warnings.warn(
                "The estimated number of cluster members is " + f">{self.N_clust_max}"
            )
            N_cluster = self.N_clust_max

        self.N_cluster = N_cluster
        self._vp("\nNumber of members estimated")
        self._vp(f"N_cluster      : {N_cluster}", 1)

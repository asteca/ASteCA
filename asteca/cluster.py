import warnings

import numpy as np


class Cluster:
    """Define a :py:class:`Cluster` object.

    This object contains the basic data required to load a group of observed stars
    that could represent a cluster or an entire field.

    :param cluster_name: Name of the cluster, only required if the data is to be
        loaded from the UCC members file. If provided, the ``UCC_file_path`` argument
        must also be provided
    :type cluster_name: str | None
    :param UCC_file_path: Path to the UCC members file containing the cluster data. If
        ``cluster_name`` is provided this argument must also be provided
    :type UCC_file_path: str | None
    :param ra: Array that contains the right ascension (RA)
    :type ra: np.ndarray | None
    :param dec: Array that contains the declination (DEC)
    :type dec: np.ndarray | None
    :param mag: Array that contains the magnitude
    :type mag: np.ndarray | None
    :param max_mag: Maximum magnitude to be considered for the cluster
    :type max_mag: float | None
    :param e_mag: Array that contains the magnitude's uncertainty
    :type e_mag: np.ndarray | None
    :param color: Array that contains the color
    :type color: np.ndarray | None
    :param e_color: Array that contains the color's uncertainty
    :type e_color: np.ndarray | None
    :param color2: Array that contains the second color
    :type color2: np.ndarray | None
    :param e_color2: Array that contains the second color's uncertainty
    :type e_color2: np.ndarray | None
    :param plx: Array that contains the parallax
    :type plx: np.ndarray | None
    :param e_plx: Array that contains the parallax uncertainty
    :type e_plx: np.ndarray | None
    :param pmra: Array that contains the RA proper motion
    :type pmra: np.ndarray | None
    :param e_pmra: Array that contains the RA proper motion's uncertainty
    :type e_pmra: np.ndarray | None
    :param pmde: Array that contains the DEC proper motion
    :type pmde: np.ndarray | None
    :param e_pmde: Array that contains the DEC proper motion's uncertainty
    :type e_pmde: np.ndarray | None
    :param N_clust_min: Minimum number of cluster members
    :type N_clust_min: int
    :param N_clust_max: Maximum number of cluster members
    :type N_clust_max: int
    :param verbose: Verbose level. A value of ``0`` hides all output
    :type verbose: int

    :raises ValueError: If ``cluster_name`` is provided but ``UCC_file_path`` is not provided

    """

    def __init__(
        self,
        cluster_name: str | None = None,
        UCC_file_path: str | None = None,
        ra: np.ndarray | None = None,
        dec: np.ndarray | None = None,
        mag: np.ndarray | None = None,
        max_mag: float | None = None,
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
        N_clust_max: int = 2000,
        verbose: int = 1,
    ) -> None:
        self.cluster_name = cluster_name
        self.UCC_file_path = UCC_file_path
        self.ra = ra
        self.dec = dec
        self.mag = mag
        self.max_mag = max_mag
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

        if self.cluster_name is not None and self.UCC_file_path is None:
            raise ValueError(
                "If a cluster name is provided, the 'UCC_file_path' must also be provided"
            )

        if self.cluster_name is not None and self.UCC_file_path is not None:
            self._vp(f"\nInstantiating cluster '{self.cluster_name}'")
            self._vp(f"Loading data from '{self.UCC_file_path}'")
            self._load_UCC_data()
        else:
            self._vp("\nInstantiating cluster from data")
            self._load_column_data()

        self._vp(f"N_stars        : {self.N_stars}", 1)
        self._vp(f"N_clust_min    : {self.N_clust_min}", 1)
        self._vp(f"N_clust_max    : {self.N_clust_max}", 1)
        self._vp("Cluster object generated")

    def _vp(self, mssg: str, level: int = 0) -> None:
        """Verbose print method"""
        if self.verbose > level:
            print(mssg)

    def _load_UCC_data(self):

        if not isinstance(self.cluster_name, str):
            raise TypeError("'cluster_name' must be a string")

        if not isinstance(self.UCC_file_path, str):
            raise TypeError("'UCC_file_path' must be a string")

        import pandas as pd

        df = pd.read_parquet(self.UCC_file_path)

        cluster_fname = (
            self.cluster_name.lower()
            .replace(" ", "")
            .replace("_", "")
            .replace("-", "")
            .replace(".", "")
            .replace("+", "p")
        )

        if cluster_fname not in df["name"].values:
            raise ValueError(
                f"Cluster name '{self.cluster_name}' not found in UCC file"
            )

        cluster_df = df[df["name"] == cluster_fname]
        if len(cluster_df) == 0:
            raise ValueError(
                f"No data found for cluster name '{self.cluster_name}' in UCC file"
            )

        if self.max_mag is not None:
            cluster_df = cluster_df[cluster_df["Gmag"] <= self.max_mag]
            if len(cluster_df) == 0:
                raise ValueError(
                    f"No stars with magnitude <= {self.max_mag} found for cluster '{self.cluster_name}'"
                )

        self.ra = np.asarray(cluster_df["RA_ICRS"], dtype=float)
        self.dec = np.asarray(cluster_df["DE_ICRS"], dtype=float)
        self.mag = np.asarray(cluster_df["Gmag"], dtype=float)
        self.e_mag = np.asarray(cluster_df["e_Gmag"], dtype=float)
        self.color = np.asarray(cluster_df["BP-RP"], dtype=float)
        self.e_color = np.asarray(cluster_df["e_BP-RP"], dtype=float)
        self.plx = np.asarray(cluster_df["Plx"], dtype=float)
        self.e_plx = np.asarray(cluster_df["e_Plx"], dtype=float)
        self.pmra = np.asarray(cluster_df["pmRA"], dtype=float)
        self.e_pmra = np.asarray(cluster_df["e_pmRA"], dtype=float)
        self.pmde = np.asarray(cluster_df["pmDE"], dtype=float)
        self.e_pmde = np.asarray(cluster_df["e_pmDE"], dtype=float)
        cols_read = (
            "RA",
            "DEC",
            "Gmag",
            "e_Gmag",
            "BP-RP",
            "e_BP-RP",
            "Plx",
            "e_Plx",
            "pmRA",
            "e_pmRA",
            "pmDE",
            "e_pmDE",
        )
        self._vp(f"Columns read   : {', '.join(cols_read)}", 1)
        self.N_stars = len(cluster_df)

    def _process_array(
        self, name, err_name=None, positive_err=False, cols_read=None, N_stars_lst=None
    ):
        """
        Convert attribute to float array and validate its uncertainty (if any).
        """
        value = getattr(self, name)
        if value is None:
            return

        arr = np.asarray(value, dtype=float)
        setattr(self, name, arr)

        if cols_read is not None:
            cols_read.append(name)
        if N_stars_lst is not None:
            N_stars_lst.append(len(arr))

        if err_name is not None:
            err = getattr(self, err_name)
            if err is None:
                raise ValueError(f"{name} uncertainty is required")

            err_arr = np.asarray(err, dtype=float)
            if positive_err and (err_arr <= 0).any():
                raise ValueError(f"{name} uncertainties must be > 0")

            setattr(self, err_name, err_arr)

            if cols_read is not None:
                cols_read.append(err_name)
            if N_stars_lst is not None:
                N_stars_lst.append(len(err_arr))

    def _load_column_data(self):
        cols_read = []
        N_stars_lst = []

        # Single arrays
        self._process_array("ra", cols_read=cols_read, N_stars_lst=N_stars_lst)
        self._process_array("dec", cols_read=cols_read, N_stars_lst=N_stars_lst)

        # Value + uncertainty pairs
        pairs = [
            ("mag", "e_mag"),
            ("color", "e_color"),
            ("color2", "e_color2"),
            ("plx", "e_plx"),
            ("pmra", "e_pmra"),
            ("pmde", "e_pmde"),
        ]

        for name, err in pairs:
            self._process_array(
                name,
                err_name=err,
                positive_err=True,
                cols_read=cols_read,
                N_stars_lst=N_stars_lst,
            )

        # Logical dependency
        if self.color is None and self.color2 is not None:
            raise ValueError("Argument 'color' must be defined if 'color2' is defined")

        if not cols_read:
            raise ValueError("No column names defined for cluster")

        self._vp(f"Columns read   : {', '.join(cols_read)}", 1)

        if len(set(N_stars_lst)) > 1:
            raise ValueError("Data arrays have different lengths")

        # Apply magnitude cut generically
        if self.max_mag is not None:
            if self.mag is None:
                raise ValueError("Magnitude data is required to apply 'max_mag' filter")

            mask = self.mag <= self.max_mag
            if not mask.any():
                raise ValueError(
                    f"No stars with magnitude <= {self.max_mag} found in provided data"
                )

            for name in [
                "ra",
                "dec",
                "mag",
                "e_mag",
                "color",
                "e_color",
                "color2",
                "e_color2",
                "plx",
                "e_plx",
                "pmra",
                "e_pmra",
                "pmde",
                "e_pmde",
            ]:
                arr = getattr(self, name)
                if arr is not None:
                    setattr(self, name, arr[mask])

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
            (``knn_5d``, ``kde_2d``)
        :type algo: str
        :param data_2d: String indicating the data to be used to estimate
            the center value, either: ``radec`` or ``pms``
        :type data_2d: str
        :param radec_c: Estimated value for the (RA, DEC) center
        :type radec_c: tuple[float, float] | None
        :param pms_c: Estimated value for the (pmRA, pmDE) center
        :type pms_c: tuple[float, float] | None
        :param plx_c: Estimated value for the plx center
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

            x_c, y_c, pmra_c, pmde_c, plx_c = cp.get_knn_5D_center(
                glon,
                glat,
                self.pmra,
                self.pmde,
                self.plx,
                xy_c=lonlat_c,
                vpd_c=pms_c,
                plx_c=plx_c,
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

        :param algo: Algorithm used to estimate the radius, one of (``field_dens, king``)
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

        :param algo: Algorithm used to estimate center values, one of (``ripley, density``)
        :type algo: str
        :param eq_to_gal: Convert ``(ra, dec)`` to ``(lon, lat)``. Useful for clusters
            with large ``dec`` values to reduce the frame's distortion
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

        if hasattr(self, "radec_c") is False:
            raise AttributeError(
                "The argument 'radec_c' must be present as a 'cluster' attribute"
            )

        # (x, y) coordinates
        xv, yv = self.ra, self.dec
        xy_center = (self.radec_c[0], self.radec_c[1])
        # Convert (RA, DEC) to (lon, lat)
        if eq_to_gal is True:
            # Add center as last element
            xv_t = np.array(list(xv) + [xy_center[0]])
            yv_t = np.array(list(yv) + [xy_center[1]])
            # Convert
            xv_t, yv_t = cp.radec2lonlat(xv_t, yv_t)
            # Extract center value
            xv, yv = xv_t[:-1], yv_t[:-1]
            xy_center = (xv_t[-1], yv_t[-1])

        if algo == "ripley":
            if (
                self.pmra is None
                or self.pmde is None
                or self.plx is None
                or hasattr(self, "pms_c") is False
                or hasattr(self, "plx_c") is False
            ):
                raise AttributeError(
                    "The arguments (pmra, pmde, plx, pms_c, plx_c) must be present as a 'cluster' attribute"
                )
            idx_selected = nm.ripley_nmembs(
                xv,
                yv,
                self.pmra,
                self.pmde,
                self.plx,
                xy_center,
                self.pms_c,
                self.plx_c,
                self.N_clust_min,
                self.N_clust_max,
            )
            # Stored for eventual use
            self.ripley_idx_selected = idx_selected
            # This is the value expected from this method
            N_cluster = len(idx_selected)

        elif algo == "density":
            if hasattr(self, "radius") is False:
                raise AttributeError(
                    "The argument 'radius' must be present as a 'cluster' attribute"
                )
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

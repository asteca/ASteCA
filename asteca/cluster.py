import warnings
import numpy as np
import pandas as pd
from .modules import cluster_priv as cp
from .modules import nmembers as nm


class Cluster:
    """Define a :py:class:`Cluster` object.

    This object contains the basic data required to load a group of observed stars
    that could represent a cluster or an entire field.

    :param obs_df: `pandas DataFrame <https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html>`__
        with the observed loaded data.
    :type obs_df: pd.DataFrame
    :param ra: Name of the DataFrame column that contains the right ascension (RA),
        defaults to ``None``
    :type ra: str | None
    :param dec: Name of the DataFrame column that contains the declination (DEC),
        defaults to ``None``
    :type dec: str | None
    :param magnitude: Name of the DataFrame column that contains the magnitude,
        defaults to ``None``
    :type magnitude: str | None
    :param e_mag: Name of the DataFrame column that contains the magnitude's
        uncertainty, defaults to ``None``
    :type e_mag: str | None
    :param color: Name of the DataFrame column that contains the color, defaults to
        ``None``
    :type color: str | None
    :param e_color: Name of the DataFrame column that contains the color's uncertainty,
        defaults to ``None``
    :type e_color: str | None
    :param color2: Name of the DataFrame column that contains the second color,
        defaults to ``None``
    :type color2: str | None
    :param e_color2: Name of the DataFrame column that contains the second color's
        uncertainty, defaults to ``None``
    :type e_color2: str | None
    :param plx: Name of the DataFrame column that contains the parallax,
        defaults to ``None``
    :type plx: str | None
    :param e_plx: Name of the DataFrame column that contains the parallax uncertainty,
        defaults to ``None``
    :type e_plx: str | None
    :param pmra: Name of the DataFrame column that contains the RA proper motion,
        defaults to ``None``
    :type pmra: str | None
    :param e_pmra: Name of the DataFrame column that contains the RA proper motion's
        uncertainty, defaults to ``None``
    :type e_pmra: str | None
    :param pmde: Name of the DataFrame column that contains the DEC proper motion,
        defaults to ``None``
    :type pmde: str | None
    :param e_pmde: Name of the DataFrame column that contains the DEC proper motion's
        uncertainty, defaults to ``None``
    :type e_pmde: str | None
    :param N_clust_min: Minimum number of cluster members, defaults to ``25``
    :type N_clust_min: int
    :param N_clust_max: Maximum number of cluster members, defaults to ``5000``
    :type N_clust_max: int

    :raises ValueError: If there are missing required attributes to generate the
        :py:class:`Cluster <asteca.cluster.Cluster>` object
    """

    def __init__(
        self,
        obs_df: pd.DataFrame,
        ra: str | None = None,
        dec: str | None = None,
        magnitude: str | None = None,
        e_mag: str | None = None,
        color: str | None = None,
        e_color: str | None = None,
        color2: str | None = None,
        e_color2: str | None = None,
        plx: str | None = None,
        e_plx: str | None = None,
        pmra: str | None = None,
        e_pmra: str | None = None,
        pmde: str | None = None,
        e_pmde: str | None = None,
        N_clust_min: int = 25,
        N_clust_max: int = 5000,
    ) -> None:
        self.obs_df = obs_df
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

        self.N_stars = len(self.obs_df)
        print("\nInstantiating cluster...")
        print(f"N_stars        : {self.N_stars}")
        print(f"N_clust_min    : {self.N_clust_min}")
        print(f"N_clust_max    : {self.N_clust_max}")

        dim_count = self._load()
        if dim_count > 0:
            print("Cluster object generated")
        else:
            raise ValueError("No column names defined for cluster")

    def _load(self):
        dim_count = 0

        if self.ra is not None:
            self.ra_v = np.array(self.obs_df[self.ra], dtype=float)
            self.dec_v = np.array(self.obs_df[self.dec], dtype=float)
            print(f"(RA, DEC)      : ({self.ra}, {self.dec})")
            dim_count += 1

        # TODO: change _p to _v
        if self.magnitude is not None:
            self.mag_p = np.array(self.obs_df[self.magnitude], dtype=float)
            if self.e_mag is not None:
                self.e_mag_p = np.array(self.obs_df[self.e_mag], dtype=float)
            else:
                raise ValueError("Magnitude uncertainty is required")
            print(f"Magnitude      : {self.magnitude} [{self.e_mag}]")
            dim_count += 1

        if self.color is not None:
            self.colors_p = [np.array(self.obs_df[self.color], dtype=float)]
            if self.e_color is not None:
                self.e_colors_p = [np.array(self.obs_df[self.e_color], dtype=float)]
            else:
                raise ValueError("Color uncertainty is required")
            print(f"Color          : {self.color} [{self.e_color}]")
            dim_count += 1
            if self.color2 is not None:
                self.colors_p.append(np.array(self.obs_df[self.color2], dtype=float))
                if self.e_color2 is not None:
                    self.e_colors_p.append(
                        np.array(self.obs_df[self.e_color2], dtype=float)
                    )
                else:
                    raise ValueError("Color2 uncertainty is required")
                print(f"Color2         : {self.color2} [{self.e_color2}]")
                dim_count += 1

        if self.plx is not None:
            self.plx_v = np.array(self.obs_df[self.plx], dtype=float)
            if self.e_plx is not None:
                self.e_plx_v = np.array(self.obs_df[self.e_plx], dtype=float)
            else:
                raise ValueError("Parallax uncertainty is required")
            print(f"plx            : {self.plx} [{self.e_plx}]")
            dim_count += 1

        if self.pmra is not None:
            self.pmra_v = np.array(self.obs_df[self.pmra], dtype=float)
            if self.e_pmra is not None:
                self.e_pmra_v = np.array(self.obs_df[self.e_pmra], dtype=float)
            else:
                raise ValueError("pmRA uncertainty is required")
            print(f"pmRA           : {self.pmra} [{self.e_pmra}]")
            dim_count += 1
        if self.pmde is not None:
            self.pmde_v = np.array(self.obs_df[self.pmde], dtype=float)
            if self.e_pmra is not None:
                self.e_pmde_v = np.array(self.obs_df[self.e_pmde], dtype=float)
            else:
                raise ValueError("pmDE uncertainty is required")
            print(f"pmDE           : {self.pmra} [{self.e_pmde}]")
            dim_count += 1

        return dim_count

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
        center coordinates as the point of maximum density. Methods:

        ``knn_5d``: Estimates the center value as the median position of the k
        (k=N_clust_min) nearest stars to an estimate of the center in proper motions
        and (ra, dec, plx), if given.

        ``kde_2d``: Estimates the center value using a Kernel Density Estimator (KDE)
        given a two dimensional array determined by the ``data_2d`` argument.

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

        if algo == "knn_5d":
            if any(
                [_ is None for _ in (self.ra, self.dec, self.pmra, self.pmde, self.plx)]
            ):
                raise ValueError(
                    f"Algorithm '{algo}' requires (ra, dec, pmra, pmde, plx) data "
                    + "to be defined"
                )

            # To galactic coordinates (not optional, always use instead of equatorial)
            glon, glat = cp.radec2lonlat(self.ra_v, self.dec_v)
            lonlat_c = None
            if radec_c is not None:
                lonlat_c = cp.radec2lonlat(radec_c[0], radec_c[1])

            # Remove nans
            X = np.array([glon, glat, self.pmra_v, self.pmde_v, self.plx_v])
            # Reject nan values and extract clean data
            _, X_no_nan = cp.reject_nans(X)
            lon, lat, pmRA, pmDE, plx = X_no_nan

            x_c, y_c, pmra_c, pmde_c, plx_c = cp.get_5D_center(
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

            print("\nCenter coordinates found:")
            print("radec_c        : ({:.4f}, {:.4f})".format(ra_c, dec_c))
            print("pms_c          : ({:.3f}, {:.3f})".format(pmra_c, pmde_c))
            print("plx_c          : {:.3f}".format(plx_c))

            self.radec_c = [ra_c, dec_c]
            self.pms_c = [pmra_c, pmde_c]
            self.plx_c = plx_c

        elif algo == "kde_2d":
            if data_2d == "radec":
                if any([_ is None for _ in (self.ra, self.dec)]):
                    raise ValueError("Data for  (ra, dec) data is required")
                c_str = "radec_c"
                x, y = self.ra_v, self.dec_v
            elif data_2d == "pms":
                if any([_ is None for _ in (self.pmra, self.pmde)]):
                    raise ValueError("Data for  (pmra, pmde) data is required")
                c_str = "pms_c"
                x, y = self.pmra_v, self.pmde_v
            else:
                raise ValueError(f"Algorithm '{data_2d}' argument not recognized")

            # Remove rows containing any NaNs
            array_2d = np.array([x, y])
            x, y = array_2d[:, ~np.isnan(array_2d).any(axis=0)]

            x_c, y_c = cp.get_2D_center(x, y)

            print("\nCenter coordinates found:")
            print("{}        : ({:.4f}, {:.4f})".format(c_str, x_c, y_c))
            if c_str == "radec_c":
                self.radec_c = [x_c, y_c]
            if c_str == "pms_c":
                self.pms_c = [x_c, y_c]

        else:
            raise ValueError(f"Selected method '{algo}' not recognized")

    # def _get_radius(self, algo: str = "field_dens") -> None:
    #     """Estimate the cluster radius

    #     - ``ripley``: Originally introduced with the ``fastMP`` membership method
    #       in `Perren et al.
    #       (2023) <https://academic.oup.com/mnras/article/526/3/4107/7276628>`__.
    #       Requires ``(ra, dec)`` and their center estimates.
    #     - ``density``: Simple algorithm that counts the number of stars within the
    #       cluster region (center+radius) and subtracts the expected number of field
    #       stars within that region. Requires the number of cluster members to be
    #       defined.

    #     :param algo: Algorithm used to estimate center values, one of
    #         (``ripley, density``); defaults to ``ripley``
    #     :type algo: str

    #     :raises ValueError: If required attributes are  missing from the
    #         :py:class:`Cluster <asteca.cluster.Cluster>` object
    #     """

    #     if algo == "field_dens":
    #         for k in ("ra", "dec", "radec_c"):
    #             if hasattr(self, k) is False:
    #                 raise ValueError(f"'{k}' must be present as a 'cluster' attribute")
    #     elif algo == "king":
    #         for k in ("ra", "dec", "radec_c"):
    #             if hasattr(self, k) is False:
    #                 raise ValueError(f"'{k}' must be present as a 'cluster' attribute")

    #     xv, yv = self.ra_v, self.dec_v
    #     xy_center = self.radec_c

    #     if algo == "field_dens":
    #         radius = cp.fdens_radius(xv, yv, list(xy_center))
    #     elif algo == "king":
    #         radius = cp.king_radius(xv, yv, xy_center, self.N_cluster)
    #     else:
    #         raise ValueError(f"Unrecognized method '{algo}")

    #     print(f"\nRadius         : {radius}")
    #     self.radius = radius

    def get_nmembers(self, algo: str = "ripley", eq_to_gal: bool = True) -> None:
        """Estimate the number of members for the cluster

        - ``ripley``: Originally introduced with the ``fastMP`` membership method
          in `Perren et al.
          (2023) <https://academic.oup.com/mnras/article/526/3/4107/7276628>`__.
          Requires ``(ra, dec, pmra, pmde, plx)`` and their center estimates.
        - ``density``: Simple algorithm that counts the number of stars within the
          cluster region (center+radius) and subtracts the expected number of field
          stars within that region. Requires the `(RA, DEC)` center and the radius of
          the cluster to be defined.

        :param algo: Algorithm used to estimate center values, one of
            (``ripley, density``); defaults to ``ripley``
        :type algo: str
        :param eq_to_gal: Convert ``(RA, DEC)`` to ``(lon, lat)``. Useful for clusters
            with large ``DEC`` values to reduce the frame's distortion,
            defaults to ``False``
        :type eq_to_gal: bool

        :raises ValueError: If required attributes are  missing from the
            :py:class:`Cluster <asteca.cluster.Cluster>` object
        """
        if algo not in ("ripley", "density"):
            raise ValueError(f"'Argument algo={algo}' not recognized")

        xv, yv = self.ra_v, self.dec_v
        xy_center = self.radec_c
        # Convert (RA, DEC) to (lon, lat)
        if eq_to_gal is True:
            xv, yv = cp.radec2lonlat(xv, yv)
            xy_center = cp.radec2lonlat(*xy_center)

        if algo == "ripley":
            for k in ("ra", "dec", "pmra", "pmde", "plx", "radec_c", "pms_c", "plx_c"):
                if hasattr(self, k) is False:
                    raise ValueError(f"'{k}' must be present as a 'cluster' attribute")
            N_cluster = nm.ripley_nmembs(
                xv,
                yv,
                self.pmra_v,
                self.pmde_v,
                self.plx_v,
                list(xy_center),
                self.pms_c,
                self.plx_c,
            )
        elif algo == "density":
            for k in ("ra", "dec", "radec_c", "radius"):
                if hasattr(self, k) is False:
                    raise ValueError(f"'{k}' must be present as a 'cluster' attribute")
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
        print(f"\nN_cluster      : {N_cluster}")

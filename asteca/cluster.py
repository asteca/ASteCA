from dataclasses import dataclass
from typing import Optional
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from .modules import cluster_priv as cp


@dataclass
class cluster:
    """Define a :class:`cluster` object.

    :param cluster_df: `pandas DataFrame <https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html>`__
        with the cluster's loaded data.
    :type cluster_df: pd.DataFrame
    :param magnitude: Name of the DataFrame column that contains the magnitude
    :type magnitude: str
    :param e_mag: Name of the DataFrame column that contains the magnitude's uncertainty
    :type e_mag: str
    :param color: Name of the DataFrame column that contains the color
    :type color: str
    :param e_color: Name of the DataFrame column that contains the color's uncertainty
    :type e_color: str
    :param ra: Name of the DataFrame column that contains the right ascension (RA),
        defaults to ``None``
    :type ra: str, optional
    :param dec: Name of the DataFrame column that contains the declination (DEC),
        defaults to ``None``
    :type dec: str, optional
    :param plx: Name of the DataFrame column that contains the parallax,
        defaults to ``None``
    :type plx: str, optional
    :param pmra: Name of the DataFrame column that contains the RA proper motion,
        defaults to ``None``
    :type pmra: str, optional
    :param pmde: Name of the DataFrame column that contains the DEC proper motion,
        defaults to ``None``
    :type pmde: str, optional
    :param color2: Name of the DataFrame column that contains the second color,
        defaults to ``None``
    :type color2: str, optional
    :param e_color2: Name of the DataFrame column that contains the second color's
        uncertainty, defaults to ``None``
    :type e_color2: str, optional
    """

    cluster_df: pd.DataFrame
    magnitude: str
    e_mag: str
    color: str
    e_color: str
    ra: Optional[str] = None
    dec: Optional[str] = None
    plx: Optional[str] = None
    pmra: Optional[str] = None
    pmde: Optional[str] = None
    e_plx: Optional[str] = None
    e_pmra: Optional[str] = None
    e_pmde: Optional[str] = None
    color2: Optional[str] = None
    e_color2: Optional[str] = None

    def __post_init__(self):
        """Load photometry."""
        self._load()

    def _load(self):
        """The photometry is store with a '_p' to differentiate from the
        self.magnitude, self.color, etc that are defined with the class is called.
        """
        print("Instantiating cluster...")

        self.mag_p = np.array(self.cluster_df[self.magnitude])
        self.e_mag_p = np.array(self.cluster_df[self.e_mag])

        self.colors_p = [np.array(self.cluster_df[self.color])]
        if self.color2 is not None:
            self.colors_p.append(np.array(self.cluster_df[self.color2]))
        self.e_colors_p = [np.array(self.cluster_df[self.e_color])]
        if self.e_color2 is not None:
            self.e_colors_p.append(np.array(self.cluster_df[self.e_color2]))

        self.ra_v, self.dec_v = None, None
        if self.ra is not None:
            self.ra_v = self.cluster_df[self.ra]
            self.dec_v = self.cluster_df[self.dec]

        self.pmra_v, self.e_pmra_v, self.pmde_v, self.pmde_v = None, None, None, None
        if self.pmra is not None:
            self.pmra_v = self.cluster_df[self.pmra].values
            self.pmde_v = self.cluster_df[self.pmde].values
            self.e_pmra_v = self.cluster_df[self.e_pmra].values
            self.e_pmde_v = self.cluster_df[self.e_pmde].values

        self.plx_v, self.e_plx_v = None, None
        if self.plx is not None:
            self.plx_v = self.cluster_df[self.plx].values
            self.e_plx_v = self.cluster_df[self.e_plx].values

        print(f"N_stars        : {len(self.mag_p)}")
        print(f"Magnitude      : {self.magnitude}")
        print(f"Color          : {self.color}")
        if self.color2 is not None:
            print(f"Color2         : {self.color2}")
        print("Cluster object generated\n")

    def radecplot(self):
        """Generate a (RA, DEC) plot.

        :return: Matplotlib axis object
        :rtype: matplotlib.axis
        """
        ra = self.ra_v  # cluster_df[self.ra].values
        dec = self.dec_v  # cluster_df[self.dec].values
        mag = self.mag_p  # self.cluster_df[self.magnitude].values

        msk = ~np.isnan(mag)
        ra = ra[msk]
        dec = dec[msk]
        mag = mag[msk]

        # Mag to flux
        sizes = 10 ** (mag / -2.5)
        # Normalize to 1
        sizes /= sizes.max()
        # Set min, max
        sizes = 1 + 75 * sizes

        f, ax = plt.subplots()
        plt.scatter(ra, dec, s=sizes, c="k", alpha=0.7)
        plt.xlabel("RA")
        plt.ylabel("DEC")
        plt.gca().invert_xaxis()

        return ax

    def clustplot(self, ax, color_idx=0, binar_probs=None):
        """Generate a color-magnitude plot.

        :param ax: Matplotlib axis where to draw the plot
        :type ax: matplotlib.axis, optional, default=None
        :param color_idx: Index of the color to plot. If ``0`` (default), plot the
            first color. If ``1`` plot the second color.
        :type color_idx: int, default=0
        :param binar_probs: Array with probabilities of being a binary system for each
            observed star
        :type binar_probs: numpy.array, optional, default=None
        :return: Matplotlib axis object
        :rtype: matplotlib.axis
        """
        if color_idx > 1:
            raise ValueError(
                f"Wrong 'color_idx' value ({color_idx}), should be one of: [0, 1]"
            )

        if binar_probs is not None:
            msk_binar = binar_probs > 0.5

        mag_col = self.magnitude
        col_col = self.color
        if color_idx == 1:
            col_col = self.color2

        if binar_probs is None:
            ax.scatter(
                self.colors_p[color_idx],
                self.mag_p,
                c="green",
                alpha=0.5,
                label=f"Observed, N={len(self.mag_p)}",
            )
        else:
            ax.scatter(
                self.colors_p[color_idx][~msk_binar],
                self.mag_p[~msk_binar],
                # c="green",
                c=binar_probs[~msk_binar],
                marker="o",
                alpha=0.5,
                label=f"Observed (single), N={len(self.mag_p[~msk_binar])}",
            )
            ax.scatter(
                self.colors_p[color_idx][msk_binar],
                self.mag_p[msk_binar],
                # c="red",
                c=binar_probs[msk_binar],
                marker="s",
                alpha=0.5,
                label=f"Observed (binary), N={len(self.mag_p[msk_binar])}",
            )

        ax.set_ylim(max(self.mag_p) + 0.5, min(self.mag_p) - 0.5)
        ax.set_xlabel(col_col)
        ax.set_ylabel(mag_col)
        ax.legend()

        return ax

    def get_center(
        self,
        method="knn_5d",
        xy_c=None,
        vpd_c=None,
        plx_c=None,
        N_cluster=None,
        N_clust_min=25,
    ):
        """Estimate center coordinates for the cluster

        Use the available data (ra, dec, pmra, pmde, plx) to estimate a cluster's
        center coordinates as the point of maximum density. Methods:

        ``knn_5d``: Estimates the center value as the median position of the k
        (k=N_clust_min) nearest stars to an estimate of the center in proper motions
        and (ra, dec, plx), if given. All 5 dimensions of data must be available.

        :param method: Method used to estimate center values, one of (``knn_5d``),
            defaults to ``knn_5d``
        :type method: str
        :param xy_c: Estimated value for the (RA, DEC) center, defaults to ``None``
        :type xy_c: tuple(float, float), optional
        :param vpd_c: Estimated value for the (pmRA, pmDE) center, defaults to ``None``
        :type vpd_c: tuple(float, float), optional
        :param plx_c: Estimated value for the plx center, defaults to ``None``
        :type plx_c: float, optional
        :param N_cluster: Estimated number of members, defaults to ``None``
        :type N_cluster: int, optional
        :param N_clust_min: Minimum number of cluster members, defaults to ``25``
        :type N_clust_min: int, optional

        :return: Dictionary of center coordinates
        :rtype: dict
        """

        if method == "knn_5d":
            if any(
                [
                    _ is None
                    for _ in (
                        self.pmra_v,
                        self.pmde_v,
                        self.plx_v,
                        self.e_pmra_v,
                        self.e_pmde_v,
                        self.e_plx_v,
                    )
                ]
            ):
                raise ValueError(
                    f"Method '{method}' requires (ra, dec, pmra, pmde, plx) data and "
                    + "their uncertainties"
                )

            # To galactic coordinates
            glon, glat = cp.radec2lonlat(self.ra_v.values, self.dec_v.values)
            X = np.array(
                [
                    glon,
                    glat,
                    self.pmra_v,
                    self.pmde_v,
                    self.plx_v,
                    self.e_pmra_v,
                    self.e_pmde_v,
                    self.e_plx_v,
                ]
            )
            # Reject nan values and extract clean data
            _, X_no_nan = cp.reject_nans(X)
            lon, lat, pmRA, pmDE, plx, _, _, _ = X_no_nan

            x_c, y_c, pmra_c, pmde_c, plx_c = cp.get_5D_center(
                lon,
                lat,
                pmRA,
                pmDE,
                plx,
                xy_c=xy_c,
                vpd_c=vpd_c,
                plx_c=plx_c,
                N_cluster=N_cluster,
                N_clust_min=N_clust_min,
            )
            ra_c, dec_c = cp.lonlat2radec(x_c, y_c)

        print("\nCenter coordinates found:")
        print("ra_c     : {:.4f}".format(ra_c))
        print("dec_c    : {:.4f}".format(dec_c))
        print("pmra_c   : {:.3f}".format(pmra_c))
        print("pmde_c   : {:.3f}".format(pmde_c))
        print("plx_c    : {:.3f}".format(plx_c))

        return {
            "ra_c": ra_c,
            "dec_c": dec_c,
            "pmra_c": pmra_c,
            "pmde_c": pmde_c,
            "plx_c": plx_c,
        }

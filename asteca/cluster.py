from dataclasses import dataclass
from typing import Optional
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


@dataclass
class cluster:
    r"""Define a ``cluster`` object.

    Parameters
    ----------
    cluster_df : pd.DataFrame
        pandas DataFrame with the cluster's loaded data
    magnitude : str
        Name of the DataFrame column that contains the magnitude
    e_mag : str
        Name of the DataFrame column that contains the magnitude's uncertainty
    color : str
        Name of the DataFrame column that contains the color
    e_color : str
        Name of the DataFrame column that contains the color's uncertainty
    ra: str, optional, default=None
        Name of the DataFrame column that contains the right ascension (RA)
    dec: str, optional, default=None
        Name of the DataFrame column that contains the declination (DEC)
    plx: str, optional, default=None
        Name of the DataFrame column that contains the parallax
    pmra: str, optional, default=None
        Name of the DataFrame column that contains the RA proper motion
    pmde: str, optional, default=None
        Name of the DataFrame column that contains the DEC proper motion
    color2: str, optional, default=None
        Name of the DataFrame column that contains the second color
    e_color2: str, optional, default=None
        Name of the DataFrame column that contains the second color's uncertainty

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
    color2: Optional[str] = None
    e_color2: Optional[str] = None

    def __post_init__(self):
        # Load photometry
        self._load()

    def _load(self):
        """
        The photometry is store with a '_p' to differentiate from the self.magnitude,
        self.color, etc that are defined with the class is called.
        """
        print("Reading and processing cluster data")

        self.mag_p = np.array(self.cluster_df[self.magnitude])
        self.e_mag_p = np.array(self.cluster_df[self.e_mag])

        self.colors_p = [np.array(self.cluster_df[self.color])]
        if self.color2 is not None:
            self.colors_p.append(np.array(self.cluster_df[self.color2]))
        self.e_colors_p = [np.array(self.cluster_df[self.e_color])]
        if self.e_color2 is not None:
            self.e_colors_p.append(np.array(self.cluster_df[self.e_color2]))

        if self.ra is not None:
            self.ra_v = self.cluster_df[self.ra]
            self.dec_v = self.cluster_df[self.dec]

        print("Cluster object generated\n")

    def radecplot(self):
        r"""Generate a (RA, DEC) plot.

        Returns
        -------
        matplotlib.axis
            Matplotlib axis object

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

    def clustplot(self, ax=None, binar_prob=None):
        r"""Generate a color-magnitude plot.

        Parameters
        ----------
        ax : matplotlib.axis, optional, default=None
            Matplotlib axis where to draw the plot
        binar_prob : numpy.array, optional, default=None
            Array with probabilities of being a binary system for each observed star

        Returns
        -------
        matplotlib.axis
            Matplotlib axis object


        """
        if ax is None:
            f, ax = plt.subplots()

        if binar_prob is not None:
            msk_binar = binar_prob > 0.5

        mag_col = self.magnitude
        col_col = self.color

        if binar_prob is None:
            ax.scatter(
                self.colors_p[0],
                self.mag_p,
                c="green",
                alpha=0.5,
                label=f"Observed, N={len(self.mag_p)}",
            )
        else:
            ax.scatter(
                self.colors_p[0][~msk_binar],
                self.mag_p[~msk_binar],
                c="green",
                alpha=0.5,
                label=f"Observed (single), N={len(self.mag_p[~msk_binar])}",
            )
            ax.scatter(
                self.colors_p[0][msk_binar],
                self.mag_p[msk_binar],
                c="red",
                alpha=0.5,
                label=f"Observed (binary), N={len(self.mag_p[msk_binar])}",
            )

        ax.set_ylim(max(self.mag_p) + .5, min(self.mag_p) - .5)
        ax.set_xlabel(col_col)
        ax.set_ylabel(mag_col)
        ax.legend()

        return ax

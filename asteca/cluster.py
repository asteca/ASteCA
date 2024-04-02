from dataclasses import dataclass
from typing import Optional
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#
from .modules import cluster_priv
from .modules import synth_cluster_priv as scp


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
    source_id: str, optional, default=None
        Name of the DataFrame column that contains the IDs for the stars
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
    # N_models: int = 1000
    source_id: Optional[str] = None
    ra: Optional[str] = None
    dec: Optional[str] = None
    plx: Optional[str] = None
    pmra: Optional[str] = None
    pmde: Optional[str] = None
    color2: Optional[str] = None
    e_color2: Optional[str] = None

    def __post_init__(self):

        # Load photometry
        cluster_priv.load(self)

    def radecplot(self):
        r"""Generate a (RA, DEC) plot.

        Returns
        -------
        matplotlib.axis
            Matplotlib axis object

        """
        ra = self.cluster_df[self.ra].values
        dec = self.cluster_df[self.dec].values
        mag = self.cluster_df[self.magnitude].values

        msk = ~np.isnan(mag)
        ra = ra[msk]
        dec = dec[msk]
        mag = mag[msk]

        # Mag to flux
        sizes = 10**(mag/-2.5)
        # Normalize to 1
        sizes /= sizes.max()
        # Set min, max
        sizes = 1 + 75*sizes

        f, ax = plt.subplots()
        plt.scatter(ra, dec, s=sizes, c="k", alpha=.7)
        plt.xlabel("RA")
        plt.ylabel("DEC")
        plt.gca().invert_xaxis()

        return ax

    def clustplot(self, ax=None):
        r"""Generate a color-magnitude plot.

        Parameters
        ----------
        ax : matplotlib.axis, optional, default=None
            Matplotlib axis where to draw the plot

        Returns
        -------
        matplotlib.axis
            Matplotlib axis object


        """
        invert_yaxis_flag = False
        if ax is None:
            f, ax = plt.subplots()
            invert_yaxis_flag = True

        mag_col = self.magnitude
        col_col = self.color

        plt.scatter(
            self.colors_p[0],
            self.mag_p,
            c="green",
            alpha=0.5,
            label=f"Observed, N={len(self.mag_p)}",
        )

        plt.xlabel(col_col)
        plt.ylabel(mag_col)
        plt.legend()

        if invert_yaxis_flag:
            plt.gca().invert_yaxis()

        return ax

    def _masses_binar(self, synthcl, fit_params, model_std):
        """
        Assign individual masses to the observed cluster's stars, along with binary
        probabilities (if binarity was estimated).

        Estimate the statistics for the mass and binary fractions (not fitted)
        """
        print("Estimating total mass, binary probabilities, and individual star masses")
        cl_ids = np.array(self.cluster_df[self.source_id])
        ra = np.median(self.cluster_df[self.ra])
        dec = np.median(self.cluster_df[self.dec])

        m_ini_idx = self.m_ini_idx
        # Extract photometry used in the best fit process
        mags_cols_cl = [self.mag] + [_ for _ in self.colors]
        obs_phot = np.array(mags_cols_cl)

        # Generate random models from the selected solution
        models = cluster_priv.ranModels(self.N_models, fit_params, model_std)

        # Identify nans in mag and colors and re-write them as -10
        nan_msk = np.full(obs_phot.shape[1], False)
        for col in obs_phot:
            nan_msk = nan_msk | np.isnan(col)
        obs_phot[:, nan_msk] = -10.0
        obs_phot = obs_phot.T
        N_obs = obs_phot.shape[0]

        # Initiate empty arrays for mean and variance
        st_mass_mean, st_mass_var = np.zeros(N_obs), np.zeros(N_obs)
        st_mass_mean_binar, st_mass_var_binar = np.zeros(N_obs), np.zeros(N_obs)
        prob_binar = np.zeros(N_obs)
        binar_vals = []
        Nm, Nm_binar = 0, 0

        # from .synth_cluster import invTrnsfSmpl
        # inv_cdf = invTrnsfSmpl(synthcl.IMF_name)
        try:
            st_dist_mass = synthcl.st_dist_mass_full
        except (KeyError, AttributeError):
            st_dist_mass = synthcl.st_dist_mass

        masses_dict = {'M_actual': [], 'M_init': []}
        for _, model in enumerate(models):
            # Generate synthetic cluster from the 'model'.
            # isoch = synthcl.generate(model, self, True)
            isoch = scp.generate(synthcl, model)
            if not isoch.any():
                continue

            # TODO
            Nm, Nm_binar, st_mass_mean, st_mass_var, st_mass_mean_binar, \
                st_mass_var_binar, binar_vals, prob_binar = cluster_priv.xxx(
                Nm, st_mass_mean, st_mass_var, Nm_binar, obs_phot, m_ini_idx, st_mass_mean_binar,
                st_mass_var_binar, prob_binar, binar_vals, synthcl.alpha, isoch)

            # Extract dm and loga
            model_comb = synthcl.fix_params | model
            loga, dm = model_comb['loga'], model_comb['dm']
            #
            masses_dict = cluster_priv.get_masses(
                masses_dict, ra, dec, m_ini_idx, st_dist_mass, isoch, loga, dm)

        # Store standard deviations instead of variances
        st_mass_std = np.sqrt(st_mass_var / Nm)
        st_mass_std_binar = np.sqrt(st_mass_var_binar / max(1, Nm_binar))

        # if synthcl.max_mass < np.median(M_init_arr) + np.std(M_init_arr):
        #     logging.warning(
        #         "The total mass is too close to the 'synth_clusters.max_mass' "
        #         + "parameter. Consider increasing this value."
        #     )

        # Use nan values for stars with initial nan photometric values
        for arr in (
            st_mass_mean,
            st_mass_std,
            st_mass_mean_binar,
            st_mass_std_binar,
            prob_binar,
        ):
            arr[nan_msk] = np.nan

        cl_masses_bfr = pd.DataFrame(
            {
                "ID": cl_ids,
                "M1": st_mass_mean,
                "M1_std": st_mass_std,
                "M2": st_mass_mean_binar,
                "M2_std": st_mass_std_binar,
                "binar_prob": prob_binar,
            }
        )

        return masses_dict, np.array(binar_vals), cl_masses_bfr

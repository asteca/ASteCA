from dataclasses import dataclass
from typing import Optional
import numpy as np
import pandas as pd
from .modules import cluster_priv
from .modules import synth_cluster_priv as scp


@dataclass
class cluster:
    source_id: str
    magnitude: str
    e_mag: str
    color: str
    e_color: str
    model_fixed: dict
    cluster_df: Optional[None] = None
    ra: Optional[str] = None
    dec: Optional[str] = None
    bin_method: str = "knuth"
    N_mag: int = 15
    N_col: int = 10
    N_models: int = 100
    color2: str | None = None
    e_color2: str | None = None

    def __post_init__(self):
        bin_methods = ("knuth", "fixed", "bayes_blocks")
        if self.bin_method not in bin_methods:
            raise ValueError(
                f"Binning '{self.bin_method}' not recognized. "
                + f"Should be one of {bin_methods}"
            )

        self.cluster_dict = cluster_priv.load(self)

    def masses_binar(self, synthcl, model_fit, model_std):
        """
        Assign individual masses to the observed cluster's stars, along with binary
        probabilities (if binarity was estimated).

        Estimate the statistics for the mass and binary fractions (not fitted)
        """
        print("Estimating total mass, binary probabilities, and individual star masses")
        cl_dict = self.cluster_dict
        ra = np.median(self.cluster_df[self.ra])
        dec = np.median(self.cluster_df[self.dec])

        model_fixed = self.model_fixed
        m_ini_idx = cl_dict["m_ini_idx"]
        # Extract photometry used in the best fit process
        mags_cols_cl = [cl_dict["mag"]] + [_ for _ in cl_dict["colors"]]
        obs_phot = np.array(mags_cols_cl)

        # Generate random models from the selected solution
        models = cluster_priv.ranModels(self.N_models, model_fit, model_std)

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
        except KeyError:
            st_dist_mass = synthcl.st_dist_mass

        masses_dict = {'M_actual': [], 'M_init': []}
        for _, model in enumerate(models):
            # Generate synthetic cluster from the 'model'.
            # isoch = synthcl.generate(model, self, True)
            isoch = scp.generate(synthcl, self, model)
            if not isoch.any():
                continue

            # TODO
            Nm, Nm_binar, st_mass_mean, st_mass_var, st_mass_mean_binar, \
                st_mass_var_binar, binar_vals, prob_binar = cluster_priv.xxx(
                Nm, st_mass_mean, st_mass_var, Nm_binar, obs_phot, m_ini_idx, st_mass_mean_binar,
                st_mass_var_binar, prob_binar, binar_vals, synthcl.alpha, isoch)

            # Extract dm and loga
            model_comb = model_fixed | model
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
                "ID": cl_dict["cl_ids"],
                "M1": st_mass_mean,
                "M1_std": st_mass_std,
                "M2": st_mass_mean_binar,
                "M2_std": st_mass_std_binar,
                "binar_prob": prob_binar,
            }
        )

        return masses_dict, np.array(binar_vals), cl_masses_bfr

    def clust_plot(self, synthcl, model_fit):
        """ """
        import matplotlib.pyplot as plt

        mag_col = self.magnitude
        col_col = self.color
        cl_dict = self.cluster_dict

        # Generate synthetic cluster.
        synth_clust = scp.generate(synthcl, self, model_fit)
        # plt.title(f"Lkl dist ={round(final_dist, 3)}")
        # y_edges, x_edges = cluster_dict["bin_edges"]
        # for xe in x_edges:
        #     plt.axvline(xe, c="grey", ls=":")
        # for ye in y_edges:
        #     plt.axhline(ye, c="grey", ls=":")

        # plt.figure(figsize=(8, 8), dpi=300)

        plt.scatter(
            cl_dict["colors"][0],
            cl_dict["mag"],
            c="green",
            alpha=0.5,
            label=f"Observed, N={cl_dict['N_obs_stars']}",
        )

        if synthcl.alpha is not None:
            binar_idx = synth_clust[-1] != 0.0
        else:
            binar_idx = np.full(synth_clust.shape[1], False)

        x_synth, y_synth = synth_clust[1], synth_clust[0]
        # Single synthetic systems
        plt.scatter(
            x_synth[~binar_idx],
            y_synth[~binar_idx],
            marker="^",
            c="#519ddb",
            alpha=0.5,
            label=f"Single, N={len(x_synth[~binar_idx])}",
        )
        # Binary synthetic systems
        plt.scatter(
            x_synth[binar_idx],
            y_synth[binar_idx],
            marker="v",
            c="#F34C4C",
            alpha=0.5,
            label=f"Binary, N={len(x_synth[binar_idx])}",
        )

        plt.xlabel(col_col)
        plt.ylabel(mag_col)

        plt.gca().invert_yaxis()
        plt.legend()
        plt.show()
        # plt.savefig(f"{fname}.png")

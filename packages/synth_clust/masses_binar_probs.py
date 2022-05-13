
import logging
import numpy as np
from scipy.spatial import cKDTree
from ..best_fit.bf_common import ranModels, getSynthClust
from ..best_fit.prep_obs_params import dataProcess
from ..aux_funcs import kde1D
from .. import update_progress


def main(clp, pd, td):
    """
    Assign masses to the (decontaminated) observed cluster, and binary
    probabilities (if binarity was estimated).

    Estimate the statistics for the mass and binary fractions (not fitted)

    Also generate the uncertainty region (for plotting).
    """
    # No best fit process was employed
    if pd['best_fit_algor'] == 'n':
        return clp

    # Generate random models from the selected solution (mean, median, mode,
    # MAP), given by 'D3_sol.
    models = ranModels(
        td['fundam_params'], pd['D3_sol'], clp['isoch_fit_params'],
        clp['isoch_fit_errors'])

    print("Estimating binary probabilities and masses")
    # Extract photometry used in the best fit process
    mags_cols_cl, _ = dataProcess(clp['cl_syn_fit'])
    # Arrange properly
    mags, cols = [np.array(_) for _ in mags_cols_cl]
    obs_phot = np.concatenate([mags, cols]).T

    # Initiate empty arrays for mean and variance
    st_mass_mean, M2 = np.zeros(obs_phot.shape[0]), np.zeros(obs_phot.shape[0])
    st_mass_mean_binar, M2_binar = np.zeros(obs_phot.shape[0]),\
        np.zeros(obs_phot.shape[0])
    prob_binar = np.zeros(obs_phot.shape[0])

    # Estimate the mean and variance for each star via recurrence.
    Nm_binar, binar_vals, MassT_vals = 0, [], []
    # Used for plotting
    synthcl_Nsigma = [[] for _ in range(sum(td['N_fc']))]
    for Nm, model in enumerate(models):

        # Generate synthetic cluster from the 'model'.
        isoch, M_total = getSynthClust(model, False, clp['syntClustArgs'])
        if not isoch.any():
            continue

        for i, photom_dim in enumerate(isoch[:sum(td['N_fc'])]):
            synthcl_Nsigma[i] += list(photom_dim)

        # Masses, binary mask
        mass_primary = isoch[td['m_ini_idx']]
        # Binaries have M2>0
        binar_idxs = isoch[-1] > 0.
        mass_secondary = isoch[-1]
        # shape: (N_stars, Ndim)
        photom = isoch[:sum(td['N_fc'])].T

        MassT_vals.append(M_total)
        binar_vals.append(binar_idxs.sum() / isoch.shape[-1])

        # For non-binary systems
        photom_single = photom[~binar_idxs]
        if photom_single.any():
            obs_mass, lkl_p = photomMatch(
                obs_phot, photom_single, mass_primary[~binar_idxs])
            # Estimate mean and variance
            st_mass_mean, M2 = recurrentStats(Nm, st_mass_mean, M2, obs_mass)

            # For binary systems
            if binar_idxs.sum() > 0:
                photom_binar = photom[binar_idxs]
                # If there are no binary systems, skip
                if photom_binar.any():
                    Nm_binar += 1
                    obs_mass, lkl_b = photomMatch(
                        obs_phot, photom_binar, mass_secondary[binar_idxs])
                    st_mass_mean_binar, M2_binar = recurrentStats(
                        Nm, st_mass_mean_binar, M2_binar, obs_mass)

                    # Bayesian probability
                    new_prob_binar = 1. / (1. + (lkl_p / lkl_b))
                    prob_binar = recurrentStats(
                        Nm, prob_binar, None, new_prob_binar)

        update_progress.updt(models.shape[0], Nm + 1)

    # Used for plotting
    synthcl_Nsigma = np.array(synthcl_Nsigma)

    # Store standard deviations
    st_mass_std = np.sqrt(M2 / Nm)
    st_mass_std_binar = np.sqrt(M2_binar / max(1, Nm_binar))

    MassT_dist_vals, binar_dist_vals = estimMassBinar(
        pd, MassT_vals, binar_vals)

    clp['st_mass_mean'], clp['st_mass_std'], clp['st_mass_mean_binar'],\
        clp['st_mass_std_binar'], clp['prob_binar'], clp['MassT_dist_vals'],\
        clp['binar_dist_vals'], clp['synthcl_Nsigma'] = st_mass_mean,\
        st_mass_std, st_mass_mean_binar, st_mass_std_binar, prob_binar,\
        MassT_dist_vals, binar_dist_vals, synthcl_Nsigma

    return clp


def photomMatch(obs_phot, photom, mass_ini):
    """
    For each observed star in 'obs_phot', find the closest synthetic star in
    the (synthetic) photometric space 'photom'
    """
    tree = cKDTree(photom)
    dd, ii = tree.query(obs_phot, k=1)

    # Assign masses to each observed star
    obs_mass = mass_ini[ii]

    # Likelihood is defined as the inverse of the distance
    lkl = 1. / dd

    return obs_mass, lkl


def recurrentStats(Nm, mean, var, newValue):
    """
    Source: en.wikipedia.org/wiki/
            Algorithms_for_calculating_variance#Welford's_online_algorithm
    """
    count = Nm + 1
    delta = newValue - mean
    mean += delta / count
    if var is None:
        return mean
    var += delta * (newValue - mean)
    return mean, var


def estimMassBinar(pd, MassT_vals, binar_vals):
    """
    Estimate the mean, etc values for the mass and binary fraction
    distributions, as these are not fitted parameters.
    """
    x_kde, M_kde = kde1D(MassT_vals)
    mode_M = x_kde[np.argmax(M_kde)]
    MassT_dist_vals = {
        'mean_sol': np.mean(MassT_vals), 'median_sol': np.median(MassT_vals),
        'mode_sol': mode_M, 'errors': (
            np.percentile(MassT_vals, 16), np.percentile(MassT_vals, 84),
            np.std(MassT_vals))
    }

    if pd['Max_mass'] < np.mean(MassT_vals) + np.std(MassT_vals):
        logging.warning(
            "The total mass is too close to the 'Max_mass' "
            + "parameter. Consider increasing this value.")

    # The binary fraction will have a small dispersion even for fixed 'beta'
    # because the probabilities depend on the masses which vary for different
    # values of (z, a)
    x_kde, bf_kde = kde1D(binar_vals)
    mode_bf = x_kde[np.argmax(bf_kde)]
    binar_dist_vals = {
        'mean_sol': np.mean(binar_vals), 'median_sol': np.median(binar_vals),
        'mode_sol': mode_bf, 'errors': (
            np.percentile(binar_vals, 16), np.percentile(binar_vals, 84),
            np.std(binar_vals))
    }

    return MassT_dist_vals, binar_dist_vals

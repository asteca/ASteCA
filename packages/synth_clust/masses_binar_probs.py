
import numpy as np
from scipy.spatial import cKDTree
from .synthClustPrep import setSynthClust
from ..best_fit.obs_clust_prepare import dataProcess
from .. import update_progress


def main(clp, pd):
    """
    Assign masses to the (decontaminated) observed cluster, and binary
    probabilities (if binarity was estimated).
    """

    # Dummy arrays
    clp['st_mass_mean'], clp['st_mass_std'],\
        clp['st_mass_mean_binar'], clp['st_mass_std_binar'],\
        clp['prob_binar'] = [np.array([]) for _ in range(5)]
    # No best fit process was employed
    if pd['best_fit_algor'] == 'n':
        return clp

    # Generate random models from the selected solution (mean, median, mode,
    # MAP), given by 'D3_sol.
    models = ranModels(
        pd['fundam_params'], pd['D3_sol'], clp['isoch_fit_params'],
        clp['isoch_fit_errors'])

    if not models.any():
        print(" WARNING: could not assign masses and binary probabilities")
        return clp

    print("Estimating binary probabilities and masses")
    # Extract photometry used in the best fit process
    mags_cols_cl, _ = dataProcess(clp['cl_max_mag'])
    # Arrange properly
    mags, cols = [np.array(_) for _ in mags_cols_cl]
    obs_phot = np.concatenate([mags, cols]).T

    # Initiate empty arrays for mean and variance
    st_mass_mean, M2 = np.zeros(obs_phot.shape[0]), np.zeros(obs_phot.shape[0])
    st_mass_mean_binar, M2_binar = np.zeros(obs_phot.shape[0]),\
        np.zeros(obs_phot.shape[0])
    prob_binar = np.zeros(obs_phot.shape[0])

    # Estimate the mean and variance for each star via recurrence.
    Nm_binar = 0
    for Nm, model in enumerate(models):

        # Generate synthetic cluster from the 'model'.
        isoch = setSynthClust(model, *clp['syntClustArgs'])
        if not isoch.any():
            continue

        # Masses, binary mask
        mass_primary = isoch[pd['m_ini_idx']]
        binar_idxs = ~(isoch[-1] == -99.)
        mass_secondary = isoch[-1]
        mass_binar = mass_primary + mass_secondary
        # shape: (N_stars, Ndim)
        photom = isoch[:sum(pd['N_fc'])].T

        # For non-binary systems
        photom_single = photom[~binar_idxs]
        if photom_single.any():
            obs_mass, lkl_p = photomMatch(
                obs_phot, photom_single, mass_primary[~binar_idxs])
            # Estimate mean and variance
            st_mass_mean, M2 = recurrentStats(Nm, st_mass_mean, M2, obs_mass)

        # For binary systems
        if pd['binar_flag']:
            photom_binar = photom[binar_idxs]
            # If there are no binary systems, skip
            if photom_binar.any():
                Nm_binar += 1
                obs_mass, lkl_b = photomMatch(
                    obs_phot, photom_binar, mass_binar[binar_idxs])
                st_mass_mean_binar, M2_binar = recurrentStats(
                    Nm, st_mass_mean_binar, M2_binar, obs_mass)

                # Bayesian probability
                new_prob_binar = 1. / (1. + (lkl_p / lkl_b))
                prob_binar = recurrentStats(
                    Nm, prob_binar, None, new_prob_binar)

        update_progress.updt(models.shape[0], Nm + 1)

    # Store standard deviations
    st_mass_std = np.sqrt(M2 / Nm)
    st_mass_std_binar = np.sqrt(M2_binar / max(1, Nm_binar))

    clp['st_mass_mean'], clp['st_mass_std'], clp['st_mass_mean_binar'],\
        clp['st_mass_std_binar'], clp['prob_binar'] = st_mass_mean,\
        st_mass_std, st_mass_mean_binar, st_mass_std_binar, prob_binar

    return clp


def ranModels(fundam_params, D3_sol, isoch_fit_params, isoch_fit_errors,
              N_models=1000):
    """
    Generate the requested models via sampling a Gaussian centered on the
    selected solution, with standard deviation given by the attached
    uncertainty.

    HARDCODED:

    N_models: number of models to generate.
    """
    # Use the selected solution values for all the parameters.
    model = isoch_fit_params[D3_sol + '_sol']

    # Extract standard deviations.
    p_vals, nancount = [], 0
    for i, p in enumerate(model):
        std = isoch_fit_errors[i][-1]
        if not np.isnan(std):
            p_vals.append([
                p, std, min(fundam_params[i]), max(fundam_params[i])])
        else:
            # The parameter has no uncertainty attached
            nancount += 1

    # Check if at least one parameter has an uncertainty attached.
    if nancount < 6:
        # Generate 'N_models' random models.
        models = []
        for par in p_vals:
            model = np.random.normal(par[0], par[1], N_models)
            model = np.clip(model, a_min=par[2], a_max=par[3])
            models.append(model)
        models = np.array(models).T
    else:
        models = np.array([])

    return models


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

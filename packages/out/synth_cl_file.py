
import numpy as np
from scipy.spatial import cKDTree
from astropy.table import Table
from ..inp import data_IO
from ..synth_clust import synth_cluster
from ..best_fit.obs_clust_prepare import dataProcess


def main(clp, npd, best_fit_algor, fundam_params, filters, colors,
         theor_tracks, R_V, m_ini_idx, binar_flag, D3_sol, **kwargs):
    """
    This function serves three purposes:

    1. Create output data file with stars in the best fit synthetic cluster
    found by the 'Best Fit' function.
    2. Generate the N-sigma uncertainty region for plotting
    3. Assign per-star masses
    """

    clp['synth_clst_plot'], clp['binar_idx_plot'], clp['shift_isoch'],\
        clp['synthcl_Nsigma'], clp['st_mass_mean'], clp['st_mass_std'],\
        clp['st_mass_mean_binar'], clp['st_mass_std_binar'],\
        clp['prob_binar'] = [np.array([]) for _ in range(9)]

    if best_fit_algor != 'n':
        isoch_moved, synth_clust, sigma, extra_pars, synthcl_Nsigma,\
            st_mass_mean, st_mass_std, st_mass_mean_binar, st_mass_std_binar,\
            prob_binar = synth_cl_plot(
                fundam_params, clp['isoch_fit_params'],
                clp['isoch_fit_errors'], theor_tracks, clp['completeness'],
                clp['max_mag_syn'], clp['st_dist_mass'], R_V, clp['ext_coefs'],
                clp['N_fc'], clp['err_pars'], m_ini_idx, binar_flag,
                clp['cl_max_mag'], D3_sol)

        # If cluster is not empty.
        if synth_clust.any():
            writeFileOut(npd, filters, colors, synth_clust, sigma, extra_pars,
                         clp['cl_max_mag'], st_mass_mean, st_mass_std,
                         st_mass_mean_binar, st_mass_std_binar, prob_binar)
            # Save for plotting.
            binar_idx = extra_pars[0]
            clp['synth_clst_plot'], clp['binar_idx_plot'],\
                clp['shift_isoch'], clp['synthcl_Nsigma'] = synth_clust,\
                binar_idx, isoch_moved, synthcl_Nsigma
        else:
            print("  ERROR: empty synthetic cluster could not be saved\n"
                  "  to file")

    return clp


def synth_cl_plot(
    fundam_params, isoch_fit_params, isoch_fit_errors, theor_tracks,
    completeness, max_mag_syn, st_dist_mass, R_V, ext_coefs, N_fc, err_pars,
        m_ini_idx, binar_flag, cl_max_mag, D3_sol):
    """
    Prepare data and generate the uncertainty region and individual masses,
    if uncertainties for at least one parameter exists.
    """

    # Pack common args.
    syntClustArgs = (
        fundam_params, isoch_fit_params['varIdxs'], theor_tracks,
        completeness, max_mag_syn, st_dist_mass, R_V, ext_coefs,
        N_fc, err_pars, m_ini_idx, binar_flag)

    # Use the selected solution values for all the parameters.
    model = isoch_fit_params[D3_sol + '_sol']

    # Generate isochrone, synthetic cluster (with uncertainties), and extra
    # parameters for the "best" fitted parameters.
    # Model with the "best" fitted parameters.
    model_var = np.array(model)[isoch_fit_params['varIdxs']]
    isoch_moved, synth_clust, sigma, extra_pars = setSynthClust(
        model_var, True, *syntClustArgs)

    # Extract standard deviations.
    p_vals, nancount = [], 0
    for i, p in enumerate(model):
        std = isoch_fit_errors[i][-1]
        if not np.isnan(std):
            vals = [p, std]
            p_vals.append(vals)
        else:
            nancount += 1

    # Check if any parameter has an uncertainty attached. Else, skip process
    synthcl_Nsigma, st_mass_mean, st_mass_std, st_mass_mean_binar,\
        st_mass_std_binar, prob_binar = [np.array([]) for _ in range(6)]
    if nancount != 6:
        synthcl_Nsigma, st_mass_mean, st_mass_std, st_mass_mean_binar,\
            st_mass_std_binar, prob_binar = ranModels(
                binar_flag, N_fc, cl_max_mag, syntClustArgs, p_vals)

    return isoch_moved, synth_clust, sigma, extra_pars, synthcl_Nsigma,\
        st_mass_mean, st_mass_std, st_mass_mean_binar, st_mass_std_binar,\
        prob_binar


def ranModels(
        binar_flag, N_fc, cl_max_mag, syntClustArgs, p_vals, N_models=1000):
    """
    Generate a synthetic cluster and a N(mu, sigma) region for plotting, from
    'N_models' models.

    Also assigns individual masses to stars in the true member subset.
    """
    print("Estimating individual masses")

    # Generate 'N_models' random models.
    models = []
    for par in p_vals:
        models.append(np.random.normal(par[0], par[1], N_models))

    # Extract photometry used in the best fit process
    mags_cols_cl, _ = dataProcess(cl_max_mag)
    # Arrange properly
    mags, cols = [np.array(_) for _ in mags_cols_cl]
    obs_phot = np.concatenate([mags, cols]).T

    # Initiate empty arrays for mean and variance
    st_mass_mean, M2 = np.zeros(obs_phot.shape[0]), np.zeros(obs_phot.shape[0])
    st_mass_mean_binar, M2_binar = np.zeros(obs_phot.shape[0]),\
        np.zeros(obs_phot.shape[0])
    prob_binar = np.zeros(obs_phot.shape[0])

    # Generate the synthetic clusters from the sampled parameters.
    synthcl_Nsigma = [[] for _ in range(sum(N_fc))]
    for Nm, model in enumerate(np.array(models).T):

        # Estimate the mean and variance for each star via recurrence.
        isoch_mv, synth_cl = setSynthClust(model, True, *syntClustArgs)[:2]

        # TODO assuming that m_ini masses are positioned last in this array
        mass_ini = isoch_mv[-1]

        # For non-binary systems
        photom = isoch_mv[:sum(N_fc)].T
        obs_mass, lkl_p = photomMatch(obs_phot, mass_ini, photom)
        # Estimate mean and variance
        st_mass_mean, M2 = recurrentStats(Nm, st_mass_mean, M2, obs_mass)

        # For binary systems
        if binar_flag:
            photom_binar = isoch_mv[sum(N_fc):2 * sum(N_fc)].T
            obs_mass, lkl_b = photomMatch(obs_phot, mass_ini, photom_binar)
            st_mass_mean_binar, M2_binar = recurrentStats(
                Nm, st_mass_mean_binar, M2_binar, obs_mass)

            # Bayesian probability
            new_prob_binar = 1. / (1. + (lkl_p / lkl_b))
            prob_binar = recurrentStats(Nm, prob_binar, None, new_prob_binar)

        # Synthetic cluster
        sc = synth_cl.T
        if sc.any():
            for i, photom_dim in enumerate(sc):
                synthcl_Nsigma[i] += list(photom_dim)
    synthcl_Nsigma = np.array(synthcl_Nsigma)

    # Store standard deviations
    st_mass_std = np.sqrt(M2 / Nm)
    st_mass_std_binar = np.sqrt(M2_binar / Nm)

    return synthcl_Nsigma, st_mass_mean, st_mass_std, st_mass_mean_binar,\
        st_mass_std_binar, prob_binar


def photomMatch(obs_phot, mass_ini, photom):
    """
    For each observed star, find the closest synthetic star in the
    photometric space 'photom'
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


def setSynthClust(
    model, extra_pars_flag, fundam_params, varIdxs, theor_tracks,
    completeness, max_mag_syn, st_dist_mass, R_V, ext_coefs, N_fc, err_pars,
        m_ini_idx, binar_flag):
    """
    Generate synthetic cluster given by 'model'.
    """
    synth_clust, sigma, extra_pars, isoch_moved, mass_dist, isoch_binar,\
        isoch_compl = synth_cluster.main(
            fundam_params, varIdxs, model, theor_tracks, completeness,
            max_mag_syn, st_dist_mass, R_V, ext_coefs, N_fc, err_pars,
            m_ini_idx, binar_flag, extra_pars_flag)

    return isoch_moved, synth_clust, sigma, extra_pars


def writeFileOut(
    npd, filters, colors, synth_clst, sigma, extra_pars, cl_max_mag, st_mean,
        st_std, st_mean_binar, st_std_binar, prob_binar):
    """
    Save best fit synthetic cluster found, and per star masses to file.
    """

    # Write synthetic cluster file
    e_mags_cols = sigma.T
    binar_idx, ini_mass = extra_pars[0], extra_pars[2]
    # Create IDs identifying binary systems
    binar_ID = []
    for i, bi in enumerate(binar_idx):
        if bi > 1.:
            binar_ID.append('2' + str(i))
        else:
            binar_ID.append('1' + str(i))
    phot_col = [f[1] for f in filters] +\
        ['(' + c[1].replace(',', '-') + ')' for c in colors]
    ephot_col = ['e_' + f[1] for f in filters] +\
        ['e_(' + c[1].replace(',', '-') + ')' for c in colors]
    col_names = ['ID'] + phot_col + ephot_col + ['Mini']
    data = [binar_ID] + [_ for _ in synth_clst.T] +\
        [_ for _ in e_mags_cols.T] + [ini_mass]
    synth_table = Table(data, names=col_names)
    synth_table.meta['comments'] = ["Binary systems ID's begin with a '2'"]
    data_IO.dataSave(synth_table, npd['synth_file_out'], 'w')

    # Write masses to file
    st_ID = np.array(list(zip(*cl_max_mag))[0])
    main_mag = np.array(list(zip(*cl_max_mag))[3]).T[0]
    first_col = np.array(list(zip(*cl_max_mag))[5]).T[0]
    mass_table = Table(
        [st_ID, main_mag, first_col, st_mean, st_std, st_mean_binar,
         st_std_binar, prob_binar],
        names=['ID', 'Mag', 'Col', 'Mass_mu', 'Mass_std', 'Mass_binar_mu',
               'Mass_binar_std', 'P_binar'])
    mass_table.meta['comments'] = ['Subset of stars selected as members']
    data_IO.dataSave(
        mass_table, npd['mass_file_out'], 'w',
        {'Mag': '%12.5f', 'Col': '%12.5f', 'Mass_mu': '%12.5f',
         'Mass_std': '%12.5f', 'Mass_binar_mu': '%12.5f',
         'Mass_binar_std': '%12.5f', 'P_binar': '%12.2f'})

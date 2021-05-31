
import numpy as np
from . import synth_cluster
from . add_errors import getSigmas


def main(clp, pd, td):
    """
    """
    if pd['best_fit_algor'] != 'n':
        # Use the selected solution values for all the parameters.
        model = clp['isoch_fit_params'][pd['D3_sol'] + '_sol']

        # Pack common args.
        clp['syntClustArgs'] = (
            td['fundam_params'], clp['isoch_fit_params']['varIdxs'],
            clp['completeness'], clp['err_lst'],
            clp['max_mag_syn'], td['ext_coefs'], td['binar_flag'],
            td['mean_bin_mr'], td['N_fc'], td['m_ini_idx'],
            td['st_dist_mass'], td['theor_tracks'], td['err_norm_rand'],
            td['binar_probs'], td['ext_unif_rand'])

        # Generate isochrone, synthetic cluster (with uncertainties), and sigma
        # values for the "best" fitted parameters.
        # Model with the "best" fitted parameters.
        model_var = np.array(model)[clp['isoch_fit_params']['varIdxs']]
        # shape: (N_dim, N_stars)
        synth_clust = setSynthClust(model_var, *clp['syntClustArgs'])
        # Get uncertainties
        sigma = []
        for i, popt_mc in enumerate(clp['err_lst']):
            sigma.append(getSigmas(synth_clust[0], popt_mc))

        # Save for plotting and storing.
        clp['synth_cl_phot'], clp['synth_cl_sigma'] = synth_clust, sigma

        # Mask that points to *binary* systems
        if td['binar_flag']:
            clp['binar_idx'] = ~(synth_clust[-1] == -99.)
        else:
            clp['binar_idx'] = np.array([
                False for _ in range(synth_clust.shape[-1])])

    return clp


def setSynthClust(
    model, fundam_params, varIdxs, completeness, err_lst,
    max_mag_syn, ext_coefs, binar_flag, mean_bin_mr, N_fc, m_ini_idx,
    st_dist_mass, theor_tracks, err_norm_rand, binar_probs,
        ext_unif_rand):
    """
    Generate synthetic cluster given by 'model'.
    """
    # This returns the non-reduced, non-transposed array
    transpose_flag = False

    return synth_cluster.main(
        fundam_params, varIdxs, model, completeness, err_lst,
        max_mag_syn, ext_coefs, binar_flag, mean_bin_mr, N_fc, m_ini_idx,
        st_dist_mass, theor_tracks, err_norm_rand, binar_probs,
        ext_unif_rand, transpose_flag)

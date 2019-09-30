
from . import move_isochrone
from . import synth_cluster
from ..best_fit.bf_common import closeSol


def main(
    best_fit_algor, fundam_params, theor_tracks, plot_isoch_data, R_V,
    e_max, isoch_fit_params, err_lst, completeness, max_mag_syn, st_dist_mass,
        ext_coefs, N_fc, cmpl_rnd, err_rnd):
    '''
    # Generate shifted isochrone and synthetic cluster for plotting.
    '''
    if best_fit_algor == 'boot+GA':
        # Use ML fit values for all parameters.
        synth_cl_params = isoch_fit_params['map_sol']
    elif best_fit_algor == 'ptemcee':
        # Use mean fit values for all parameters.
        synth_cl_params = isoch_fit_params['mean_sol']
    # Grid values for (z, a)
    synth_cl_params_grid = closeSol(fundam_params, synth_cl_params, [0, 1])

    # Find indexes for metallicity and age. If indexes are not found due
    # to some difference in the significant figures, use the indexes
    # [0, 0] to prevent the code from halting.
    try:
        m_i = fundam_params[0].index(synth_cl_params_grid[0])
    except Exception:
        m_i = 0
        print("  WARNING: metallicity for the best match synthetic\n"
              "  cluster not found.")
    try:
        a_i = fundam_params[1].index(synth_cl_params_grid[1])
    except Exception:
        a_i = 0
        print("  WARNING: age for the best match synthetic\n"
              "  cluster not found.")

    # Generate shifted best fit isochrone.
    isochrone = plot_isoch_data[m_i][a_i]
    shift_isoch = move_isochrone.main(
        isochrone, synth_cl_params[2], synth_cl_params[3], R_V, ext_coefs,
        N_fc)

    # Generate best fit synthetic cluster.
    isochrone = theor_tracks[m_i][a_i]
    synth_clst = synth_cluster.main(
        e_max, err_lst, completeness, max_mag_syn, st_dist_mass, isochrone,
        R_V, ext_coefs, N_fc, cmpl_rnd, err_rnd, synth_cl_params_grid)

    return shift_isoch, synth_clst


import move_isochrone
import synth_cluster


def main(e_max, fundam_params, theor_tracks, isoch_fit_params, err_lst,
         completeness, max_mag_syn, st_dist_mass, R_V, ext_coefs, N_fc):
    '''
    # Generate shifted isochrone and synthetic cluster for plotting.
    '''
    # Read best fit values for all parameters.
    synth_cl_params = isoch_fit_params[0]
    m, a, e, d, mass, binar_f = synth_cl_params
    # Find indexes for metallicity and age. If indexes are not found due
    # to some difference in the significant figures, use the indexes
    # [0, 0] to prevent the code from halting.
    try:
        m_i, a_i = fundam_params[0].index(m), fundam_params[1].index(a)
    except Exception:
        m_i, a_i = [0, 0]
        print("  WARNING: metallicity and age for best match synthetic\n"
              "  cluster not found.")
    # Generate shifted best fit isochrone.
    isochrone = theor_tracks[m_i][a_i]
    shift_isoch = move_isochrone.main(isochrone, e, d, R_V, ext_coefs, N_fc)
    # Generate best fit synthetic cluster.
    synth_clst = synth_cluster.main(
        e_max, err_lst, completeness, max_mag_syn, st_dist_mass, isochrone,
        R_V, ext_coefs, N_fc, synth_cl_params)

    return shift_isoch, synth_clst

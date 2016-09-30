
import move_isochrone
import synth_cluster


def main(ip_list, isoch_fit_params, err_lst, completeness,
         st_dist_mass, e_max, bin_mass_ratio, cmd_sel):
    '''
    # Generate shifted isochrone and synthetic cluster for plotting.
    '''
    # Get list of stored isochrones and their parameters.
    isoch_list, param_values = ip_list[0], ip_list[1]
    # Read best fit values for all parameters.
    m, a, e, d, mass, binar_f = isoch_fit_params[0]
    # Find indexes for metallicity and age. If indexes are not found due
    # to some difference in the significant figures, use the indexes
    # [0, 0] to prevent the code from halting.
    try:
        m_i, a_i = param_values[0].index(m), param_values[1].index(a)
    except:
        m_i, a_i = [0, 0]
    # Generate shifted best fit isochrone.
    shift_isoch = move_isochrone.main(isoch_list[m_i][a_i][:2], e, d, cmd_sel)
    # Generate best fit synthetic cluster.
    synth_clst = synth_cluster.main(
        err_lst, completeness, st_dist_mass, isoch_list[m_i][a_i],
        [-1., -1., e, d, mass, binar_f], e_max, bin_mass_ratio, cmd_sel)

    return shift_isoch, synth_clst

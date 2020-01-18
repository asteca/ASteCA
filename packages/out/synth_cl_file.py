
import numpy as np
from ..synth_clust import synth_cluster
from ..synth_clust import zaWAverage
from ..synth_clust import move_isochrone
from ..synth_clust import cut_max_mag
from ..synth_clust import mass_distribution
from ..synth_clust import mass_interp
from ..synth_clust import binarity
from ..synth_clust import completeness_rm
from ..synth_clust import add_errors


def main(clp, npd, bf_flag, best_fit_algor, fundam_params, filters, colors,
         theor_tracks, R_V, m_ini_idx, binar_flag, **kwargs):
    """
    Create output data file with stars in the best fit synthetic cluster found
    by the 'Best Fit' function.
    """

    clp['synth_clst'], clp['shift_isoch'] = [], []
    if bf_flag:

        shift_isoch, synth_clst = synth_cl_plot(
            best_fit_algor, fundam_params, clp['isoch_fit_params'],
            theor_tracks, clp['completeness'], clp['max_mag_syn'],
            clp['st_dist_mass'], R_V, clp['ext_coefs'], clp['N_fc'],
            clp['err_pars'], m_ini_idx, binar_flag)  # clp['err_pars_old']

        # If cluster is not empty.
        if synth_clst[0].any():
            # Prepare data.
            mags_cols, e_mags_cols = synth_clst[0].T, synth_clst[1].T
            binar_idx, extra_pars = synth_clst[2][0].T, synth_clst[2][2:].T
            # Prepare header.
            hdr = ['#ID  '] + [f[1] + '   ' for f in filters]
            hdr += ['(' + c[1].replace(',', '-') + ')   ' for c in colors]
            hdr += ['e_' + f[1] + '   ' for f in filters]
            hdr += ['e_(' + c[1].replace(',', '-') + ')   ' for c in colors]
            hdr = ''.join(hdr) + 'Mini\n'
            # int_IMF  Mass  logL  logTe  logg  label  mbolmag\n
            # Save best fit synthetic cluster found to file.
            with open(npd['synth_file_out'], "w") as f_out:
                f_out.write(hdr)
                for i, st in enumerate(mags_cols):
                    if binar_idx[i] > 1.:
                        ID = '2' + str(i)
                    else:
                        ID = '1' + str(i)
                    f_out.write("{:<8}".format(ID))
                    for _ in st:
                        f_out.write("{:<8.4f}".format(_))
                    for _ in e_mags_cols[i]:
                        f_out.write("{:<8.4f}".format(_))
                    f_out.write("{:<8.4f}\n".format(*extra_pars[i]))

            # Save for plotting.
            clp['synth_clst'], clp['shift_isoch'] = synth_clst, shift_isoch

        else:
            print("  ERROR: empty synthetic cluster could not be saved\n"
                  "  to file")

    return clp


def synth_cl_plot(
    best_fit_algor, fundam_params, isoch_fit_params, theor_tracks,
    completeness, max_mag_syn, st_dist_mass, R_V, ext_coefs, N_fc, err_pars,
        m_ini_idx, binar_flag):
    """
    Generate shifted isochrone and synthetic cluster for plotting.
    """

    if best_fit_algor == 'boot+GA':
        # Use ML fit values for all parameters.
        model = isoch_fit_params['map_sol']
    elif best_fit_algor in ('ptemcee', 'emcee'):
        # Use mean fit values for all parameters.
        model = isoch_fit_params['mean_sol']

    # Generate a model with the "best" fitted parameters.
    model_var = np.array(model)[isoch_fit_params['varIdxs']]

    # Generate best fit synthetic cluster.
    model_proper, z_model, a_model, ml, mh, al, ah = synth_cluster.properModel(
        fundam_params, model_var, isoch_fit_params['varIdxs'])
    isochrone = zaWAverage.main(
        theor_tracks, fundam_params, z_model, a_model, ml, mh, al, ah)
    e, d, M_total, bin_frac = model_proper
    isoch_moved = move_isochrone.main(
        isochrone, e, d, R_V, ext_coefs, N_fc, binar_flag, m_ini_idx)
    isoch_cut = cut_max_mag.main(isoch_moved, max_mag_syn)
    synth_clust = np.array([])
    if isoch_cut.any():
        mass_dist = mass_distribution.main(st_dist_mass, M_total)
        isoch_mass = mass_interp.main(isoch_cut, mass_dist, m_ini_idx)
        if isoch_mass.any():
            isoch_binar = binarity.main(isoch_mass, bin_frac, m_ini_idx, N_fc)
            isoch_compl = completeness_rm.main(isoch_binar, completeness)
            if isoch_compl.any():
                # This is what makes this call different from 'synth_cluster'
                synth_clust = add_errors.main(
                    isoch_compl, err_pars, binar_flag, m_ini_idx, True)

    return isoch_moved, synth_clust

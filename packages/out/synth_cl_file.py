
from ..synth_clust import synth_cl_plot


def main(clp, npd, bf_flag, best_fit_algor, fundam_params, filters, colors,
         theor_tracks, R_V, **kwargs):
    '''
    Create output data file with stars in the best fit synthetic cluster found
    by the 'Best Fit' function.
    '''
    clp['synth_clst'], clp['shift_isoch'] = [], []
    if bf_flag:

        # Index of m_ini (theoretical initial mass), stored in the theoretical
        # isochrones.
        m_ini = 2 * clp['N_fc'][0] + 2 * clp['N_fc'][1] + 2
        # Pack synthetic cluster arguments.
        synthcl_args = [
            theor_tracks, clp['em_float'], clp['err_lst'], clp['completeness'],
            clp['max_mag_syn'], clp['st_dist_mass'], R_V, clp['ext_coefs'],
            clp['N_fc'], m_ini, clp['err_rnd']]

        shift_isoch, synth_clst = synth_cl_plot.main(
            best_fit_algor, fundam_params, clp['isoch_fit_params'],
            synthcl_args)

        # If cluster is not empty.
        if synth_clst:
            # Prepare data.
            mags_cols, e_mags_cols = list(zip(*synth_clst[0][0])),\
                list(zip(*synth_clst[0][1]))
            binar_idx = synth_clst[1][0]
            extra_pars = list(zip(*synth_clst[1][2:]))
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

            print("Best fit synthetic cluster saved to file")
            # Save for plotting.
            clp['synth_clst'], clp['shift_isoch'] = synth_clst, shift_isoch

        else:
            print("  ERROR: empty synthetic cluster could not be saved\n"
                  "  to file")

    return clp

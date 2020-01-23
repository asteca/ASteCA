
import numpy as np
from ..inp.get_data import flatten


def main(
    npd, pd, flag_3pk_conver, flag_memb_par, flag_num_memb_low, K_memb_num,
    K_conct_par, cont_index, n_memb, memb_par, n_memb_da, frac_cl_area,
    x_offset, y_offset, kde_cent, clust_rad, e_rad, core_rad, e_core,
        tidal_rad, e_tidal, isoch_fit_params, isoch_fit_errors, **kwargs):
    '''
    Add data obtained to the 'asteca_output.dat' file.
    '''

    # Unpack data.
    out_file_name, write_name = npd['out_file_name'], npd['write_name']
    # Invert flag.
    flag_3pk_no_conver = not flag_3pk_conver

    # Create list containing all the flags.
    flags_list = [flag_3pk_no_conver, flag_num_memb_low, flag_memb_par]

    # Convert True & False flag values to 1 and 0 respectively.
    int_flags = [1 if flg else 0 for flg in flags_list]

    # Sum all flags to obtain the FC value (flags count) and append to the
    # end of the list. Do not count the manual flags, hence the [2:].
    int_flags.append(sum(int_flags))

    # Round structure parameters.
    frmt = "{:.6f}" if pd['coords'] == 'deg' else "{:.0f}"
    if pd['coords'] == 'deg' and pd['project']:
        x_cent = (kde_cent[0] / np.cos(np.deg2rad(kde_cent[1] + y_offset))) +\
            x_offset
    else:
        x_cent = kde_cent[0]
    cr_r = [frmt.format(_) for _ in
            [x_cent, kde_cent[1] + y_offset, clust_rad, core_rad, tidal_rad]]
    cr_e = [frmt.format(_) for _ in [e_rad[0], e_rad[1], e_core, e_tidal]]

    # Interwine these lists.
    cre_r = cr_r[:2] + [
        cr_r[2], cr_e[0], cr_e[1], cr_r[3], cr_e[2], cr_r[4], cr_e[3]]

    # Cluster parameters and errors.
    cpe_r = [
        item for t in zip(
            isoch_fit_params['mean_sol'], isoch_fit_params['map_sol'],
            isoch_fit_params['median_sol'], isoch_fit_params['mode_sol'],
            list(zip(*isoch_fit_errors))[0], list(zip(*isoch_fit_errors))[1],
            list(zip(*isoch_fit_errors))[2], isoch_fit_params['param_r2'])
        for item in t]

    # Store all parameter values in list.
    # Using main magnitude only
    line = [write_name, cre_r, K_conct_par, cont_index, K_memb_num,
            n_memb, n_memb_da, memb_par, frac_cl_area, cpe_r,
            isoch_fit_params['N_total']]
    # Flatten list.
    line_f = list(flatten(line))

    # TDOD not sure if this is really an improvement
    # with open(out_file_name, mode='a') as f:
    #     # Some platforms don't automatically seek to end when files opened
    #     # in append mode
    #     f.seek(0, os.SEEK_END)
    #     t2 = Table(zip(*[line_f + int_flags]))
    #     t2.write(f, format='ascii.no_header', formats={'col0': '%-16s'})

    frmts = (
        '{:<30} {:>10} {:>10} {:>10} {:>10} {:>10} ' +
        '{:>10} {:>10} {:>10} {:>10} {:>8.2f} {:>7.2f} {:>10.0f} ' +
        '{:>10.0f} {:>10.0f} {:>9.2f} {:>7.2f} ' +
        '{:>10.5} {:>10.5} {:>10.5} {:>10.5} {:>10.5} {:>10.5} {:>10.5} {:>6.2} ' +
        '{:>10.4} {:>10.4} {:>10.4} {:>10.4} {:>10.4} {:>10.4} {:>10.4} {:>6.2} ' +
        '{:>10.3} {:>10.3} {:>10.3} {:>10.3} {:>10.3} {:>10.3} {:>10.3} {:>6.2} ' +
        '{:>10.5} {:>10.5} {:>10.5} {:>10.5} {:>10.5} {:>10.5} {:>10.5} {:>6.2} ' +
        '{:>10.5} {:>10.5} {:>10.5} {:>10.5} {:>10.5} {:>10.5} {:>10.5} {:>6.2} ' +
        '{:>10.2} {:>10.2} {:>10.2} {:>10.2} {:>10.2} {:>10.2} {:>10.2} {:>6.2} ' +
        '{:>10.2E}')
    # Write values to file.
    with open(out_file_name, "a") as f_out:
        f_out.write(frmts.format(*line_f))
        # Flags.
        #               f1    f2    f3    FC
        f_out.write('''{:>8} {:>2} {:>2} {:>3}'''.format(
            *int_flags))
        f_out.write('\n')

    print("Analysis results added to output file")

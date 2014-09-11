"""
@author: gabriel
"""

from error_round import round_sig_fig as rsf
from compiler.ast import flatten


def add_data_output(out_file_name, sub_dir, output_dir, clust_name,
    center_params, radius_params, kp_params, cont_index, n_c, prob_cl_kde,
    integr_return, err_flags, flag_num_memb_low, bf_return):
    '''
    Add data obtained to the 'data_output.dat' file.
    '''

    # Unpack data.
    center_cl, e_cent = center_params[1][0], center_params[2]
    flag_center_med, flag_center_std, flag_center_manual = center_params[-3:]
    clust_rad, e_rad = radius_params[:2]
    flag_delta_total, flag_not_stable, flag_delta, flag_radius_manual = \
    radius_params[-4:]
    rc, e_rc, rt, e_rt, n_c_k = kp_params[:5]
    err_all_fallback, err_max_fallback = err_flags
    # Unpack KP flags.
    flag_2pk_conver, flag_3pk_conver = kp_params[-2:]
    # Invert flag.
    flag_3pk_no_conver = not flag_3pk_conver

    # Get parameter from list.
    if integr_return:
        integ_mag, integ_col = integr_return[2], integr_return[5]
    else:
        integ_mag, integ_col = 1., -98.

    # Create list containing all the flags.
    flags_list = [flag_center_manual, flag_radius_manual, flag_center_med,
        flag_center_std, flag_delta_total, flag_not_stable, flag_delta,
        flag_3pk_no_conver, err_all_fallback, err_max_fallback,
        flag_num_memb_low]

    # Converty True & False flag values to 1 and 0 respectively.
    int_flags = [1 if flg else 0 for flg in flags_list]

    # Sum all flags to obtain the FC value (flags count) and append to the
    # end of the list. Do not count the manual flags, hence the [2:].
    int_flags.append(sum(int_flags[2:]))

    # Round structure parameters.
    # If 3-P King profile converged.
    if flag_3pk_no_conver is False:
        rtt, e_rtt = rt, e_rt
    elif flag_2pk_conver is True:
        rtt, e_rtt = -1., -1.
    cr_r, cr_e = rsf([center_cl[0], center_cl[1], clust_rad, rc, rtt],
        [e_cent[0], e_cent[1], e_rad, e_rc, e_rtt])
    # Interwine these lists.
    cre_r = [item for t in zip(cr_r, cr_e) for item in t]

    # Unpack cluster param values and their errors.
    m, a, e, d = bf_return[0][0]
    e_m, e_a, e_e, e_d = bf_return[1]
    # Round cluster parameters.
    cp_r, cp_e = rsf([m, a, e, d], [e_m, e_a, e_e, e_d])
    # Interwine these lists.
    cpe_r = [item for t in zip(cp_r, cp_e) for item in t]

    # Store all parameter values in list.
    line = [str(sub_dir) + '/' + str(clust_name), cre_r,
        cont_index, n_c, n_c_k, prob_cl_kde, (integ_col - integ_mag),
        cpe_r]
    # Flatten list.
    line_f = flatten(line)

    # Write values to file.
    with open(output_dir + out_file_name, "a") as f_out:
        f_out.write('''{:<16} {:>8} {:>8} {:>8} {:>8} {:>8} \
{:>8} {:>8} {:>8} {:>8} {:>8} \
{:>8.2f} {:>8.0f} {:>8.0f} {:>8.0f} {:>8.2f} \
{:>8} {:>8} {:>8} {:>8} \
{:>8} {:>8} {:>8} {:>8}'''.format(*line_f))
        # Flags.
        f_out.write('''{:>8} {:>2} {:>3} {:>2} {:>2} {:>2} {:>2} {:>2} {:>2} \
{:>2} {:>2} {:>3}'''.format(*int_flags))
        f_out.write('\n')

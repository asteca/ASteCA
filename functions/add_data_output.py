"""
@author: gabriel
"""


def add_data_output(out_file_name, sub_dir, output_dir, clust_name,
    center_params, radius_params, kp_params, cont_index, n_c, prob_cl_kde,
    ccc, integr_return, rjct_errors_fit, flag_num_memb_low, bf_return):
    '''
    Add data obtained to the 'data_output.dat' file.
    '''

    # Unpack data.
    center_cl = center_params[5][0]
    flag_center_med, flag_center_std, flag_center_manual = center_params[-3:]
    clust_rad = radius_params[0]
    flag_delta_total, flag_not_stable, flag_delta, flag_radius_manual = \
    radius_params[-4:]
    rc, e_rc, rt, e_rt, n_c_k = kp_params[:5]
    # Invert flag.
    flag_3pk_no_conver = not kp_params[-1]

    # Get parameter from list.
    integ_mag, integ_col = integr_return[2], integr_return[5]

    # Create list containing all the flags.
    flags_list = [flag_center_manual, flag_radius_manual, rjct_errors_fit,
                  flag_center_med, flag_center_std, flag_delta_total,
                  flag_not_stable, flag_delta,
                  flag_3pk_no_conver, flag_num_memb_low]

    # Converty True & False values to 1 and 0 respectively.
    int_flags = [1 if flg else 0 for flg in flags_list]

    # Sum all flags to obtain the FC value (flags count) and append to the
    # end of the list. Do not count the manual flags, hence the [3:].
    int_flags.append(sum(int_flags[3:]))

    isoch_fit_params, isoch_fit_errors = bf_return[0], bf_return[1]
    m, a, e, d = isoch_fit_params[0]
    e_m, e_a, e_e, e_d = isoch_fit_errors
    # str(sub_dir)+'_'+str(clust_name)
    line = [str(sub_dir) + '/' + str(clust_name),
        str('%.1f' % round(center_cl[0], 1)),
        str('%.1f' % round(center_cl[1], 1)),
        str('%.1f' % round(clust_rad, 1)), str('%.1f' % round(rc, 1)),
        str('%.1f' % round(e_rc, 1)),
        str('%.1f' % round(rt, 1)),
        str('%.1f' % round(e_rt, 1)),
        str(round(cont_index, 2)),
        str(int(n_c)), str(n_c_k), str('%0.2f' % prob_cl_kde),
        str('%0.2f' % ccc), str('%0.2f' % (integ_col - integ_mag)),
        str('%0.4f' % m), str('%0.4f' % e_m), str('%0.2f' % a),
        str('%0.2f' % e_a), str('%0.2f' % e), str('%0.2f' % e_e),
        str('%0.2f' % d), str('%0.2f' % e_d)]

    # "a" opens the file for appending
    with open(output_dir + out_file_name, "a") as f_out:
        f_out.write('{:<16} {:>7} {:>7} {:>8} {:>7} {:>8} {:>7} {:>8} {:>8} \
{:>4} {:>6} {:>7} {:>5} {:>7} {:>7} {:>7} {:>5} {:>5} {:>7} {:>5} {:>6} \
{:>5}'.format(*line))
        # Flags.
        f_out.write('{:>4} {:>2} {:>2} {:>3} {:>2} {:>2} {:>2} {:>2} {:>2} \
{:>2} {:>3}'.format(*int_flags))
        f_out.write('\n')

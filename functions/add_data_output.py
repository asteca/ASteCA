"""
@author: gabriel
"""


def add_data_output(sub_dir, output_dir, clust_name, center_cl, clust_rad,
                    k_prof, n_c_k, flag_king_no_conver, cont_index, n_c,
                    prob_cl_kde, ccc, integr_return,
                    flag_center, flag_center_manual,
                    flag_radius_manual, rjct_errors_fit, flag_bin_count,
                    radius_params, flag_num_memb_low, bf_return):
    '''Add data obtained to the 'data_output.dat' file.'''

    # Parameters from get_radius function.
    flag_delta_total, flag_not_stable, flag_rad_500, flag_delta,\
    flag_delta_points = radius_params

    # Get parameter from list.
    integ_mag = integr_return[2]

    # Create list containing all the flags.
    flags_list = [flag_center_manual, flag_radius_manual, rjct_errors_fit,
                  flag_center, flag_bin_count, flag_delta_total,
                  flag_not_stable, flag_rad_500, flag_delta, flag_delta_points,
                  flag_king_no_conver, flag_num_memb_low]

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
        str('%0.f' % round(center_cl[0])), str('%0.f' % round(center_cl[1])),
        str('%0.f' % round(clust_rad)), str('%0.f' % round(k_prof[0])),
        str('%0.f' % round(k_prof[1])), str(round(cont_index, 2)),
        str(int(n_c)), str(n_c_k), str('%0.2f' % prob_cl_kde),
        str('%0.2f' % ccc), str('%0.2f' % integ_mag),
        str('%0.4f' % m), str('%0.4f' % e_m), str('%0.2f' % a),
        str('%0.2f' % e_a), str('%0.2f' % e), str('%0.2f' % e_e),
        str('%0.2f' % d), str('%0.2f' % e_d)]

    # "a" opens the file for appending
    with open(output_dir + 'data_output.dat', "a") as f_out:
        f_out.write('{:<16} {:>7} {:>7} {:>8} {:>7} {:>7} {:>8} {:>4} {:>6} \
{:>7} {:>5} {:>7} {:>7} {:>7} {:>5} {:>5} {:>7} {:>5} {:>6} \
{:>5}'.format(*line))
        # Flags.
        f_out.write('{:>4} {:>2} {:>2} {:>3} {:>2} {:>2} {:>2} {:>2} {:>2} \
{:>2} {:>2} {:>2} {:>3}'.format(*int_flags))
        f_out.write('\n')

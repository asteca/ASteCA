"""
@author: gabriel
"""

''' Add data obtained to the 'data_output' file.'''

def add_data_output(sub_dir, output_dir, clust_name, center_cl, clust_rad,
                    k_prof, n_c_k, flag_king_no_conver, cont_index, n_c, p_value,
                    flag_center, flag_std_dev, flag_center_manual, 
                    flag_radius_manual, flag_errors_manual, flag_bin_count, 
                    flag_delta_total, flag_not_stable, flag_rad_500, flag_delta,
                    flag_delta_points, flag_num_memb_low, flag_no_memb):

    # Create list containing all the flags.  
    flags_list = [flag_center_manual, flag_radius_manual, flag_errors_manual,
                  flag_center, flag_std_dev, flag_bin_count, flag_delta_total,
                  flag_not_stable, flag_rad_500, flag_delta, flag_delta_points,
                  flag_king_no_conver, flag_num_memb_low, flag_no_memb]
    
    # Converty True & False values to 1 and 0 respectively.
    int_flags = [1 if flg else 0 for flg in flags_list]

    # Sum all flags to obtain CI value (cluster index) and append to the end of
    # the list. Do not count the manual flags, hence the [3:].
    int_flags.append(sum(int_flags[3:]))

    # str(sub_dir)+'_'+str(clust_name)
    line = [str(sub_dir)+'/'+str(clust_name), str('%0.f' % round(center_cl[0])),
            str('%0.f' % round(center_cl[1])), str('%0.f' % round(clust_rad[0])),\
            str('%0.f' % round(k_prof[0])), str('%0.f' % round(k_prof[1])), \
            str(round(cont_index, 2)), str(int(n_c)), str(n_c_k), \
            str('%0.2f' % p_value)]
    
    # "a" opens the file for appending
    with open(output_dir+'data_output', "a") as f_out:
        f_out.write('{:<16} {:>7} {:>7} {:>8} {:>7} {:>7} {:>8} {:>4} {:>6} {:>7}'.format(*line))
        f_out.write('{:>4} {:>2} {:>2}  {:>2} {:>2} {:>2} {:>2} {:>2} {:>2} {:>2}\
 {:>2} {:>2} {:>3} {:>3} {:>3}'.format(*int_flags))
        f_out.write('\n')

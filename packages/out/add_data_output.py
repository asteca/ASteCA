
import collections


def flatten(l):
    '''
    Flatten list.
    '''
    for el in l:
        if isinstance(el, collections.Iterable) and not \
                isinstance(el, basestring):
            for sub in flatten(el):
                yield sub
        else:
            yield el


def main(npd, pd, flag_center_std, flag_center_manual,
         flag_delta_total, flag_not_stable, flag_delta, flag_radius_manual,
         flag_2pk_conver, flag_3pk_conver, flag_memb_par, flag_num_memb_low,
         err_flags, K_memb_num, K_conct_par, cont_index, n_memb, memb_par,
         n_memb_da, frac_cl_area, pval_test_params, integ_mag, kde_cent,
         clust_rad, e_rad, core_rad, e_core, tidal_rad, e_tidal,
         fit_params_r, fit_errors_r, **kwargs):
    '''
    Add data obtained to the 'data_output.dat' file.
    '''

    # Unpack data.
    out_file_name, write_name = npd['out_file_name'], npd['write_name']
    # Invert flag.
    flag_3pk_no_conver = not flag_3pk_conver
    err_all_fallback, err_max_fallback = err_flags

    # Create list containing all the flags.
    flags_list = [flag_center_manual, flag_radius_manual,
                  flag_center_std, flag_delta_total, flag_not_stable,
                  flag_delta, flag_3pk_no_conver, err_all_fallback,
                  err_max_fallback, flag_num_memb_low, flag_memb_par]

    # Convert True & False flag values to 1 and 0 respectively.
    int_flags = [1 if flg else 0 for flg in flags_list]

    # Sum all flags to obtain the FC value (flags count) and append to the
    # end of the list. Do not count the manual flags, hence the [2:].
    int_flags.append(sum(int_flags[2:]))

    # Round structure parameters.
    cr_r = ["{:.0f}".format(_) for _ in
            [kde_cent[0], kde_cent[1], clust_rad, core_rad, tidal_rad]]
    cr_e = ["{:.0f}".format(_) for _ in [e_rad, e_core, e_tidal]]
    # Interwine these lists.
    cre_r = cr_r[:2] + [item for t in zip(cr_r[2:], cr_e) for item in t]
    # Rounded cluster parameters and errors.
    cpe_r = [item for t in zip(fit_params_r, fit_errors_r) for item in t]

    # Store all parameter values in list.
    # TODO using main magnitude only
    line = [write_name, cre_r, K_conct_par, cont_index, K_memb_num,
            n_memb, n_memb_da, memb_par, frac_cl_area, pval_test_params[0],
            integ_mag[0], cpe_r]
    # Flatten list.
    line_f = list(flatten(line))

    # Write values to file.
    with open(out_file_name, "a") as f_out:
        f_out.write('''{:<16} {:>8} {:>8} {:>8} \
{:>8} {:>8} {:>8} {:>8} {:>8} {:>8.2f} {:>7.2f} {:>10.0f} \
{:>10.0f} {:>10.0f} {:>9.2f} {:>7.2f} {:>8.2f} {:>8.2f} \
{:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} \
{:>8} {:>8} {:>8} {:>8}'''.format(*line_f))
        # Flags.
        f_out.write('''{:>8} {:>2} {:>3} {:>2} {:>2} {:>2} {:>2} {:>2} \
{:>2} {:>2} {:>2} {:>3}'''.format(*int_flags))
        f_out.write('\n')

    print("Analysis results added to output file.")


import collections
from ..errors.error_round import round_sig_fig


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


def main(npd, pd, clust_cent, e_cent, flag_center_med, flag_center_std,
         flag_center_manual, clust_rad, e_rad, flag_delta_total,
         flag_not_stable, flag_delta, flag_radius_manual, core_rad,
         e_core, tidal_rad, e_tidal, K_memb_num, K_conct_par, flag_2pk_conver,
         flag_3pk_conver, cont_index, n_memb, memb_par, n_memb_da,
         flag_memb_par, frac_cl_area, pval_test_params, integr_return,
         err_flags, flag_num_memb_low, isoch_fit_params, isoch_fit_errors,
         **kwargs):
    '''
    Add data obtained to the 'data_output.dat' file.
    '''

    # Unpack data.
    out_file_name, write_name = npd['out_file_name'], npd['write_name']
    # Invert flag.
    flag_3pk_no_conver = not flag_3pk_conver
    err_all_fallback, err_max_fallback = err_flags

    # Construct integrated color.
    m_ord = pd['axes_params'][2]
    if integr_return:
        integ_mag1, integ_mag2 = integr_return[2], integr_return[5]
        sig = 1. if m_ord == 21 else -1.
        integ_col = sig * (integ_mag2 - integ_mag1)
    else:
        integ_col = -99.

    # Create list containing all the flags.
    flags_list = [flag_center_manual, flag_radius_manual, flag_center_med,
                  flag_center_std, flag_delta_total, flag_not_stable,
                  flag_delta, flag_3pk_no_conver, err_all_fallback,
                  err_max_fallback, flag_num_memb_low, flag_memb_par]

    # Convert True & False flag values to 1 and 0 respectively.
    int_flags = [1 if flg else 0 for flg in flags_list]

    # Sum all flags to obtain the FC value (flags count) and append to the
    # end of the list. Do not count the manual flags, hence the [2:].
    int_flags.append(sum(int_flags[2:]))

    # Round structure parameters.
    cr_r, cr_e = round_sig_fig(
        [clust_cent[0], clust_cent[1], clust_rad, core_rad, tidal_rad],
        [e_cent[0], e_cent[1], e_rad, e_core, e_tidal])
    # Interwine these lists.
    cre_r = [item for t in zip(cr_r, cr_e) for item in t]

    # Round cluster parameters.
    # See if bootstrap process was applied.
    N_b = pd['bf_params'][-1]
    if N_b >= 2:
        # Round cluster params using the g format.
        cp_r = ['{:g}'.format(_) for _ in isoch_fit_params[0]]
        cp_e = isoch_fit_errors
    else:
        # Round cluster parameters.
        cp_r, cp_e = round_sig_fig(isoch_fit_params[0], isoch_fit_errors)
    # Interwine these lists.
    cpe_r = [item for t in zip(cp_r, cp_e) for item in t]

    # Store all parameter values in list.
    line = [write_name, cre_r, K_conct_par, cont_index, K_memb_num,
            n_memb, n_memb_da, memb_par, frac_cl_area, pval_test_params[0],
            integ_col, cpe_r]
    # Flatten list.
    line_f = list(flatten(line))

    # Write values to file.
    with open(out_file_name, "a") as f_out:
        f_out.write('''{:<16} {:>8} {:>8} {:>8} {:>8} {:>8} \
{:>8} {:>8} {:>8} {:>8} {:>8} {:>8.2f} {:>7.2f} {:>10.0f} \
{:>10.0f} {:>10.0f} {:>9.2f} {:>7.2f} {:>8.2f} {:>8.2f} \
{:>8} {:>8} {:>8} {:>8} \
{:>8} {:>8} {:>8} {:>8} \
{:>8} {:>8} {:>8} {:>8}'''.format(*line_f))
        # Flags.
        f_out.write('''{:>8} {:>2} {:>3} {:>2} {:>2} {:>2} {:>2} {:>2} {:>2} \
{:>2} {:>2} {:>3} {:>3}'''.format(*int_flags))
        f_out.write('\n')

    print("Cluster's parameters and flags added to output file.")

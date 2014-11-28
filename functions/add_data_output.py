"""
@author: gabriel
"""

import numpy as np
from error_round import round_sig_fig as rsf
import collections


def flatten(l):
    '''
    Flatten list.
    '''
    for el in l:
        if isinstance(el, collections.Iterable) and not isinstance(el,
            basestring):
            for sub in flatten(el):
                yield sub
        else:
            yield el


def add_data_output(out_file_name, write_name,
    center_params, radius_params, kp_params, cont_index, n_c, prob_cl_kde,
    integr_return, axes_params, err_flags, flag_num_memb_low, bf_return):
    '''
    Add data obtained to the 'data_output.dat' file.
    '''

    # Unpack data.
    center_cl, e_cent = center_params[1], center_params[2]
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

    # Construct integrated color.
    m_ord = axes_params[2]
    if integr_return:
        integ_mag1, integ_mag2 = integr_return[2], integr_return[5]
        sig = 1. if m_ord == 21 else -1.
        integ_col = sig * (integ_mag2 - integ_mag1)
    else:
        integ_col = -99.

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
    cr_r, cr_e = rsf([center_cl[0], center_cl[1], clust_rad, rc, rt],
        [e_cent[0], e_cent[1], e_rad, e_rc, e_rt])
    # Interwine these lists.
    cre_r = [item for t in zip(cr_r, cr_e) for item in t]

    # Round cluster parameters.
    # See if bootstrap process was applied.
    cp_e = bf_return[1]
    if np.array([_ == -1. for _ in cp_e]).all():
        # Round cluster params using the g format.
        cp_r = ['{:g}'.format(_) for _ in bf_return[0][0]]
    else:
        # Round cluster parameters.
        cp_r, cp_e = rsf(bf_return[0][0], bf_return[1])
    # Interwine these lists.
    cpe_r = [item for t in zip(cp_r, cp_e) for item in t]

    # Store all parameter values in list.
    line = [write_name, cre_r,
        cont_index, n_c, n_c_k, prob_cl_kde, integ_col, cpe_r]
    # Flatten list.
    line_f = list(flatten(line))

    # Write values to file.
    with open(out_file_name, "a") as f_out:
        f_out.write('''{:<16} {:>8} {:>8} {:>8} {:>8} {:>8} \
{:>8} {:>8} {:>8} {:>8} {:>8} \
{:>8.2f} {:>8.0f} {:>8.0f} {:>8.2f} {:>8.2f} \
{:>8} {:>8} {:>8} {:>8} \
{:>8} {:>8} {:>8} {:>8}'''.format(*line_f))
        # Flags.
        f_out.write('''{:>8} {:>2} {:>3} {:>2} {:>2} {:>2} {:>2} {:>2} {:>2} \
{:>2} {:>2} {:>3}'''.format(*int_flags))
        f_out.write('\n')

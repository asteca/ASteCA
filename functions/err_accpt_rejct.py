"""
@author: gabriel
"""

import sys
#import matplotlib.pyplot as plt
#from functions.display_errors import disp_errors as d_e
from functions.err_accpt_rejct_lowexp import err_a_r_lowexp as e_a_r_le
from functions.err_accpt_rejct_eyefit import err_a_r_eyefit as e_a_r_ef
from functions.err_accpt_rejct_max import err_a_r_m as e_a_r_m


def err_sel_stars(acpt_indx, rjct_indx, phot_data):
    '''
    Go through the lists of indexes to select accepted/rejected stars.
    '''
    # unpack data.
    id_star, x_data, y_data, mag, e_mag, col1, e_col1 = phot_data
    # Initialize empty lists.
    acpt_stars, rjct_stars = [], []

    for st_ind in acpt_indx:
        acpt_stars.append([id_star[st_ind], x_data[st_ind], y_data[st_ind],
            mag[st_ind], e_mag[st_ind], col1[st_ind], e_col1[st_ind]])

    for st_ind in rjct_indx:
        rjct_stars.append([id_star[st_ind], x_data[st_ind], y_data[st_ind],
            mag[st_ind], e_mag[st_ind], col1[st_ind], e_col1[st_ind]])

    return acpt_stars, rjct_stars


def err_accpt_rejct(phot_data, axes_params, er_params, mode, semi_return):
    """
    Accept and reject stars in and out of the cluster's boundaries according to
    a given criteria based on their photometric errors.
    """

    # Unpack data.
    mag, e_mag, e_col1 = phot_data[3], phot_data[4], phot_data[6]
    er_mode, e_max, be = er_params[:3]

    # Get value of brightest and dimmest stars.
    min_mag, max_mag = min(mag), max(mag)
    # Define max limit for the box that holds the brightest stars.
    bright_end = (min_mag + be)
    # Create a segmented list in magnitude.
    # Magnitude range.
    delta_mag = max_mag - bright_end
    # Width of the intervals in magnitude.
    interv_mag = 0.5
    # Number of intervals.
    n_interv = int(delta_mag / interv_mag)
    # Define list of points spanning the magnitude range starting from the
    # bright end.
    mag_value = [bright_end + interv_mag * (q + 0.5) for q in
    range(n_interv - 1)]

    # Pack params to pass.
    err_pck = [er_params, bright_end, n_interv, interv_mag, mag_value]

    # Check selected mode.
    if mode in ('auto', 'semi'):
        # If 'semi' is set, check for the flag that indicates which method to
        # use and override the one in the input params file.
        if mode == 'semi':
            # Unpack semi flag
            err_flag_semi = semi_return[4]
            if err_flag_semi != 0:
                print 'Semi: using method selected %d.' % int(err_flag_semi)
                # Set error mode to use.
                er_mode == float(err_flag_semi)

        # Flag indicates that the function had to fall back to the
        # 'e_max'-based rejection method since the selected one failed.
        err_max_fallback = False

        # Check which error rejection algorithm was selected in the input
        # file.
        if er_mode == 'emax':
            # Call function to reject stars with errors > e_max.
            acpt_indx, rjct_indx, err_plot = e_a_r_m(e_mag, e_col1, err_pck)

        elif er_mode in ('lowexp', 'eyefit'):
            try:
                if er_mode == 'lowexp':
                    # Call N sigma exp function.
                    acpt_indx, rjct_indx, err_plot = e_a_r_le(mag, e_mag,
                        e_col1, err_pck)

                elif er_mode == 'eyefit':
                    # Call 'eyefit' function.
                    acpt_indx, rjct_indx, err_plot = e_a_r_ef(mag, e_mag,
                        e_col1, err_pck)
            except RuntimeError:
                print ('  WARNING: function could not be fitted. Falling back\n'
                '  to e_max function.')
                # Call function to reject stars with errors > e_max.
                acpt_indx, rjct_indx, err_plot = e_a_r_m(e_mag, e_col1, err_pck)
                err_max_fallback = True

        # Flag indicates that no error rejection was possible.
        err_all_fallback = False
        # If list of accepted stars is empty, fall back to e_max limit.
        if not acpt_indx and er_mode != 'emax':
            print '  WARNING: No stars accepted based on their errors.'
            print '  Using all stars with errors < %0.2f.' % e_max
            # Call function to reject stars with errors > e_max.
            acpt_indx, rjct_indx, err_plot = e_a_r_m(e_mag, e_col1, err_pck)
            err_max_fallback = True
        # If the method used was e_max, use all stars.
        elif not acpt_indx and er_mode == 'emax':
            print '  WARNING: No stars accepted based on their errors.'
            print '  Using all stars.'
            # Store all indexes.
            acpt_indx, rjct_indx = [i for i in range(len(mag))], \
            [i for i in range(len(mag))]
            err_all_fallback = True

    # If 'manual' mode is set, display errors distributions and ask the user
    # to accept it or else use all stars except those with errors > e_max in
    # either the magnitude or the color.
    elif mode == 'manual':
        print 'Plot error distributions.'
        ## Display automatic errors rejection.
        #d_e(phot_data[3], popt_mag, popt_col1, acpt_stars, rjct_stars,
            #err_plot, er_params, axes_params)
        #plt.show()
        ## Ask if keep or reject.
        #wrong_answer = True
        #while wrong_answer:
            #answer_rad = raw_input('Accept fit for errors (otherwise use \
#all stars with photom errors < %0.2f)? (y/n) ' % e_max)
            #if answer_rad == 'y':
                #print 'Fit accepted.'
                #wrong_answer = False
            #elif answer_rad == 'n':
                #print 'Using *all* stars with errors < %0.2f.' % e_max
                ## Call function to reject stars w errors > e_max.
                #popt_mag, popt_col1, acpt_stars, rjct_stars = \
                #e_a_r_m(phot_data, er_params)
                #emax_only = True
                #wrong_answer = False
            #else:
                #print 'Wrong input. Try again.\n'

    if len(acpt_indx) == 0:
        # Halt code if no stars are accepted.
        sys.exit('FATAL: no stars left after error rejection.'
        ' Halting code.')
    else:
        # Call function to store stars according to the returned indexes.
        acpt_stars, rjct_stars = err_sel_stars(acpt_indx, rjct_indx, phot_data)

    print 'Stars accepted/rejected based on their errors.'

    err_flags = [err_all_fallback, err_max_fallback]
    return acpt_stars, rjct_stars, err_plot, err_flags, err_pck
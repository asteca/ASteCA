"""
@author: gabriel
"""

from functions.err_medians import err_med
from functions.exp_function import exp_func
from scipy.optimize import curve_fit
import numpy as np


def exp_fit(mag, bright_end, mag_value, err_val):
    '''
    Fit exponential envelope.
    '''
    popt, pcov = curve_fit(exp_func, mag_value, err_val)
    # Fit polynomial of grade 3 envelope.
    pol = np.polyfit(mag_value, err_val, 3)
    # Find point where curves intersect.
    # Initialize value in case no intersect value is found.
    intersec = 0.
    mag_x = np.linspace(bright_end, max(mag), 100)
    for x_val in mag_x:
        if np.polyval(pol, (x_val)) > exp_func(x_val, *popt):
            intersec = x_val
            break

    return popt, pol, intersec


def fit_curves(mag, mag_value, bright_end, e_mag_value, e_col_value):
    '''
    Get best fitting curves for the median+sigma values obtained for errors
    in magnitude and color.
    '''
    # Process magnitudes.
    popt_mag, pol_mag, intersec_mag = [], [], []
    for err_val in e_mag_value:
        # Fit curves for errors in magnitude.
        pop, pol, inter = exp_fit(mag, bright_end, mag_value, err_val)
        intersec_mag.append(inter)
        popt_mag.append(pop)
        pol_mag.append(pol)

    # Process colors.
    popt_col, pol_col, intersec_col = [], [], []
    for err_val in e_col_value:
        # Fit curves for errors in color.
        pop, pol, inter = exp_fit(mag, bright_end, mag_value, err_val)
        intersec_col.append(inter)
        popt_col.append(pop)
        pol_col.append(pol)

    return intersec_mag, intersec_col, popt_mag, pol_mag, popt_col, pol_col


def separate_stars(mag, e_mag, e_col, e_max, bright_end, be_e,
    intersec_mag, intersec_col, popt_mag, pol_mag, popt_col, pol_col):
    '''
    Use the curves obtained above to accept or reject stars in the
    magnitude range beyond the (brightest star + be) limit.
    '''

    # Initialize empty lists.
    acpt_indx, rjct_indx = [], []

    # Create list of combined errors for all stars.
    em_z, ec_z = zip(*e_mag), zip(*e_col)
    e_mc = [em_z[i] + ec_z[i] for i in range(len(em_z))]

    # Iterate through all stars and accept or reject those beyond
    # the (brightest star + be mag) limit according to the curve
    # obtained for the errors in magnitudes and colors.
    for st_ind, st_mag in enumerate(mag):

        # Reject stars with at least one error >= e_max.
        if any(e >= e_max for e in e_mc[st_ind]):
            rjct_indx.append(st_ind)
        else:
            # For stars brighter than the bright end.
            if mag[st_ind] <= bright_end:
                # For values in this range accept all stars with all errors
                # < be_e.
                if all(e < be_e for e in e_mc[st_ind]):
                    # Accept star.
                    acpt_indx.append(st_ind)
                else:
                    # Reject star.
                    rjct_indx.append(st_ind)

            else:
            # For the reminder of stars, we check to see if they are located
            # above or below the upper envelope for all errors. If they are
            # above in any one, we reject them, otherwise we accept them.

                # Check if star is located to the left or right of the intersect
                # value for each error and compare with the corresponding curve.
                mag_rjct = False
                for j, e_m in enumerate(e_mag):
                    if mag[st_ind] <= intersec_mag[j]:
                        # Compare with exponential.
                        if e_m[st_ind] > exp_func(mag[st_ind], *popt_mag[j]):
                            # Reject star.
                            mag_rjct = True
                    else:
                        # Compare with polynomial.
                        if e_m[st_ind] > np.polyval(pol_mag[j], mag[st_ind]):
                            # Reject star.
                            mag_rjct = True

                col_rjct = False
                for j, e_c in enumerate(e_col):
                    if mag[st_ind] <= intersec_col[j]:
                        # Compare with exponential.
                        if e_c[st_ind] > exp_func(mag[st_ind], *popt_col[j]):
                            # Reject star.
                            col_rjct = True
                    else:
                        # Compare with polynomial.
                        if e_c[st_ind] > np.polyval(pol_col[j], mag[st_ind]):
                            # Reject star.
                            col_rjct = True

                if mag_rjct or col_rjct:
                    # Reject star.
                    rjct_indx.append(st_ind)
                else:
                    # Accept star.
                    acpt_indx.append(st_ind)

    return acpt_indx, rjct_indx


def divide(mag_value, intersec_mag, intersec_col, diag_axis):
    '''
    Divide magnitude interval in two, the first fitted with the
    exponential and the second with the polynomial. Do this for the errors in
    magnitude and in color.
    '''

    # Separate mag_values between those to the left and to the right of this
    # intersect value. The exponential will be the first part of the upper
    # envelope and the polynomial will be the second.
    top_val_left, top_val_right = [], []
    bot_val_left, bot_val_right = [], []

    # Identify the mag/colors to be plotted.
    if diag_axis[0] == 0:
        # A magnitude is used for the top error plot.
        m_ind = diag_axis[1]
        top_intersec = intersec_mag[m_ind]

        if diag_axis[2] == 0:
            # A magnitude is used for the bottom error plot.
            m_ind = diag_axis[3]
            bot_intersec = intersec_mag[m_ind]
        else:
            # A color is used for the bottom error plot.
            c_ind = diag_axis[3]
            bot_intersec = intersec_col[c_ind]
    else:
        # A color is used for the top error plot.
        c_ind = diag_axis[1]
        top_intersec = intersec_col[c_ind]

        if diag_axis[2] == 0:
            # A magnitude is used for the bottom error plot.
            m_ind = diag_axis[3]
            bot_intersec = intersec_mag[m_ind]
        else:
            # A color is used for the bottom error plot.
            c_ind = diag_axis[3]
            bot_intersec = intersec_col[c_ind]

    # Get and store values.
    for item in mag_value:

        if item <= top_intersec:
            top_val_left.append(item)
        else:
            top_val_right.append(item)

        if item <= bot_intersec:
            bot_val_left.append(item)
        else:
            bot_val_right.append(item)

    return top_val_left, top_val_right, bot_val_left, bot_val_right


def err_a_r_eyefit(mag, e_mag, e_col, params, diag_axis):
    '''
    Accept/reject stars based on an algorithm that attempts to imitate
    an 'eye fit' curve on a photometric error diagram.
    '''

    # Unpack params.
    er_params, bright_end, n_interv, interv_mag, mag_value = params
    e_max, be, be_e = er_params[1:-1]

    # Call function to obtain the median+sigmas points for magnitude
    # and color errors to fit the curves below.
    e_mag_value, e_col_value = err_med('eyefit', mag_value, e_max,
        bright_end, n_interv, interv_mag, mag, e_mag, e_col)

    # Fit polynomial + exponential curves.
    intersec_mag, intersec_col, popt_mag, pol_mag, popt_col, pol_col = \
    fit_curves(mag, mag_value, bright_end, e_mag_value, e_col_value)

    # Use the fitted curves to identify accepted/rejected stars and store
    # their indexes.
    acpt_indx, rjct_indx = separate_stars(mag, e_mag, e_col, e_max,
        bright_end, be_e, intersec_mag, intersec_col, popt_mag, pol_mag,
        popt_col, pol_col)

    # Values are used for plotting purposes only.
    top_val_left, top_val_right, bot_val_left, bot_val_right = \
    divide(mag_value, intersec_mag, intersec_col, diag_axis)
    # This list holds all the values necessary for plotting.
    popt_mag_col = [popt_mag, popt_col]
    pol_mag_col = [pol_mag, pol_col]
    err_plot = [popt_mag_col, pol_mag_col, top_val_left, top_val_right,
        bot_val_left, bot_val_right]

    return acpt_indx, rjct_indx, err_plot
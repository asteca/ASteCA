
import numpy as np
from ..inp import input_params as g
import err_medians


def fit_curves(mag, mag_value, bright_end, e_mag_value, e_col1_value):
    '''
    Get best fitting curves for the median+sigma values obtained for errors
    in magnitude and color.
    '''

    # Fit curves for errors in magnitude and color.
    for i, err_val in enumerate([e_mag_value, e_col1_value]):

        # Fit polynomial envelope of grade 2.
        pol = np.polyfit(mag_value, err_val, 2)
        # Set half magnitude value between the bright end and the maximum
        # as the point of intersection.
        intersec = (max(mag_value) + bright_end) / 2.
        # Define left max error value as the value of the polynomial fit in
        # the intersection magnitude.
        popt = np.polyval(pol, intersec)

        if i == 0:
            val_mag, pol_mag, intersec_mag = popt, pol, intersec
        elif i == 1:
            val_col1, pol_col1, intersec_col1 = popt, pol, intersec

    return intersec_mag, intersec_col1, val_mag, pol_mag, val_col1, pol_col1


def separate_stars(mag, e_mag, e_col1, be_m, intersec_mag, intersec_col1,
                   val_mag, pol_mag, val_col1, pol_col1):
    '''
    Use the curves obtained above to accept or reject stars in the
    magnitude range beyond the (brightest star + be) limit.
    '''

    e_max, be_e = g.er_params[1], g.er_params[3]

    # Initialize empty lists.
    acpt_indx, rjct_indx = [], []

    # Iterate through all stars and accept or reject those beyond
    # the (brightest star + be mag) limit according to the curve
    # obtained for the errors in magnitude and color.
    for st_ind, st_mag in enumerate(mag):

        # Reject stars with at least one error >= e_max.
        if e_mag[st_ind] >= e_max or e_col1[st_ind] >= e_max:
            rjct_indx.append(st_ind)
        else:
            # For stars brighter than the bright end.
            if mag[st_ind] <= be_m:
                # For values in this range accept all stars with both errors
                # < be_e.
                if e_mag[st_ind] < be_e and e_col1[st_ind] < be_e:
                    # Accept star.
                    acpt_indx.append(st_ind)
                else:
                    # Reject star.
                    rjct_indx.append(st_ind)

            else:
                # For the reminder of stars, we check to see if they are
                # located above or below the upper envelope for both errors.
                # If they are above in either one, we reject them, otherwise
                # we accept them.

                # Check if star is located to the left or right of the
                # intersect value for each error and compare with the
                # corresponding curve.
                mag_rjct = False
                if mag[st_ind] <= intersec_mag:
                    # Compare with linear value.
                    if e_mag[st_ind] > val_mag:
                        # Reject star.
                        mag_rjct = True
                else:
                    # Compare with polynomial.
                    if e_mag[st_ind] > np.polyval(pol_mag, (mag[st_ind])):
                        # Reject star.
                        mag_rjct = True

                col1_rjct = False
                if mag[st_ind] <= intersec_col1:
                    # Compare with linear value.
                    if e_col1[st_ind] > val_col1:
                        # Reject star.
                        col1_rjct = True
                else:
                    # Compare with polynomial.
                    if e_col1[st_ind] > np.polyval(pol_col1, (mag[st_ind])):
                        # Reject star.
                        col1_rjct = True

                if mag_rjct or col1_rjct:
                    # Reject star.
                    rjct_indx.append(st_ind)
                else:
                    # Accept star.
                    acpt_indx.append(st_ind)

    return acpt_indx, rjct_indx


def divide(mag_value, intersec_mag, intersec_col1):
    '''
    Divide magnitude interval in two, the first fitted with the linear value
    and the second with the polynomial. Do this for the errors in magnitude
    and in color.
    '''

    # Separate mag_values between those to the left and to the right of this
    # intersect value. The exponential will be the first part of the upper
    # envelope and the polynomial will be the second.
    mag_val_left, mag_val_right = [], []
    for item in mag_value:
        if item <= intersec_mag:
            mag_val_left.append(item)
        else:
            mag_val_right.append(item)

    col1_val_left, col1_val_right = [], []
    for item in mag_value:
        if item <= intersec_col1:
            col1_val_left.append(item)
        else:
            col1_val_right.append(item)

    return mag_val_left, mag_val_right, col1_val_left, col1_val_right


def main(mag, e_mag, e_col1, err_pck):
    '''
    Accept/reject stars based on an algorithm that attempts to imitate
    an 'eye fit' curve on a photometric error diagram.
    '''

    # Unpack params.
    bright_end, mag_value = err_pck[0], err_pck[2]

    # Call function to obtain the median+sigmas points for magnitude
    # and color errors to fit the curves below.
    e_mag_value, e_col1_value = err_medians.main('eyefit', err_pck, mag,
                                                 e_mag, e_col1)

    # Fit polynomial + exponential curves.
    intersec_mag, intersec_col1, val_mag, pol_mag, val_col1, pol_col1 = \
        fit_curves(mag, mag_value, bright_end, e_mag_value, e_col1_value)

    # Use the fitted curves to identify accepted/rejected stars and store
    # their indexes.
    acpt_indx, rjct_indx = separate_stars(mag, e_mag, e_col1, bright_end,
                                          intersec_mag, intersec_col1, val_mag,
                                          pol_mag, val_col1, pol_col1)

    # Values are used for plotting purposes only.
    top_val_left, top_val_right, bot_val_left, bot_val_right = \
        divide(mag_value, intersec_mag, intersec_col1)
    # This list holds all the values necessary for plotting.
    err_plot = [val_mag, pol_mag, val_col1, pol_col1, top_val_left,
                top_val_right, bot_val_left, bot_val_right]

    return acpt_indx, rjct_indx, err_plot


import numpy as np
import err_medians


def fit_curves(mag_value, bright_end, e_mag_value, e_col1_value):
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


def separate_stars(mmag, em, ec, er_params, be_m, intersec_mag, intersec_col1,
                   val_mag, pol_mag, val_col1, pol_col1):
    '''
    Use the curves obtained above to accept or reject stars in the
    magnitude range beyond the (brightest star + be) limit.
    '''

    e_max, be_e = er_params[1], er_params[3]

    # Initialize empty lists.
    acpt_indx = []

    # Iterate through all stars and accept or reject those beyond
    # the (brightest star + be mag) limit according to the curve
    # obtained for the errors in magnitude and color.
    for st_ind, st_mag in enumerate(mmag):

        # Reject stars with at least one error >= e_max.
        if em[st_ind] >= e_max or ec[st_ind] >= e_max:
            pass
        else:
            # For stars brighter than the bright end.
            if mmag[st_ind] <= be_m:
                # For values in this range accept all stars with both errors
                # < be_e.
                if em[st_ind] < be_e and ec[st_ind] < be_e:
                    # Accept star.
                    acpt_indx.append(st_ind)

            else:
                # For the reminder of stars, we check to see if they are
                # located above or below the upper envelope for both errors.
                # If they are above in either one, we reject them, otherwise
                # we accept them.

                # Check if star is located to the left or right of the
                # intersect value for each error and compare with the
                # corresponding curve.
                mag_rjct = False
                if mmag[st_ind] <= intersec_mag:
                    # Compare with linear value.
                    if em[st_ind] > val_mag:
                        # Reject star.
                        mag_rjct = True
                else:
                    # Compare with polynomial.
                    if em[st_ind] > np.polyval(pol_mag, (mmag[st_ind])):
                        # Reject star.
                        mag_rjct = True

                col1_rjct = False
                if mmag[st_ind] <= intersec_col1:
                    # Compare with linear value.
                    if ec[st_ind] > val_col1:
                        # Reject star.
                        col1_rjct = True
                else:
                    # Compare with polynomial.
                    if ec[st_ind] > np.polyval(pol_col1, (mmag[st_ind])):
                        # Reject star.
                        col1_rjct = True

                if mag_rjct or col1_rjct:
                    # Reject star.
                    pass
                else:
                    # Accept star.
                    acpt_indx.append(st_ind)

    return acpt_indx


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


def main(err_pck, cld, er_params, **kwargs):
    '''
    Accept/reject stars based on an algorithm that attempts to imitate
    an 'eye fit' curve on a photometric error diagram.
    '''

    # Unpack params.
    mmag, em, ec = cld['mags'][0], cld['em'][0], cld['ec'][0]
    bright_end, mag_value = err_pck[0], err_pck[2]

    # Call function to obtain the median+sigmas points for magnitude
    # and color errors to fit the curves below.
    e_mag_value, e_col1_value = err_medians.main('eyefit', err_pck, cld,
                                                 er_params)

    # Fit polynomial + exponential curves.
    intersec_mag, intersec_col1, val_mag, pol_mag, val_col1, pol_col1 = \
        fit_curves(mag_value, bright_end, e_mag_value, e_col1_value)

    # Use the fitted curves to identify accepted/rejected stars and store
    # their indexes.
    acpt_indx = separate_stars(
        mmag, em, ec, er_params, bright_end, intersec_mag, intersec_col1,
        val_mag, pol_mag, val_col1, pol_col1)

    # Values are used for plotting purposes only.
    top_val_left, top_val_right, bot_val_left, bot_val_right = \
        divide(mag_value, intersec_mag, intersec_col1)
    # This list holds all the values necessary for plotting.
    err_plot = [val_mag, pol_mag, val_col1, pol_col1, top_val_left,
                top_val_right, bot_val_left, bot_val_right]

    return acpt_indx, err_plot

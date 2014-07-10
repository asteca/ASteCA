"""
@author: gabriel
"""

from scipy.optimize import curve_fit
import numpy as np


# Define exponential function.
def exp_func(x, a, b, c):
    '''
    Exponential function.
    '''
    return a * np.exp(b * x) + c


def mag_median(method, mag_value, e_max, bright_end, n_interv, interv_mag,
    mag, e_mag, e_col1):
    '''
    1- Separate stars in magnitude intervals for errors of magnitude and of
    color.

    2- Obtain histograms for each interval and get median of the
    interval and sigma value using the histogram. Do this for magnitude and
    color errors.

    '''

    # Each list within the 'mag_interv' list holds all the magnitude error
    # values for all the stars in the interval 'q' which corresponds to the
    # mag value:
    # [bright_end+(interv_mag*q) + bright_end+(interv_mag*(q+1))]/2.
    # where 'q' is the index that points to the interval being filled. Idem
    # for 'col1_interv' but with color errors.
    mag_interv = [[] for _ in range(n_interv - 1)]
    col1_interv = [[] for _ in range(n_interv - 1)]

    # Iterate through all stars
    for st_ind, st_mag in enumerate(mag):

        # Use only stars below the e_max limit.
        if e_mag[st_ind] < e_max and e_col1[st_ind] < e_max and \
        st_mag > bright_end:
            # Store each star in its corresponding
            # interval in the segmented mag list which will be used to calculate
            # the curve fits for both the mag and the color.

            # Iterate through all intervals in magnitude.
            for q in range(n_interv - 1):
                # Store star's errors in corresponding interval.
                if st_mag < (bright_end + (interv_mag * (q + 1))) and \
                st_mag >= (bright_end + (interv_mag * (q))):
                    # Star falls in this interval, store its e_mag value
                    mag_interv[q].append(e_mag[st_ind])
                    # Star falls in this interval, store its e_col1 value
                    col1_interv[q].append(e_col1[st_ind])
                    break

    # We have the magnitude errors of stars beyond the (brightest
    # star + 2.) limit stored in magnitude intervals in the 'mag_interv' list
    # and the same for color errors in the 'col1_interv'. We need to find
    # the histogram for each interval and store its mean and mean+sigma in the
    # lists 'e_mag_value' and 'e_col1_value'. The values in these lists will be
    # used to fit the curve for the magnitude and color photometric errors.

    # 'e_mag_value' will hold two lists: the first one for the mean of the
    # fitted gaussians for the stars in the interval corresponding to the
    # mag_value and the second one corresponding to the mean plus one standard
    # deviation (sigma) for the same mag interval (for the upper curve). Idem
    # for the 'e_col1_value' but for the color.
    e_mag_value, e_col1_value = [], []

    # Initialize lists for median and sigma.
    median, sigma = [0.01, 0.005], [0.01, 0.005]
    # Initial value for sigma_prev.
    sigma_prev = [0.05, 0.05]

    # Iterate through all intervals (lists) in the magnitude range.
    for indx in range(n_interv - 1):

        # Iterate first for magnitude errors and then for color errors.
        for indx2, interv in enumerate([mag_interv[indx], col1_interv[indx]]):

            # Get median and standard deviation. We use the median instead of
            # the mean to protect against outliers.
            if sum(interv) == 0:
                # If there are no stars in the interval, set some default
                # values.
                median[indx2], sigma[indx2] = 0.01, 0.005
            else:
                median[indx2] = np.median(interv)

                # Generate interval's histogram.
                hist, bin_edges = np.histogram(interv, bins=50)
                bin_centres = (bin_edges[:-1] + bin_edges[1:]) / 2
                # Get number of stars in the interval and initialize the stars
                # counter.
                hist_total, hist_count = sum(hist), 0.

                # Set value for percentage of stars required to obtain 'sigma'.
                perc = 95.
                # Get the sigma value in the hist that holds perc% of the stars
                # in the interval by iterating through all bins in this hist.
                for bin_ind, bin_num in enumerate(hist):
                    # Set sigma value as centre of this bin.
                    sigma[indx2] = bin_centres[bin_ind]
                    # If condition is met, jump out of for loop.
                    if hist_count >= perc * hist_total / 100.:
                        break
                    else:
                        hist_count = hist_count + bin_num

                # Condition to keep sigma values low for low values of mag.
                if mag_value[indx] < (bright_end + 1.) and (median[indx2] +
                sigma[indx2]) > 0.05:
                    sigma[indx2] = 0.025
                elif mag_value[indx] >= bright_end + 1. and mag_value[indx] < \
                    bright_end + 2.:
                    if median[indx2] + sigma[indx2] > 0.1:
                        sigma[indx2] = 0.075

                # This condition prevents abrupt jumps in the value of sigma.
                if sigma[indx2] > sigma_prev[indx2] * 1.25:
                    sigma[indx2] = sigma_prev[indx2] * 1.25

                # Update sigma_prev value.
                sigma_prev[indx2] = sigma[indx2]

            # We obtained the median and sigma value for this interval.
            # Store just median OR median+sigma values depending on the
            # method selected.
            if method == 'lowexp':
                if indx2 == 0:
                    e_mag_value.append(median[indx2])
                else:
                    e_col1_value.append(median[indx2])

            elif method == 'eyefit':
                if indx2 == 0:
                    e_mag_value.append(median[indx2] + sigma[indx2])
                else:
                    e_col1_value.append(median[indx2] + sigma[indx2])

    return e_mag_value, e_col1_value


def fit_curves(mag, mag_value, bright_end, e_mag_value, e_col1_value):
    '''
    Get best fitting curves for the median+sigma values obtained for errors
    in magnitude and color.
    '''

    # Fit curves for errors in magnitude and color.
    for i, err_val in enumerate([e_mag_value, e_col1_value]):

        # Fit exponential envelope.
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

        if i == 0:
            popt_umag, pol_mag, intersec_mag = popt, pol, intersec
        elif i == 1:
            popt_ucol1, pol_col1, intersec_col1 = popt, pol, intersec

    return intersec_mag, intersec_col1, popt_umag, pol_mag, popt_ucol1, pol_col1


def separate_stars(mag, e_mag, e_col1, e_max, bright_end, be_e,
    intersec_mag, intersec_col1, popt_umag, pol_mag, popt_ucol1, pol_col1):
    '''
    Use the curves obtained above to accept or reject stars in the
    magnitude range beyond the (brightest star + be) limit.
    '''

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
            if mag[st_ind] <= bright_end:
                # For values in this range accept all stars with both errors
                # < be_e.
                if e_mag[st_ind] < be_e and e_col1[st_ind] < be_e:
                    # Accept star.
                    acpt_indx.append(st_ind)
                else:
                    # Reject star.
                    rjct_indx.append(st_ind)

            else:
            # For the reminder of stars, we check to see if they are located
            # above or below the upper envelope for both errors. If they are
            # above in either one, we reject them, otherwise we accept them.

                # Check if star is located to the left or right of the intersect
                # value for each error and compare with the corresponding curve.
                mag_rjct = False
                if mag[st_ind] <= intersec_mag:
                    # Compare with exponential.
                    if e_mag[st_ind] > exp_func(mag[st_ind], *popt_umag):
                        # Reject star.
                        mag_rjct = True
                else:
                    # Compare with polynomial.
                    if e_mag[st_ind] > np.polyval(pol_mag, (mag[st_ind])):
                        # Reject star.
                        mag_rjct = True

                col1_rjct = False
                if mag[st_ind] <= intersec_col1:
                    # Compare with exponential.
                    if e_col1[st_ind] > exp_func(mag[st_ind], *popt_ucol1):
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
    Divide magnitude interval in two, the first fitted with the
    exponential and the second with the polynomial. Do this for the errors in
    magnitude and in color.
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


def err_a_r_eyefit(mag, e_mag, e_col1, params):
    '''
    Accept/reject stars based on an algorithm that attempts to imitate
    an 'eye fit' curve on a photometric error diagram.
    '''

    # Unpack params.
    er_params, bright_end, n_interv, interv_mag, mag_value = params
    e_max, be, be_e = er_params[1:-1]

    # Call function to obtain the median+sigmas points for magnitude
    # and color errors to fit the curves below.
    e_mag_value, e_col1_value = mag_median('eyefit', mag_value, e_max,
        bright_end, n_interv, interv_mag, mag, e_mag, e_col1)

    # Fit polynomial + exponential curves.
    intersec_mag, intersec_col1, popt_umag, pol_mag, popt_ucol1, pol_col1 = \
    fit_curves(mag, mag_value, bright_end, e_mag_value, e_col1_value)

    # Use the fitted curves to identify accepted/rejected stars and store
    # their indexes.
    acpt_indx, rjct_indx = separate_stars(mag, e_mag, e_col1, e_max,
        bright_end, be_e, intersec_mag, intersec_col1, popt_umag, pol_mag,
        popt_ucol1, pol_col1)

    # Values are used for plotting purposes only.
    mag_val_left, mag_val_right, col1_val_left, col1_val_right = \
    divide(mag_value, intersec_mag, intersec_col1)
    # This list holds all the values necessary for plotting.
    err_plot = [popt_umag, pol_mag, popt_ucol1, pol_col1,
    mag_val_left, mag_val_right, col1_val_left, col1_val_right]

    return acpt_indx, rjct_indx, err_plot
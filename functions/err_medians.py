"""
@author: gabriel
"""

import numpy as np


def err_med(method, mag_value, e_max, bright_end, n_interv, interv_mag,
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
            if method == 'synth_clust':
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
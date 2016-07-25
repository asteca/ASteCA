
import numpy as np
from ..inp import input_params as g


def median_sigma(interv, be_m, sigma_prev):
    '''
    Get median and standard deviation. We use the median instead of
    the mean to protect against outliers.
    '''

    median = np.median(interv)

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
        sigma = bin_centres[bin_ind]
        # If condition is met, jump out of for loop.
        if hist_count >= perc * hist_total / 100.:
            break
        else:
            hist_count = hist_count + bin_num

    # This condition prevents abrupt jumps in the value of sigma.
    sigma = min(sigma, sigma_prev * 1.25)

    # Update sigma_prev value.
    sigma_prev = sigma

    return median, sigma, sigma_prev


def main(method, err_pck, mag, e_mag, e_col):
    '''
    1- Separate stars in magnitude intervals for errors of magnitude and of
    color.

    2- Obtain histograms for each interval and get median of the
    interval and sigma value using the histogram. Do this for magnitude and
    color errors.

    '''

    e_max = g.er_params[1]
    be_m, interv_mag, mag_value = err_pck

    # Each list within the 'mag_interv' list holds all the magnitude error
    # values for all the stars in the interval 'q' which corresponds to the
    # mag value:
    # [bright_end+(interv_mag*q) + bright_end+(interv_mag*(q+1))]/2.
    # where 'q' is the index that points to the interval being filled. Idem
    # for 'col1_interv' but with color errors.
    mag_interv = [[] for _ in mag_value]
    col_interv = [[] for _ in mag_value]

    # Iterate through all stars
    for st_ind, st_mag in enumerate(mag):

        # Use only stars below the e_max limit.
        if e_mag[st_ind] < e_max and e_col[st_ind] < e_max and st_mag > be_m:
            # Store each star in its corresponding
            # interval in the segmented mag list which will be used to
            # calculate the curve fits for both the mag and the color.

            # Iterate through all intervals in magnitude.
            for q in range(len(mag_value)):
                # Store star's errors in corresponding interval.
                if (be_m + interv_mag * q) <= st_mag < (be_m + interv_mag *
                                                        (q + 1)):
                    # Star falls in this interval, store its e_mag value
                    mag_interv[q].append(e_mag[st_ind])
                    # Star falls in this interval, store its e_col1 value
                    col_interv[q].append(e_col[st_ind])
                    break

    # We have the magnitude errors of stars beyond the (brightest
    # star + 2.) limit stored in magnitude intervals in the 'mag_interv' list
    # and the same for color errors in the 'col_interv'. We need to find
    # the histogram for each interval and store its mean and mean+sigma in the
    # lists 'e_mag_value' and 'e_col_value'. The values in these lists will be
    # used to fit the curve for the magnitude and color photometric errors.

    # 'e_mag_value' will hold two lists: the first one for the mean of the
    # fitted gaussians for the stars in the interval corresponding to the
    # mag_value and the second one corresponding to the mean plus one standard
    # deviation (sigma) for the same mag interval (for the upper curve). Idem
    # for the 'e_col_value' but for the color.
    e_mag_value, e_col_value = [], []

    # Initial values for median, sigma and sifma_prev.
    median, sigma, sigma_prev = 0.01, 0.005, 0.05

    # Iterate through all intervals (lists) in the magnitude range.
    for indx in range(len(mag_value)):
        # Check that list is not empty.
        if sum(mag_interv[indx]) != 0:
            median, sigma, sigma_prev = median_sigma(mag_interv[indx],
                                                     be_m, sigma_prev)

        # We obtained the median and sigma value for this interval.
        # Store just median OR median+sigma values depending on the
        # method selected.
        if method == 'synth_clust':
            e_mag_value.append(median)
        elif method == 'eyefit':
            e_mag_value.append(median + sigma)

    # Now for colors.
    for indx in range(len(mag_value)):
        if sum(col_interv[indx]) != 0:
            median, sigma, sigma_prev = median_sigma(col_interv[indx],
                                                     be_m, sigma_prev)

        if method == 'synth_clust':
            e_col_value.append(median)
        elif method == 'eyefit':
            e_col_value.append(median + sigma)

    return e_mag_value, e_col_value

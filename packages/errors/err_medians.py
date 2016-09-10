
import numpy as np


def median_sigma(interv, sigma_prev):
    '''
    Calculate the 'sigma', used by the 'eyefit' error function.
    '''
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
        # Set sigma value as center of this bin.
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

    return sigma, sigma_prev


def main(err_pck, e_max, mmag, e_mc, s_factor=0.):
    '''
    1- Separate stars in magnitude intervals for errors of magnitude and of
    color.

    2- Obtain histograms for each interval and get median of the
    interval and sigma value using the histogram. Do this for magnitude and
    color errors.
    '''

    be_m, interv_mag, mmag_interv_pts = err_pck

    # Each list within the 'mc_interv' list holds all the photometric error
    # values for all the stars in the interval 'q' which corresponds to the
    # mag value:
    # [bright_end+(interv_mag*q) + bright_end+(interv_mag*(q+1))]/2.
    # where 'q' is the index that points to the interval being filled.
    mc_interv = [[] for _ in mmag_interv_pts]

    # Iterate through all stars.
    for st_ind, st_mag in enumerate(mmag):

        # Use only stars above the bright end, and below the e_max limit.
        if be_m < st_mag and e_mc[st_ind] < e_max:
            # Store each star in its corresponding interval in the segmented
            # mag list. Will be used to calculate the curve fit.

            # Iterate through all intervals in magnitude.
            for q in range(len(mmag_interv_pts)):
                # Store star's errors in corresponding interval.
                if (be_m + interv_mag * q) <= st_mag < (be_m + interv_mag *
                                                        (q + 1)):
                    # Star falls in this interval, store its error value.
                    mc_interv[q].append(e_mc[st_ind])
                    break

    # We have the photometric errors of stars within the (be_m, e_max) range,
    # stored in magnitude intervals (from the main magnitude) in the
    # 'mc_interv' list.

    # 'e_mc_value' will hold the median photometric error for each interval
    # of the main magnitude. If the 'eyefit' module called, we add a 'sigma'
    # factor to the median.
    e_mc_value = []

    # Initial values for median, sigma and sigma_prev.
    median, sigma, sigma_prev = 0.01, 0.005, 0.05

    # Iterate through all intervals (lists) in the main magnitude range.
    for interv in mc_interv:
        # Check that list is not empty.
        if interv:
            median = median = np.median(interv)
            if s_factor != 0.:
                sigma, sigma_prev = median_sigma(interv, sigma_prev)

        # Store just median OR median+sigma values depending on the
        # module that called.
        e_mc_value.append(median + s_factor * sigma)

    return e_mc_value

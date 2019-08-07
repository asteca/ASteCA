
import numpy as np


def main(mmag, e_mc, be_m, interv_mag, n_interv):
    '''
    Store median of photometric errors for each main magnitude interval
    Do this for magnitude and color errors.
    '''
    # Each list within the 'mc_interv' list holds all the photometric error
    # values for all the stars in the interval 'q' which corresponds to the
    # mag value:
    # [bright_end+(interv_mag*q) + bright_end+(interv_mag*(q+1))]/2.
    # where 'q' is the index that points to the interval being filled.
    mc_interv = [[] for _ in range(n_interv)]

    # Iterate through all stars.
    for st_ind, st_mag in enumerate(mmag):

        # Use only stars above the bright end. All stars are already below
        # the err_max limit.
        if be_m <= st_mag:
            # Store each star in its corresponding interval in the segmented
            # mag list. Will be used to calculate the curve fit.

            # Iterate through all intervals in magnitude.
            for q in range(n_interv):
                # Store star's errors in corresponding interval.
                if (be_m + interv_mag * q) <= st_mag < (be_m + interv_mag *
                                                        (q + 1)):
                    # Star falls in this interval, store its error value.
                    mc_interv[q].append(e_mc[st_ind])
                    break

    # We have the photometric errors of stars within the (be_m, err_max) range,
    # stored in magnitude intervals (from the main magnitude) in the
    # 'mc_interv' list.

    # 'e_mc_value' will hold the median photometric error for each interval
    # of the main magnitude.
    e_mc_value = []
    # Initial value for the median.
    median = 0.0001
    # Iterate through all intervals (lists) in the main magnitude range.
    for interv in mc_interv:
        # Check that list is not empty.
        if interv:
            median = np.median(interv)
        e_mc_value.append(median)

    return e_mc_value

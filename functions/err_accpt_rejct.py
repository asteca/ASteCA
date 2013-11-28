"""
@author: gabriel
"""

from scipy.optimize import curve_fit
import numpy as np

def err_accpt_rejct(id_star, x_data, y_data, mag, e_mag, col1, e_col1):
    """    
    Accept and reject stars in and out of the cluster's boundaries according to
    a given criteria based on their photometric errors.
    """

################################################################################
# First part: separate stars in magnitud intervals for errors of magnitude and
# of color.

    # Get value of brightest and dimmest stars.
    min_mag, max_mag = min(mag), max(mag)
        
    # Define left side max limit for the box that holds the brightest stars.
    bright_end = (min_mag+2.)
        
    # Create a segmented list in magnitude.
    # Magnitude range.
    delta_mag = max_mag - bright_end
    # Width if the intervals in magnitude.
    interv_mag = 0.5
    # Number of intervals.
    n_interv = int(delta_mag/interv_mag)

    # Each list within the 'mag_interv' list holds all the magnitude error
    # values for all the stars in the interval 'q' which corresponds to the
    # mag value:
    # [bright_end+(interv_mag*q) + bright_end+(interv_mag*(q-1))]/2.
    # where 'q' is the index that points to the interval being filled. Idem
    # for 'col1_interv' but with color errors.
    mag_interv = [[] for _ in range(n_interv-1)]
    col1_interv = [[] for _ in range(n_interv-1)]
    mag_value = [((bright_end+(interv_mag*(q+1))) + \
    (bright_end+(interv_mag*(q))))/2. for q in range(n_interv-1)]
        
    # Initialize empty list to hold accepted/rejected stars.
    acpt_stars, rjct_stars = [], []
        
    # Iterate through all stars
    for st_ind, star_id in enumerate(id_star):

        # Reject stars with at least one error >= 0.3.
        if e_mag[st_ind] >= 0.3 or e_col1[st_ind] >= 0.3:

            rjct_stars.append([star_id, x_data[st_ind], y_data[st_ind], 
                               mag[st_ind], e_mag[st_ind], col1[st_ind],\
                               e_col1[st_ind]])
        
        else:
            # Accept star.
            
            # For stars brighter than the brightest plus 2 magnitudes.
            if mag[st_ind] <= bright_end:
                # For values in this range accept all stars with both errors
                # <0.1
                if e_mag[st_ind] < 0.1 and e_col1[st_ind] < 0.1:
                    # Accept star.
                    acpt_stars.append([star_id, x_data[st_ind], y_data[st_ind], 
                                      mag[st_ind], e_mag[st_ind], col1[st_ind],\
                                      e_col1[st_ind]])
                else:
                    # Reject star.
                    rjct_stars.append([star_id, x_data[st_ind], y_data[st_ind], 
                                      mag[st_ind], e_mag[st_ind], col1[st_ind],\
                                      e_col1[st_ind]])
            else:
            # For the reminder of stars, we store them in its corresponding
            # interval in the segmented mag list which will be used to calculate
            # the lower and upper curve fits for both the magnit. and the color.
                    
                # Iterate through all intervals in magnitude.
                for q in range(n_interv-1):
                    
                    if mag[st_ind] < (bright_end+(interv_mag*(q+1))) and \
                    mag[st_ind] >= (bright_end+(interv_mag*(q))):
                        # Star falls in this interval, store its e_mag value
                        mag_interv[q].append(e_mag[st_ind])
                        # Star falls in this interval, store its e_col1 value
                        col1_interv[q].append(e_col1[st_ind])
                        break
################################################################################


################################################################################
# Second part: obtain histograms for each interval and get median of the
# interval and sigma value using the histogram. Do this for magnitude and color
# errors.
                          
    # At this point I have the magnitude errors of stars beyond the (brightest
    # star + 2.) limit stored in magnitude intervals in the 'mag_interv' list
    # and the same for color errors in the 'col1_interv'. We need to find
    # the histogram for each interval and store its mean and mean+sigma in the
    # lists 'e_mag_value' and 'e_col1_value'. The values in these lists will be
    # used to fit the lower and upper curves for the magnitude and color 
    # photometric errors.
        
    # 'e_mag_value' will hold two lists: the first one for the mean of the
    # fitted gaussians for the stars in the interval corresponding to the
    # mag_value and the second one corresponding to the mean plus one standard
    # deviation (sigma) for the same mag interval (for the upper curve). Idem
    # for the 'e_col1_value' but for the color.
    e_mag_value, e_col1_value = [[], []], [[], []]        
        
    # Initialize lists for median and sigma.
    median, sigma = [0.01, 0.005], [0.01, 0.005]
    # Initial value for sigma_prev.
    sigma_prev = [0.05, 0.05]
        
    # Iterate through all intervals (lists) in the magnitude range.
    for indx in range(n_interv-1):

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
                bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
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
                    if hist_count >= perc*hist_total/100.:
                        break
                    else:
                        hist_count = hist_count + bin_num
            
                # Condition to keep sigma values low for low values of T1.
                if mag_value[indx] < (bright_end+1.) and (median[indx2] + \
                sigma[indx2]) > 0.05:
                    sigma[indx2] = 0.025
                elif mag_value[indx] >= bright_end+1. and mag_value[indx] < \
                    bright_end+2.:
                    if median[indx2] + sigma[indx2] > 0.1:
                        sigma[indx2] = 0.075
        
                # This condition prevents abrupt jumps in the value of sigma.
                if sigma[indx2] > sigma_prev[indx2]*1.25:
                    sigma[indx2] = sigma_prev[indx2]*1.25
                
                # Update sigma_prev value.
                sigma_prev[indx2] = sigma[indx2]
        
            # We obtained the median and sigma value for this interval.
            # Store median and median+sigma values.
            if indx2 == 0:
                e_mag_value[0].append(median[indx2])
                e_mag_value[1].append(median[indx2]+sigma[indx2])
            else:
                e_col1_value[0].append(median[indx2])
                e_col1_value[1].append(median[indx2]+sigma[indx2])
################################################################################


################################################################################
# Third part: get best fitting curves for the median and median+sigma values
# obtained previously for errors in magnitude and color.
            
    # After iterating though all intervals in magnitude, we fit an exponential
    # and a polynomial to these points (median and median+sigma).
    
    # Define exponential function.
    def func(x, a, b, c):
        '''
        Exponential function.
        '''
        return a * np.exp(b * x) + c
            
    # Fit curves for errors in magnitude.            
            
    # Fit lower (exponential) curve.
    popt_mag, pcov_mag = curve_fit(func, mag_value, e_mag_value[0])
        
    # Fit upper envelope. Exponential:
    popt_umag, pcov_mag = curve_fit(func, mag_value, e_mag_value[1])
    # Polynomial of grade 3:
    pol_mag = np.polyfit(mag_value, e_mag_value[1], 3)
        
    # Find point where curves intersect.
    # Initialize value in case no intersect value is found.
    intersec_mag = 0.   
    mag_x = np.linspace(bright_end, max(mag), 100)
    for x_val in mag_x:
        if np.polyval(pol_mag, (x_val)) > func(x_val, *popt_umag):
            intersec_mag = x_val
            break
        
        
    # Fit curves for errors in color.
            
    # Fit lower (exponential) curve.
    popt_col1, pcov_col1 = curve_fit(func, mag_value, e_col1_value[0])
        
    # Fit upper envelope. Exponential:
    popt_ucol1, pcov_col1 = curve_fit(func, mag_value, e_col1_value[1])
    # Polynomial of grade 3:
    pol_col1 = np.polyfit(mag_value, e_col1_value[1], 3)
        
    # Find point where curves intersect.
    # Initialize value in case no intersect value is found.
    intersec_col1 = 0.   
    mag_x = np.linspace(bright_end, max(mag), 100)
    for x_val in mag_x:
        if np.polyval(pol_col1, (x_val)) > func(x_val, *popt_ucol1):
            intersec_col1 = x_val
            break
################################################################################


################################################################################
# Fourth part: divide magnitude interval in two, the first fitted with the
# exponential and the second with the polynomial. Do this for the errors in
# magnitude and in color.
                
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
################################################################################


################################################################################
# Fifth part: use the curves obtained above to accept or reject stars in the
# magnitude range beyond the (brightest star + 2. mag) limit.
    
    # Iterate through all stars once again and now accept or reject those beyond
    # the (brightest star + 2. mag) limit according to the upper curves we
    # obtained for the errors in magnitude and color.
    for st_ind, star_id in enumerate(id_star):
        
        # Star already rejected and stored.
        if e_mag[st_ind] >= 0.3 or e_col1[st_ind] >= 0.3:
            pass
        else:
            # Stars brighter than the brightest plus 2 magnitudes also already
            # accepted or rejected.
            if mag[st_ind] <= bright_end:
                pass
            else:
            # For the reminder of stars, we check to see if they are located
            # above or below the upper envelope for both errors. If they are
            # above either one, we reject them, otherwise we accept them.
            
                # Check if star is located to the left or right of the intersect
                # value for each error and compare with the corresponding curve.
                mag_rjct = False
                if mag[st_ind] <= intersec_mag:
                    # Compare with exponential.
                    if e_mag[st_ind] > func(mag[st_ind], *popt_umag):
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
                    if e_col1[st_ind] > func(mag[st_ind], *popt_ucol1):
                        # Reject star.
                        col1_rjct = True
                else:
                    # Compare with polynomial.
                    if e_col1[st_ind] > np.polyval(pol_col1, (mag[st_ind])):
                        # Reject star.
                        col1_rjct = True
            
            
                if mag_rjct or col1_rjct:
                    # Reject star.
                    rjct_stars.append([star_id, x_data[st_ind], y_data[st_ind], 
                                      mag[st_ind], e_mag[st_ind], col1[st_ind],\
                                      e_col1[st_ind]])
                else:
                    # Accept star.
                    acpt_stars.append([star_id, x_data[st_ind], y_data[st_ind], 
                                      mag[st_ind], e_mag[st_ind], col1[st_ind],\
                                      e_col1[st_ind]])

                
    return bright_end, popt_mag, popt_umag, pol_mag, popt_col1, popt_ucol1, \
    pol_col1, mag_val_left, mag_val_right, col1_val_left, col1_val_right, \
    acpt_stars, rjct_stars

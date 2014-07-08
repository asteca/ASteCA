# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 17:41:28 2013

@author: gabriel
"""


def err_a_r_m(e_mag, e_col1, params):
    """
    Accept stars with photom errors < e_max both in mag and in color.
    """

    # Unpack params.
    e_max = params[0][1]
    print e_max

    # Initialize empty list to hold accepted/rejected stars' indexes.
    acpt_indx, rjct_indx = [], []

    # Iterate through all stars
    for st_ind in range(len(e_mag)):

        # Reject stars with at least one error >= e_max.
        if e_mag[st_ind] >= e_max or e_col1[st_ind] >= e_max:
            rjct_indx.append(st_ind)
        else:
            # Accept star.
            acpt_indx.append(st_ind)

    # Fit exponential curve for errors in magnitude.
    #mag_value = [min(zip(*acpt_stars)[3]), max(zip(*acpt_stars)[3])]
    #e_mag_value = [0.01, e_max]
    #popt_mag0, pcov_mag = curve_fit(exp_func, mag_value, e_mag_value)
    ## Use the same values for color error.
    #popt_col10 = popt_mag0
    ## Insert 'b' value into exponential function (not fitted here because
    ## otherwise the number of variables would be larger than the data points)
    #popt_mag = np.insert(popt_mag0, 1, 1.)
    #popt_col1 = np.insert(popt_col10, 1, 1.)

    # Pass empty list.
    err_plot = []
    return acpt_indx, rjct_indx, err_plot
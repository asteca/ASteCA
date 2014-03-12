# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 17:41:28 2013

@author: gabriel
"""

from scipy.optimize import curve_fit
import numpy as np


# Define exponential function.
def func(x, a, b):
    '''
    Exponential function.
    '''
    return a * np.exp(x) + b


def err_a_r_m(phot_data, er_params):
    """
    Accept stars with photom errors < e_max both in mag and in color.
    """

    id_star, x_data, y_data, mag, e_mag, col1, e_col1 = phot_data
    e_max = er_params[2]

    # Initialize empty list to hold accepted/rejected stars.
    acpt_stars, rjct_stars = [], []

    # Iterate through all stars
    for st_ind, star_id in enumerate(id_star):

        # Reject stars with at least one error >= e_max.
        if e_mag[st_ind] >= e_max or e_col1[st_ind] >= e_max:

            rjct_stars.append([star_id, x_data[st_ind], y_data[st_ind],
                               mag[st_ind], e_mag[st_ind], col1[st_ind],
                               e_col1[st_ind]])

        else:
            # Accept star.
            acpt_stars.append([star_id, x_data[st_ind], y_data[st_ind],
                               mag[st_ind], e_mag[st_ind], col1[st_ind],
                               e_col1[st_ind]])

    # Fit exponential curve for errors in magnitude.
    mag_value = [min(zip(*acpt_stars)[3]), max(zip(*acpt_stars)[3])]
    e_mag_value = [0.01, e_max]
    popt_mag0, pcov_mag = curve_fit(func, mag_value, e_mag_value)
    # Use the same values for color error.
    popt_col10 = popt_mag0
    # Insert 'b' value into exponential function (not fitted here because
    # otherwise the number of variables would be larger than the data points)
    popt_mag = np.insert(popt_mag0, 1, 1.)
    popt_col1 = np.insert(popt_col10, 1, 1.)

    return popt_mag, popt_col1, acpt_stars, rjct_stars
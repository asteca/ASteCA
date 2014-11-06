"""
@author: gabriel
"""

from functions.err_medians import err_med
from functions.exp_function import exp_func
from scipy.optimize import curve_fit
import numpy as np


# Define exponential function.
def exp_2p(x, a, b):
    '''
    Two-params exponential function.
    '''
    return a * np.exp(x) + b


def synth_clust_err(phot_data, err_pck, bf_params, da_params):
    '''
    Generate exponential error function parameters to feed the synthetic
    cluster generation function.
    '''

    da_flag, bf_flag = da_params[0], bf_params[0]

    # Check if function should run. If DA was executed, we obtain these
    # values for ploting purposes only.
    if da_flag != 'skip' or bf_flag:

        # Unpack params.
        mag, e_mag, e_col1 = phot_data[3], phot_data[4], phot_data[6]
        er_params, bright_end, n_interv, interv_mag, mag_value = err_pck
        e_max = er_params[1]

        # Call function to obtain the median points for magnitude
        # and color errors to fit the exponential curve.
        e_mag_value, e_col1_value = err_med('synth_clust', mag_value, e_max,
            bright_end, n_interv, interv_mag, mag, e_mag, e_col1)

        try:
            # Fit 3-param exponential curve.
            popt_mag, pcov_mag = curve_fit(exp_func, mag_value, e_mag_value)
            popt_col1, pcov_col1 = curve_fit(exp_func, mag_value, e_col1_value)

        # If the 3-param exponential fitting process fails.
        except RuntimeError:

            try:
                # Fit simple 2-params exponential curve.
                popt_mag, pcov_mag = curve_fit(exp_2p, mag_value, e_mag_value)
                popt_col1, pcov_col = curve_fit(exp_2p, mag_value, e_col1_value)
                # Insert 'b' value into exponential function (not fitted here
                # because otherwise the number of variables would be larger than
                # the data points)
                popt_mag = np.insert(popt_mag, 1., 1.)
                popt_col1 = np.insert(popt_col1, 1., 1.)

            # If the 2-param exponential fitting process also fails, try with a
            # 2P exp but using only two magnitude values, ie: a min and a max.
            except RuntimeError:

                # Fit simple 2-params exponential curve.
                mag_value = [min(mag), max(mag) - (max(mag) - min(mag)) / 20.]
                e_mag_value = [0.01, e_max]
                popt_mag, pcov_mag = curve_fit(exp_2p, mag_value, e_mag_value)
                # Use the same values for color error.
                popt_col1, pcov_col = curve_fit(exp_2p, mag_value, e_mag_value)
                # Insert 'b' value into exponential function.
                popt_mag = np.insert(popt_mag, 1., 1.)
                popt_col1 = np.insert(popt_col1, 1., 1.)

        err_lst = [popt_mag, popt_col1, e_max]
    else:
        err_lst = []

    return err_lst
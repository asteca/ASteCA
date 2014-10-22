"""
@author: gabriel
"""

from functions.err_medians import err_med
from functions.exp_function import exp_func
from scipy.optimize import curve_fit
import numpy as np
import get_in_params as g


# Define exponential function.
def exp_2p(x, a, b):
    '''
    Two-params exponential function.
    '''
    return a * np.exp(x) + b


def get_m_c_errors(mag, mag_value, e_max, e_mc_v):
    '''
    '''
    try:
        # Fit 3-param exponential curve.
        popt_mc, dummy = curve_fit(exp_func, mag_value, e_mc_v)

    # If the 3-param exponential fitting process fails.
    except RuntimeError:

        try:
            # Fit simple 2-params exponential curve.
            popt_mc, dummy = curve_fit(exp_2p, mag_value, e_mc_v)
            # Insert 'b' value into exponential function (not fitted here
            # because otherwise the number of variables would be larger than
            # the data points)
            popt_mc = np.insert(popt_mc, 1., 1.)

        # If the 2-param exponential fitting process also fails, try with a
        # 2P exp but using only two magnitude values, ie: a min and a max.
        except RuntimeError:

                # Fit simple 2-params exponential curve.
                mag_value = [min(mag), max(mag) - (max(mag) - min(mag)) / 20.]
                e_mc_r = [0.01, e_max]
                popt_mc, dummy = curve_fit(exp_2p, mag_value, e_mc_r)
                # Insert 'b' value into exponential function.
                popt_mc = np.insert(popt_mc, 1., 1.)

    return popt_mc


def synth_clust_err(phot_data, err_pck):
    '''
    Generate exponential error function parameters to feed the synthetic
    cluster generation function.
    '''

    da_flag, bf_flag = g.da_params[0], g.bf_params[0]
    err_lst = []

    # Check if function should run. If DA was executed, we obtain these
    # values for ploting purposes only.
    if da_flag != 'skip' or bf_flag:

        # Unpack params. Use *main* magnitude.
        mag, e_mag, e_col = phot_data[0][0], phot_data[1], phot_data[3]
        e_max, mag_value = g.er_params[1], err_pck[3]

        # Call function to obtain the median points for magnitude
        # and color errors to fit the exponential curve.
        e_mag_value, e_col_value = err_med('synth_clust', err_pck, mag, e_mag,
            e_col)

        err_lst = [[], [], e_max]

        for e_mag_v in e_mag_value:
            popt_mc = get_m_c_errors(mag, mag_value, e_max, e_mag_v)
            err_lst[0].append(popt_mc)

        for e_col_v in e_col_value:
            popt_mc = get_m_c_errors(mag, mag_value, e_max, e_col_v)
            err_lst[1].append(popt_mc)

    return err_lst
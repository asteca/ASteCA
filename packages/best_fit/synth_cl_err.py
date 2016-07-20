
from scipy.optimize import curve_fit
import numpy as np
from ..inp import input_params as g
from ..errors import err_medians
from ..math_f import exp_function


def get_m_c_errors(mag, mag_value, e_mc_v):
    '''
    Fit 3P or 2P exponential curve.
    '''

    e_max = g.er_params[1]

    try:
        # Fit 3-param exponential curve.
        popt_mc, dummy = curve_fit(exp_function.exp_3p, mag_value, e_mc_v)

    # If the 3-param exponential fitting process fails.
    except RuntimeError:

        try:
            # Fit simple 2-params exponential curve.
            popt_mc, dummy = curve_fit(exp_function.exp_2p, mag_value, e_mc_v)
            # Insert empty 'c' value to be fed later on to the 3P exponential
            # function. This makes the 2P exp function equivalent with the 3P
            # exp function, with the 'c' param 0.
            popt_mc = np.insert(popt_mc, 2, 0.)

        # If the 2-param exponential fitting process also fails, try with a
        # 2P exp but using only two magnitude values, ie: a min and a max.
        except RuntimeError:

                # Fit simple 2-params exponential curve.
                mag_value = [min(mag), max(mag) - (max(mag) - min(mag)) / 20.]
                e_mc_r = [0.01, e_max]
                popt_mc, dummy = curve_fit(exp_function.exp_2p, mag_value,
                                           e_mc_r)
                # Insert 'c' value into exponential function param list.
                popt_mc = np.insert(popt_mc, 2, 0.)

    return popt_mc


def main(phot_data, err_pck):
    '''
    Generate exponential error function parameters to feed the synthetic
    cluster generation function.
    '''

    # Unpack params.
    mag, e_mag, e_col = phot_data[3], phot_data[4], phot_data[6]
    mag_value = err_pck[2]

    # Call function to obtain the median points for magnitude
    # and color errors to fit the exponential curve.
    e_mag_value, e_col_value = err_medians.main('synth_clust', err_pck, mag,
                                                e_mag, e_col)

    err_lst = []

    for e_mag_v in [e_mag_value]:
        popt_mc = get_m_c_errors(mag, mag_value, e_mag_v)
        err_lst.append(popt_mc)

    for e_col_v in [e_col_value]:
        popt_mc = get_m_c_errors(mag, mag_value, e_col_v)
        err_lst.append(popt_mc)

    return err_lst

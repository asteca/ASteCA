
from scipy.optimize import curve_fit
import numpy as np
from ..errors import err_medians
from ..math_f import exp_function


def get_m_c_errors(mags, err_pck, e_mc_v, er_params):
    '''
    Fit 3P or 2P exponential curve.
    '''
    # List of points spanning the magnitude range starting from the bright end.
    mag_value = err_pck[2]
    # Maximum accepted photometric error value.
    e_max = er_params[1]

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
                mag_value = [min(mags), max(mags) - (max(mags) -
                             min(mags)) / 20.]
                e_mc_r = [0.01, e_max]
                popt_mc, dummy = curve_fit(exp_function.exp_2p, mag_value,
                                           e_mc_r)
                # Insert 'c' value into exponential function param list.
                popt_mc = np.insert(popt_mc, 2, 0.)

    return popt_mc


def main(cld, clp, er_params, **kwargs):
    '''
    Generate exponential error function parameters to feed the synthetic
    cluster generation function.
    '''

    # Unpack params.
    err_pck = clp['err_pck']
    mags = cld['mags']

    # Call function to obtain the median points for magnitude
    # and color errors to fit the exponential curve.
    e_mag_value, e_col_value = err_medians.main('synth_clust', err_pck, cld,
                                                er_params)

    err_lst = []

    for e_mag_v in [e_mag_value]:
        popt_mc = get_m_c_errors(mags, err_pck, e_mag_v, er_params)
        err_lst.append(popt_mc)

    for e_col_v in [e_col_value]:
        popt_mc = get_m_c_errors(mags, err_pck, e_col_v, er_params)
        err_lst.append(popt_mc)

    clp['err_lst'] = err_lst
    return clp

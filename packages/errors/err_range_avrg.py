
from scipy.optimize import curve_fit
import numpy as np
from ..errors import err_medians
from ..math_f import exp_function


def get_m_c_errors(mags, err_pck, e_mc_v, e_max):
    '''
    Fit 3P or 2P exponential curve.
    '''
    # List of points spanning the magnitude range starting from the bright end.
    mag_value = err_pck[2]
    # Maximum accepted photometric error value.

    try:
        # Fit 3-param exponential curve.
        popt_mc, dummy = curve_fit(exp_function.exp_3p, mag_value, e_mc_v)

    # If the 3-param exponential fitting process fails.
    except RuntimeError:
        print("  3P exponential error function fit failed. Attempt 2P fit.")
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
            print("  2P exponential error function fit failed."
                  " Perform min-max magnitude fit.")
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
    Fit an exponential function to the errors in each photometric dimension,
    using the main magnitude as the x coordinate.
    This data is used to display the error bars, and more importantly, to
    generate the synthetic clusters in the best match module.
    '''
    mmag, e_max = cld['mags'][0], er_params[1]

    # Call function to obtain the median points for photometric errors,
    e_mc_medians = []
    for e_mc in cld['em'].tolist() + cld['ec'].tolist():
        e_mc_medians.append(
            err_medians.main(clp['err_pck'], e_max, mmag, e_mc))

    # Fit exponential curve for each photometric error dimension.
    err_lst = []
    for e_mc_v in e_mc_medians:
        popt_mc = get_m_c_errors(mmag, clp['err_pck'], e_mc_v, e_max)
        err_lst.append(popt_mc)

    clp['err_lst'] = err_lst
    return clp

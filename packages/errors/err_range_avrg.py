
from scipy.optimize import curve_fit
import numpy as np
from ..errors import err_medians
from ..math_f import exp_function


def errorData(mmag):
    """
    Use the main magnitude to determine the error parameters that will be
    used in the synthetic cluster generation.
    """
    # Define max limit for the box that holds the brightest stars.
    be_m = min(mmag) + 2.
    # Create a segmented list in magnitude.
    # Magnitude range.
    delta_mag = max(mmag) - be_m
    # Width of the intervals in magnitude.
    interv_mag = 0.5
    # Number of intervals.
    n_interv = int(round(delta_mag / interv_mag))
    # Define list of points spanning the magnitude range starting from the
    # bright end. The '+ interv_mag' term is intentional so that the
    # err_medians function defines ranges around these values and they get
    # positioned in the middle of the magnitude interval.
    mmag_interv_pts = [
        be_m + interv_mag * (q + interv_mag) for q in range(n_interv)]

    return be_m, interv_mag, n_interv, mmag_interv_pts


def get_m_c_errors(mags, e_mc_v, err_max, mmag_interv_pts):
    '''
    Fit 3P or 2P exponential curve.
    '''
    try:
        # Fit 3-param exponential curve.
        popt_mc, dummy = curve_fit(
            exp_function.exp_3p, mmag_interv_pts, e_mc_v)

    # If the 3-param exponential fitting process fails.
    except RuntimeError:
        print("  3P exponential error function fit failed. Attempt 2P fit.")
        try:
            # Fit simple 2-params exponential curve.
            popt_mc, dummy = curve_fit(
                exp_function.exp_2p, mmag_interv_pts, e_mc_v)
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
            mmag_interv_pts = [min(mags), max(mags) - (max(mags) -
                               min(mags)) / 20.]
            e_mc_r = [0.01, err_max]
            popt_mc, dummy = curve_fit(exp_function.exp_2p, mmag_interv_pts,
                                       e_mc_r)
            # Insert 'c' value into exponential function param list.
            popt_mc = np.insert(popt_mc, 2, 0.)

    return popt_mc


def main(cld, clp, err_max, **kwargs):
    '''
    Fit an exponential function to the errors in each photometric dimension,
    using the main magnitude as the x coordinate.
    This data is used to display the error bars, and more importantly, to
    generate the synthetic clusters in the best match module.
    '''
    # Use the main magnitude.
    mmag = cld['mags'][0]

    be_m, interv_mag, n_interv, mmag_interv_pts = errorData(mmag)

    # Obtain the median points for photometric errors. Append magnitude
    # values first, and colors after.
    e_mc_medians = []
    for e_mc in cld['em'].tolist() + cld['ec'].tolist():
        e_mc_medians.append(err_medians.main(
            mmag, err_max, e_mc, be_m, interv_mag, n_interv))

    # Fit exponential curve for each photometric error dimension.
    err_lst = []
    for e_mc_v in e_mc_medians:
        popt_mc = get_m_c_errors(mmag, e_mc_v, err_max, mmag_interv_pts)
        err_lst.append(popt_mc)

    clp['err_lst'] = err_lst
    return clp


from scipy.optimize import curve_fit
import numpy as np
from ..errors import err_medians
from ..math_f import exp_function


def errorData(mmag):
    """
    Use the main magnitude to determine the error parameters that will be
    used in the synthetic cluster generation.
    """
    # Define 'bright end' leaving out the brightest stars.
    m_s = np.sort(mmag)
    # Magnitude of the .5% brightest star
    m_i = int(m_s.size * 0.005)
    be_m = max(min(mmag) + 1., m_s[m_i])
    # Magnitude range.
    delta_mag = max(mmag) - be_m
    # Width of the intervals in magnitude.
    interv_mag = 0.5
    # Number of intervals.
    n_interv = int(round(delta_mag / interv_mag))
    # Create a segmented list in magnitude.
    # Define list of points spanning the magnitude range starting from the
    # bright end. The '+ interv_mag' term is intentional so that the
    # err_medians function defines ranges around these values and they get
    # positioned in the middle of the magnitude interval.
    if n_interv < 2:
        print("  WARNING: main magnitude range is very small.")
    mmag_interv_pts = [
        be_m + interv_mag * (q + interv_mag) for q in range(n_interv)]

    return be_m, interv_mag, n_interv, mmag_interv_pts


def get_m_c_errors(mags, e_mc_v, mmag_interv_pts):
    '''
    Fit 3P or 2P exponential curve.
    '''
    try:
        if len(mmag_interv_pts) >= 3:
            # Fit 3-param exponential curve.
            popt_mc, dummy = curve_fit(
                exp_function.exp_3p, mmag_interv_pts, e_mc_v)
        else:
            # If the length of this list is 2, it means that the main
            # magnitude length is too small. If this is the case, do not
            # attempt to fit a 3 parameter exp function since it will fail.
            raise RuntimeError

    # If the 3-param exponential fitting process fails.
    except RuntimeError:
        print("  3P exponential error function fit failed. Attempt 2P fit.")
        try:
            # Fit simple 2-params exponential curve.
            popt_mc, dummy = curve_fit(
                exp_function.exp_2p, mmag_interv_pts, e_mc_v)
            # Insert empty 'c' value to be fed later on to the 3P exponential
            # function used to obtain the plotted error bars. This makes the
            # 2P exp function equivalent with the 3P exp function, with the
            # 'c' parameter equal to 0.
            popt_mc = np.insert(popt_mc, 2, 0.)

        # If the 2-param exponential fitting process also fails, try with a
        # 2P exp but using only min and max error values.
        except RuntimeError:
            print("  2P exponential error function fit failed."
                  " Perform min-max magnitude fit.")
            # Fit simple 2-params exponential curve.
            mmag_interv_pts = [
                min(mags), max(mags) - (max(mags) - min(mags)) / 20.]
            e_mc_r = [min(e_mc_v), max(e_mc_v)]
            popt_mc, dummy = curve_fit(exp_function.exp_2p, mmag_interv_pts,
                                       e_mc_r)
            # Insert 'c' value into exponential function param list.
            popt_mc = np.insert(popt_mc, 2, 0.)

    return popt_mc


def main(clp):
    '''
    Fit an exponential function to the errors in each photometric dimension,
    using the main magnitude as the x coordinate.
    This data is used to display the error bars, and more importantly, to
    generate the synthetic clusters in the best match module.
    '''
    # Use the main magnitude after max error rejection.
    mmag = np.array(list(zip(*(list(zip(*clp['acpt_stars_c']))[3])))[0])
    be_m, interv_mag, n_interv, mmag_interv_pts = errorData(mmag)

    # Obtain the median points for photometric errors. Append magnitude
    # values first, and colors after.
    e_mags = list(zip(*(list(zip(*clp['acpt_stars_c']))[4])))
    e_cols = list(zip(*(list(zip(*clp['acpt_stars_c']))[6])))
    e_mc_medians = []
    for e_mc in e_mags + e_cols:
        e_mc_medians.append(err_medians.main(
            mmag, e_mc, be_m, interv_mag, n_interv))

    # Fit exponential curve for each photometric error dimension.
    err_lst = []
    for e_mc_v in e_mc_medians:
        popt_mc = get_m_c_errors(mmag, e_mc_v, mmag_interv_pts)
        err_lst.append(popt_mc)

    clp['err_lst'] = err_lst
    return clp

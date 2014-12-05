"""
@author: gabriel
"""

from functions.exp_function import exp_3p
from scipy.optimize import curve_fit
import numpy as np
from .._in import get_in_params as g


def separate_stars(mag, e_mag, e_col1, be_m, popt_mag, popt_col1):
    '''
    Use the exponential curve obtained to accept/reject stars in the
    magnitude range beyond the (brightest star + be) limit.
    '''

    e_max, be_e = g.er_params[1], g.er_params[3]

    # Initialize empty lists.
    acpt_indx, rjct_indx = [], []

    # Iterate through all stars and accept or reject those beyond
    # the (brightest star + be mag) limit according to the curve
    # obtained for the errors in magnitude and color.
    for st_ind, st_mag in enumerate(mag):

        # Reject stars with at least one error >= e_max.
        if e_mag[st_ind] >= e_max or e_col1[st_ind] >= e_max:
            rjct_indx.append(st_ind)
        else:
            # For stars brighter than the bright end.
            if mag[st_ind] <= be_m:
                # For values in this range accept all stars with both errors
                # < be_e.
                if e_mag[st_ind] < be_e and e_col1[st_ind] < be_e:
                    # Accept star.
                    acpt_indx.append(st_ind)
                else:
                    # Reject star.
                    rjct_indx.append(st_ind)

            else:
            # For the reminder of stars, we check to see if they are located
            # above or below the exp envelope for both errors. If they are
            # above in either one, we reject them, otherwise we accept them.

                # Compare with exponential curve.
                mag_rjct, col1_rjct = False, False
                if e_mag[st_ind] > exp_3p(mag[st_ind], *popt_mag):
                    # Reject star.
                    mag_rjct = True

                if e_col1[st_ind] > exp_3p(mag[st_ind], *popt_col1):
                    # Reject star.
                    col1_rjct = True

                if mag_rjct or col1_rjct:
                    # Reject star.
                    rjct_indx.append(st_ind)
                else:
                    # Accept star.
                    acpt_indx.append(st_ind)

    return acpt_indx, rjct_indx


def err_a_r_lowexp(mag, e_mag, e_col1, be_m):
    '''
    Find the exponential fit to the photometric errors in mag and color
    and reject stars beyond the N*sigma limit.
    '''

    # Fit exponential curve for the magnitude.
    popt_mag, pcov_mag = curve_fit(exp_3p, mag, e_mag)
    # Fit exponential curve for the color.
    popt_col1, pcov_col1 = curve_fit(exp_3p, mag, e_col1)

    # Add a number of sigmas to one or more parameters of the exponential.
    sigmas_m = np.sqrt(np.diag(pcov_mag))
    sigmas_c = np.sqrt(np.diag(pcov_col1))
    N_sig = g.er_params[4]
    for i in [0, 1, 2]:
        popt_mag[i] = popt_mag[i] + float(N_sig) * sigmas_m[i]
        popt_col1[i] = popt_col1[i] + float(N_sig) * sigmas_c[i]

    # Use the fitted curves to identify accepted/rejected stars and store
    # their indexes.
    acpt_indx, rjct_indx = separate_stars(mag, e_mag, e_col1, be_m, popt_mag,
        popt_col1)

    err_plot = [popt_mag, popt_col1]

    return acpt_indx, rjct_indx, err_plot

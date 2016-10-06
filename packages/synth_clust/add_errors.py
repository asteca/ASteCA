
import numpy as np
from ..math_f import exp_function


def gauss_error(col, e_col, mag, e_mag):
    '''
    Randomly move mag and color through a Gaussian function.
    '''
    col_gauss = col + np.random.normal(0, 1, len(col)) * e_col
    mag_gauss = mag + np.random.normal(0, 1, len(col)) * e_mag

    return col_gauss, mag_gauss


def main(isoch_compl, err_lst, e_max):
    '''
    Randomly move stars according to given error distributions.
    '''

    popt_mag, popt_col = err_lst
    sigma_mag = np.array(exp_function.exp_3p(isoch_compl[1], *popt_mag))
    sigma_col = np.array(exp_function.exp_3p(isoch_compl[1], *popt_col))
    # Replace all error values greater than e_max with e_max.
    sigma_mag[sigma_mag > e_max] = e_max
    sigma_col[sigma_col > e_max] = e_max

    ###################################################################
    # # Generate errors that depend only on the theoretical isochrone.
    # b, c, max_err_mag, max_err_col = 0.25, 0.015, 0.1, 0.25
    # a1 = (max_err_mag - c) / np.exp(b * max(isoch_compl[1]))
    # a2 = (max_err_col - c) / np.exp(b * max(isoch_compl[1]))
    # sigma_mag = a1 * np.exp(b * isoch_compl[1]) + c
    # sigma_col = a2 * np.exp(b * isoch_compl[0]) + c
    ###################################################################

    # Call function to shift stars around these errors.
    col_gauss, mag_gauss = gauss_error(isoch_compl[0], sigma_col,
                                       isoch_compl[1], sigma_mag)

    isoch_error = [col_gauss, sigma_col, mag_gauss, sigma_mag]

    return isoch_error

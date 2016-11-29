
import numpy as np
from ..math_f import exp_function


def gauss_error(mc, e_mc):
    '''
    Randomly move mag and color through a Gaussian function.
    '''
    mc_gauss = mc + np.random.normal(0, 1, len(mc)) * e_mc

    return mc_gauss


def main(isoch_compl, err_lst, e_max, N_fc):
    '''
    Randomly move stars according to given error distributions.
    '''
    # Fit parameters for the 3 params exponential function, are stored in
    # err_lst with magnitude values first, and colors later.
    isoch_error = []
    for i, popt_mc in enumerate(err_lst):
        # isoch_compl[0] is the main magnitude.
        sigma_mc = np.array(exp_function.exp_3p(isoch_compl[0], *popt_mc))
        # Replace all error values greater than e_max with e_max.
        sigma_mc[sigma_mc > e_max] = e_max

        ###################################################################
        # # Generate errors that depend only on the theoretical isochrone.
        # b, c, max_err_mc = 0.25, 0.015, 0.1
        # a1 = (max_err_mc - c) / np.exp(b * max(isoch_compl[i]))
        # sigma_mc = a1 * np.exp(b * isoch_compl[i]) + c
        ###################################################################

        # Call function to shift stars around these errors.
        mc_gauss = gauss_error(isoch_compl[i], sigma_mc)

        # Create list with: 1st photometric dimension, associated error,
        # 2nd photom dimension, associated error, etc.
        isoch_error.append([mc_gauss, sigma_mc])

    # Append extra information.
    synth_clust = np.array(isoch_error + [isoch_compl[(N_fc[0] + N_fc[1]):]])

    return synth_clust

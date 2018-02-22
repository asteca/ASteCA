
import numpy as np
if __name__ == 'add_errors':
    import sys, os
    sys.path.append(os.path.abspath(os.path.join('..', '')))
    from math_f import exp_function
else:
    from ..math_f import exp_function


def gauss_error(rnd, mc, e_mc):
    '''
    Randomly move mag and color through a Gaussian function.
    '''
    mc_gauss = mc + rnd * e_mc

    return mc_gauss


def main(isoch_compl, binar_idx, err_lst, err_max, N_fc):
    """
    Randomly move stars according to given error distributions.

    Fit parameters for the 3 parameters exponential function are stored in
    'err_lst', with magnitude values first and colors second.
    """

    isoch_error = [[], []]
    i_s, i_l = 0, len(isoch_compl[0])

    # Generate once a random array with enough length to apply for all
    # photometric dimensions.
    rnd = np.random.normal(0, 1, len(isoch_compl[0]) * sum(N_fc))
    for i, popt_mc in enumerate(err_lst):
        # isoch_compl[0] is the main magnitude.
        sigma_mc = np.array(exp_function.exp_3p(isoch_compl[0], *popt_mc))
        # Replace all error values greater than err_max with err_max.
        sigma_mc[sigma_mc > err_max] = err_max

        ###################################################################
        # # Generate errors that depend only on the theoretical isochrone.
        # b, c, max_err_mc = 0.25, 0.015, 0.1
        # a1 = (max_err_mc - c) / np.exp(b * max(isoch_compl[i]))
        # sigma_mc = a1 * np.exp(b * isoch_compl[i]) + c
        ###################################################################

        # Randomly move stars around these errors.
        mc_gauss = gauss_error(
            rnd[i_s:i_l * (i + 1)], isoch_compl[i], sigma_mc)
        i_s = i_l * (i + 1)

        # Create list with photometric dimensions in first sub-list, and
        # associated errors in the second.
        isoch_error[0].append(mc_gauss)
        isoch_error[1].append(sigma_mc)

    # Append indexes that identify binaries, and extra information.
    # isoch_compl = [f1, f2, .., c1, c2, .., fc1, fc2, .., extra_pars(6)]
    synth_clust = [isoch_error] + [list(binar_idx)] +\
        [list(isoch_compl[(N_fc[0] + N_fc[1] + 2 * N_fc[1]):])]

    return synth_clust


from ..core_imp import np
if __name__ in ['add_errors', '__main__']:
    import sys, os
    sys.path.append(os.path.abspath(os.path.join('..', '')))
    from math_f import exp_function
else:
    from ..math_f import exp_function


def gauss_error(rnd, mc, e_mc):
    '''
    Randomly move mag and color through a Gaussian function.

    mc  : magnitude or color dimension
    rnd : random array of floats (0., 1.)
    e_mc: fitted observational uncertainty value
    '''
    mc_gauss = mc + rnd[:len(mc)] * e_mc

    return mc_gauss


def main(isoch_compl, err_lst, err_max, m_ini, err_rnd):
    """
    Randomly move stars according to given error distributions.

    Fit parameters for the 3 parameters exponential function are stored in
    'err_lst', with magnitude values first and colors second.
    """

    isoch_error = [[], []]

    for i, popt_mc in enumerate(err_lst):
        # isoch_compl[0] is the main magnitude.
        sigma_mc = np.array(exp_function.exp_3p(isoch_compl[0], *popt_mc))

        # Replace all photometric error values greater than err_max with
        # err_max.
        sigma_mc[sigma_mc > err_max[i]] = err_max[i]

        ###################################################################
        # # Generate errors that depend only on the theoretical isochrone.
        # b, c, max_err_mc = 0.25, 0.015, 0.1
        # a1 = (max_err_mc - c) / np.exp(b * max(isoch_compl[i]))
        # sigma_mc = a1 * np.exp(b * isoch_compl[i]) + c
        ###################################################################

        # Randomly move stars around these errors.
        mc_gauss = gauss_error(err_rnd, isoch_compl[i], sigma_mc)

        # Create list with photometric dimensions in first sub-list, and
        # associated errors in the second.
        isoch_error[0].append(mc_gauss)
        isoch_error[1].append(sigma_mc)

    # Append extra information (binary prob, binary mass, 6 extra params).
    synth_clust = [isoch_error] + [list(isoch_compl[(m_ini - 2):])]

    return synth_clust

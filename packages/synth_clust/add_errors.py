
from ..core_imp import np
from ..math_f import exp_function


def main(isoch_compl, err_pars):
    """
    Randomly move stars according to given error distributions.

    Fit parameters for the 3 parameters exponential function are stored in
    'err_lst', with magnitude values first and colors second.
    """

    err_lst, err_max, err_rnd = err_pars

    photom = []
    for i, popt_mc in enumerate(err_lst):
        # isoch_compl[0] is the main magnitude.
        sigma_mc = np.array(exp_function.exp_3p(isoch_compl[0], *popt_mc))

        # Replace all photometric error values greater than err_max with
        # err_max.
        sigma_mc[sigma_mc > err_max[i]] = err_max[i]

        # Randomly move stars around these errors.
        mc_gauss = gauss_error(err_rnd, isoch_compl[i], sigma_mc)

        # Create list with photometric dimensions in first sub-list, and
        # associated errors in the second.
        photom.append(mc_gauss)

    # Transposing is necessary for np.histogramdd()
    return np.array(photom).T


def errPlot(isoch_compl, err_pars, binar_flag, m_ini_idx):

    err_lst, err_max, err_rnd = err_pars

    photom, sigma = [], []
    for i, popt_mc in enumerate(err_lst):
        # isoch_compl[0] is the main magnitude.
        sigma_mc = np.array(exp_function.exp_3p(isoch_compl[0], *popt_mc))

        # Replace all photometric error values greater than err_max with
        # err_max.
        sigma_mc[sigma_mc > err_max[i]] = err_max[i]

        # Randomly move stars around these errors.
        mc_gauss = gauss_error(err_rnd, isoch_compl[i], sigma_mc)

        # Create list with photometric dimensions in first sub-list, and
        # associated errors in the second.
        photom.append(mc_gauss)
        sigma.append(sigma_mc)

    # The '-2' is there to include the binary probabilities and masses.
    if binar_flag:
        extra_pars = isoch_compl[m_ini_idx - 2:]
    else:
        # Placeholders for binarity data
        binar = np.zeros((2, isoch_compl.shape[1]))
        extra_pars = np.concatenate((binar, isoch_compl[m_ini_idx:]))

    # Append extra information (binary prob, binary mass, 6 extra params).
    return np.array(photom), np.array(sigma), extra_pars


def gauss_error(rnd, mc, e_mc):
    """
    Randomly move mag and color through a Gaussian function.

    mc  : magnitude or color dimension
    rnd : random array of floats normally distributed around 0. with stddev 1.
    e_mc: fitted observational uncertainty value
    """
    mc_gauss = mc + rnd[:len(mc)] * e_mc

    return mc_gauss

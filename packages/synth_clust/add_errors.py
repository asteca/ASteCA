
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


def new(isoch_compl, err_pars):
    """
    Scatter photometry given the observed uncertanties distribution.
    """

    mag, mc_gauss = err_pars[:2]

    # Find the closest indexes in 'mag' pointing to elements in
    # 'isoch_compl[0]' (ie: magnitudes)
    closest_idxs = np.searchsorted(mag, isoch_compl[0])

    # Add uncertanties to all the photometric dimensions
    # sigma_mc = sigma_mc[:, closest_idxs]
    photom = isoch_compl[:mc_gauss.shape[0]] + mc_gauss[:, closest_idxs]

    # Transposing is necessary for np.histogramdd()
    return photom.T


def newerrPlot(isoch_compl, err_pars, binar_flag, m_ini_idx):
    """
    Return more values, used for plotting and storing the synthetic cluster
    to file.
    """

    mag, mc_gauss, sigma_mc = err_pars

    # Find the closest indexes in 'mag' pointing to elements in
    # 'isoch_compl[0]' (ie: magnitudes)
    closest_idxs = np.searchsorted(mag, isoch_compl[0])

    # Add uncertanties to all the photometric dimensions
    sigma_mc = sigma_mc[:, closest_idxs]
    photom = isoch_compl[:mc_gauss.shape[0]] + mc_gauss[:, closest_idxs]

    # The '-2' is there to include the binary probabilities and masses.
    if binar_flag:
        extra_pars = isoch_compl[m_ini_idx - 2:]
    else:
        # Placeholders for binarity data
        binar = np.zeros((2, isoch_compl.shape[1]))
        extra_pars = np.concatenate((binar, isoch_compl[m_ini_idx:]))

    return photom, sigma_mc, extra_pars


def prep(cl_max_mag, err_max, err_lst, Nmag=100000):
    """
    Assigns uncertainties to magnitude values for all the photometric
    dimensions defined (magnitude + colors)

    Nmag: controls the magnitude grid. A larger value means more accurate
    assignation in 'main()', but it also has a toll on the performace. The
    value 1e3 seems ilke a reasonable compromise.
    """

    # HARDCODED +- 1. mags added to the observed magnitude range
    mags = np.array(cl_max_mag)[:, 3]
    mmin, mmax = mags.min()[0] - 1., mags.max()[0] + 1.
    mag = np.linspace(mmin, mmax, Nmag)

    mc_gauss, sigma_mc = [], []
    for i, popt_mc in enumerate(err_lst):
        # Uncertainties over the main magnitude.
        errors = np.array(exp_function.exp_3p(mag, *popt_mc))

        # Clip at 0. and err_max.
        errors = np.clip(errors, a_min=0., a_max=err_max[i])

        # Randomly move stars around these errors.
        mc_gauss.append(np.random.normal(0., 1., Nmag) * errors)
        sigma_mc.append(errors)

    return mag, np.array(mc_gauss), np.array(sigma_mc)

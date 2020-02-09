
import numpy as np
from ..math_f import exp_function


def main(isoch_compl, err_pars, binar_flag=None, m_ini_idx=None):
    """
    The 'flag' parameter is there to be used by the 'synt_cl_file.py' file.
    """

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
        # Called from 'synth_cl_file'
        if binar_flag is not None:
            sigma.append(sigma_mc)

    # If this flag is not set, it means this is being called from the
    # likelihood function and thus the 'extra_pars' are not required.
    if binar_flag is None:
        # Transposing is necessary for np.histogramdd()
        return np.array(photom).T
    else:
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


def randIdxs(lkl_method, N_errors=1000000):
    """
    Generate random indexes to use in the error function when generating the
    synthetic clusters.

    HARDCODED 'N_errors': this assumes that there will never be more than 1e6
    stars in a synthetic cluster
    """
    from .set_rand_seed import np

    if lkl_method == 'tolstoy':
        # Tolstoy likelihood considers uncertainties, there's no need to
        # add it to the synthetic clusters.
        err_rand = np.zeros(N_errors)
    else:
        err_rand = np.random.normal(0., 1., N_errors)

    return err_rand


from ..math_f import exp_function


def main(isoch_compl, err_lst, err_norm_rand):
    """
    Add random synthetic uncertainties to the magnitude and color(s)
    """

    rnd = err_norm_rand[:isoch_compl.shape[-1]]

    for i, popt_mc in enumerate(err_lst):
        # isoch_compl[0] is the main magnitude.
        sigma_mc = getSigmas(isoch_compl[0], popt_mc)

        # Randomly move stars around these errors.
        isoch_compl[i] = gauss_error(rnd, isoch_compl[i], sigma_mc)

    return isoch_compl


def getSigmas(main_mag, popt_mc):
    """
    Uncertainties for each photometric dimension
    """
    return exp_function.exp_3p(main_mag, *popt_mc)


def gauss_error(rnd, mc, e_mc):
    """
    Randomly move mag and color through a Gaussian function.

    mc  : magnitude or color dimension
    rnd : random array of floats normally distributed around 0. with stddev 1.
    e_mc: fitted observational uncertainty value
    """
    mc_gauss = mc + rnd * e_mc

    return mc_gauss

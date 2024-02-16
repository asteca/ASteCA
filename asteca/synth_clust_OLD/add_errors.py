import numpy as np


def main(isoch_compl, err_lst, rand_norm_vals):
    """
    Add random synthetic uncertainties to the magnitude and color(s)
    """

    rnd = rand_norm_vals[: isoch_compl.shape[-1]]
    main_mag = isoch_compl[0]

    for i, (a, b, c) in enumerate(err_lst):
        # isoch_compl[0] is the main magnitude.
        # sigma_mc = getSigmas(isoch_compl[0], popt_mc)

        # Three-parameters exponential function for the uncertainty
        sigma_mc = a * np.exp(b * main_mag) + c

        # Randomly move stars around these errors.
        # isoch_compl[i] = gauss_error(rnd, isoch_compl[i], sigma_mc)

        # Randomly perturb photometry with a Gaussian distribution
        isoch_compl[i] += rnd * sigma_mc

    return isoch_compl


# def getSigmas(main_mag, popt_mc):
#     """
#     Uncertainties for each photometric dimension

#     Three-parameters exponential function.
#     """
#     a, b, c = popt_mc
#     return a * np.exp(b * main_mag) + c


# def gauss_error(rnd, mc, e_mc):
#     """
#     Randomly move mag and color through a Gaussian function.

#     mc  : magnitude or color dimension
#     rnd : random array of floats normally distributed around 0. with stddev 1.
#     e_mc: fitted observational uncertainty value
#     """
#     mc_gauss = mc + rnd * e_mc

#     return mc_gauss

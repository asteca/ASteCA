import numpy as np


def main(binarity_flag, isoch_mass, alpha, beta, m_ini_idx, rand_unif_vals):
    """
    Select a fraction of stars to be binaries, given a chosen method.
    """
    # No binarity process defined
    if binarity_flag is False:
        return isoch_mass

    Ns = isoch_mass.shape[-1]
    mass = isoch_mass[m_ini_idx]

    b_p = binarProbsF(mass, alpha, beta)
    # Stars (masses) with the largest binary probabilities are selected
    # proportional to their probability
    bin_indxs = b_p > rand_unif_vals[:Ns]

    # Index of the binary magnitude: mag_binar = m_ini_idx + 1
    # Update array with new values of magnitudes and colors for the binary
    # systems.
    if bin_indxs.any():
        # for i in range(N_fc[0] + N_fc[1]):
        for i in range(m_ini_idx):
            isoch_mass[i][bin_indxs] = isoch_mass[m_ini_idx + 1 + i][bin_indxs]

    # Update the binary systems' masses so that the secondary masses for
    # SINGLE systems are identified with a '0.' value.
    isoch_mass[-1][~bin_indxs] = 0.0

    return isoch_mass


def binarProbsF(x, alpha=None, beta=None):
    """
    Distribution of probability of binarity (multiplicity fraction) versus
    primary masses.

    D&K: from Duchene and Kraus 2013.
    Source: https://stackoverflow.com/a/29359275/1391441
    """
    # if bp_vs_mass == "D&K":
    #     # xx = (0.08, .29, .96, 2.4, 7.65, 28.5, 151)
    #     # yy = (.22, .26, .46, .52, .64, .82, 1)
    #     # logx = np.log10(xx)
    #     # logy = np.log10(yy)
    #     logx = (-1.097, -0.538, -0.018, 0.380, 0.884, 1.455, 2.18)
    #     logy = (-0.658, -0.585, -0.337, -0.284, -0.194, -0.086, 0.)
    #     lin_interp = interp1d(logx, logy)
    #     b_p = np.power(10.0, lin_interp(np.log10(x)))

    # elif bp_vs_mass == "logfit":
    b_p = alpha + beta * np.log(x + 1)
    b_p = np.clip(b_p, a_min=0, a_max=1)

    return b_p

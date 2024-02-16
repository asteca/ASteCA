import numpy as np


def main(alpha, beta, m_ini_idx, rand_unif_vals, isoch_mass):
    """
    Select a fraction of stars to be binaries, given a chosen method.
    """
    # No binarity process defined
    if alpha == 0.:
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


def binarProbsF(x, alpha, beta):
    """
    Distribution of probability of binarity (multiplicity fraction) versus
    primary masses.

    D&K: from Duchene and Kraus 2013.
    """
    b_p = alpha + beta * np.log(x)
    b_p = np.clip(b_p, a_min=0, a_max=1)

    return b_p

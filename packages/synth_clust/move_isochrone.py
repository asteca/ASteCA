
import numpy as np


def main(isochrone, e, d, ext_coefs, N_fc):
    '''
    Receives an isochrone of a given age and metallicity and modifies
    its color and magnitude values according to given values for the extinction
    E(B-V) (e) and distance modulus (d).

    N_fc is the number of filters (N_fc[0]), and colors defined (N_fc[1]).

                 |------Nf----|  |------Nc-----|  |---2*Nc-----|
    isochrone = [f1, f2, .., fN, c1, c2, .., cNc, fc1, .., fc2Nc, m_ini, ..]

                 |----Nf----|  |-----Nc-----|  |-----2*Nc-----|
    ext_coefs = [cf1, .., cfNf, cc1, .., ccNc, cfc1, .., cfc2Nc]
    where:
    cfX = [a, b]  ; ccX = [[a1, b1], [a2, b2]]  ; cfcX = [a, b]
    and:
    ccm_coef = a + b / Rv = ext_coefs[0] + ext_coefs[1] / Rv
    '''
    iso_moved = []

    Rv = 3.1  # TODO fix in #170
    Av = Rv * e
    Nf, Nc = N_fc

    # Move filters.
    for fi, mag in enumerate(isochrone[:Nf]):
        # mx = Mx + dist_mod + Ax
        # Ax = cx * Av
        #
        Ax = (ext_coefs[fi][0] + ext_coefs[fi][1] / Rv) * Av
        iso_moved.append(np.array(mag) + d + Ax)

    # Move colors.
    for ci, col in enumerate(isochrone[Nf:(Nf + Nc)]):
        # (m1 - m2)o = (m1 - m2)i + E(m1 - m2)
        # E(m1 - m2) = A_m1 - A_m2
        # A_x = ef * Av ; ef = a + b/Rv (CCM model)
        # E(m1 - m2) = (ef_m1 - ef_m2) * Av
        # E(m1 - m2) = (ef_m1 - ef_m2) * Rv * E(B-V)
        #
        Ex = ((ext_coefs[Nf + ci][0][0] + ext_coefs[Nf + ci][0][1] / Rv) -
              (ext_coefs[Nf + ci][1][0] + ext_coefs[Nf + ci][1][1] / Rv)) * Av
        iso_moved.append(np.array(col) + Ex)

    # Move filters that make up the colors.
    for fci, fcol in enumerate(isochrone[(Nf + Nc):((Nf + Nc) + (Nc * 2))]):
        Ax = (ext_coefs[(Nf + Nc) + fci][0] +
              ext_coefs[(Nf + Nc) + fci][1] / Rv) * Av
        iso_moved.append(np.array(fcol) + d + Ax)

    # Append the extra parameters, not affected by distance/reddening.
    iso_moved = iso_moved + isochrone[((Nf + Nc) + (Nc * 2)):]

    return iso_moved

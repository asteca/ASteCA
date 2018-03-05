
import numpy as np


def main(isochrone, e, d, R_V, ext_coefs, N_fc):
    '''
    Receives an isochrone of a given age and metallicity and modifies
    its color and magnitude values according to given values for the extinction
    E(B-V) (e) and distance modulus (d).
    N_fc is the number of filters (N_fc[0]), and colors defined (N_fc[1]).
                 |------Nf-----|  |------Nc-----|
    isochrone = [f1, f2, .., fNf, c1, c2, .., cNc,
                 |-----Nf----|  |--------Nc------|
                 f1b, .., fNfb, c1b, c2b, .., cNcb, bp, mb, m_ini,.., m_bol
    ]
                 |----Nf----|  |-----Nc-----|
    ext_coefs = [cf1, .., cfNf, cc1, .., ccNc]
    where:
    cfX = [a, b]  ; ccX = [[a1, b1], [a2, b2]]
    and:
    ccm_coef = a + b / Rv = ext_coefs[0] + ext_coefs[1] / Rv

    Ax = ef * Av
    m_obs = M_int + Ax + dist_mod
          = M_int + ef * R_V * E(B-V) + dist_mod
          = M_int + (a + b / Rv) * R_V * E(B-V) + dist_mod

    E(m1 - m2) = A_m1 - A_m2
               = (ef_m1 - ef_m2) * Av
               = [(a1 + b1/Rv) - (a2 + b2/Rv)] * Av
               = [(a1 - a2) + (b1 - b2)/Rv] * Av
               = (a12 + b12/Rv) * Av
               = (a12 + b12/Rv) * R_V * E(B-V)
    (m1 - m2)_obs = (m1 - m2)_int + E(m1 - m2)
    (m1 - m2)_obs = (m1 - m2)_int + (a12 + b12/Rv) * R_V * E(B-V)

    '''
    Av = R_V * e
    Nf, Nc = N_fc
    iso_moved = []
    # Move filters.
    for fi, mag in enumerate(isochrone[:Nf]):
        Ax = (ext_coefs[fi][0] + ext_coefs[fi][1] / R_V) * Av
        iso_moved.append(np.array(mag) + d + Ax)

    # Move colors.
    for ci, col in enumerate(isochrone[Nf:(Nf + Nc)]):
        Ex = ((ext_coefs[Nf + ci][0][0] + ext_coefs[Nf + ci][0][1] / R_V) -
              (ext_coefs[Nf + ci][1][0] + ext_coefs[Nf + ci][1][1] / R_V)) * Av
        iso_moved.append(np.array(col) + Ex)

    # Move filters with binary data.
    for fi, mag in enumerate(isochrone[(Nf + Nc):(Nf + Nc + Nf)]):
        Ax = (ext_coefs[fi][0] + ext_coefs[fi][1] / R_V) * Av
        iso_moved.append(np.array(mag) + d + Ax)

    # Move colors with binary data.
    for ci, col in enumerate(isochrone[(Nf + Nc + Nf):(Nf + Nc + Nf + Nc)]):
        Ex = ((ext_coefs[Nf + ci][0][0] + ext_coefs[Nf + ci][0][1] / R_V) -
              (ext_coefs[Nf + ci][1][0] + ext_coefs[Nf + ci][1][1] / R_V)) * Av
        iso_moved.append(np.array(col) + Ex)

    # Append the extra parameters, not affected by distance/reddening.
    iso_moved = np.array(iso_moved + list(isochrone[(2 * Nf + 2 * Nc):]))

    return iso_moved

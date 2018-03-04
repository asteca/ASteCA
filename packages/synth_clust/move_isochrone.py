
import numpy as np


def isochMv(isochrone, e, d, R_V, ext_coefs, N_fc):
    """
    Move isochrone. Similar to the manin() function except this one only
    acts on the first sum(N_fc) sub-arrays (since the remaining extra
    parameters are not stored for the isochrone plotting process)
    """
    Av = R_V * e
    mv = (ext_coefs[0] + ext_coefs[1] / R_V) * Av + ext_coefs[2] * d
    iso_moved = isochrone + np.hstack(mv[:sum(N_fc)])

    return iso_moved


def main(isochrone, e, d, R_V, ext_coefs):
    """
    Receives an isochrone of a given age and metallicity and modifies
    its color and magnitude values according to given values for the extinction
    E(B-V) (e) and distance modulus (d).

    N_fc is the number of filters (N_fc[0]), and colors defined (N_fc[1]).

                 |------Nf-----|  |------Nc-----|
    isochrone = [f1, f2, .., fNf, c1, c2, .., cNc,

                 |-----Nf----|  |--------Nc------|
                 f1b, .., fNfb, c1b, c2b, .., cNcb, bp, mb, m_ini,.., m_bol
    ]

    ext_coefs = [a_coefs, b_coefs, dist_coeff]
    where:
               |----Nf----|    |----Nc----|   |----Nf----|   |----Nc----|
    a_coefs = [af1, .., afNf, ac1, .., acNc, af1, .., afNf, ac1, .., acNc
               |--------N_extra_pars-------|
               0., 0., 0., 0., 0., 0., 0., 0.]
    with: N_extra_pars is the number of parameters for each isochrone minus
          the magnitudes and colors times two (to account fot the binarity
          duplication)
               |----Nf----|    |----Nc----|   |----Nf----|   |----Nc----|
    b_coefs = [bf1, .., bfNf, bc1, .., bcNc, bf1, .., bfNf, bc1, .., bcNc
               |--------N_extra_pars-------|
               0., 0., 0., 0., 0., 0., 0., 0.]
                      |-(Nf-1)+Nc-|       |-(Nf-1)+Nc+N_extra_pars-|
    dist_coeff = [1., 0., ......, 0., 1., 0.,0., ................, 0.]

    such that:
    ef = a + b / Rv

    Ax = ef * Av
    m_obs = M_int + Ax + dist_mod
    m_obs = M_int + ef * R_V * E(B-V) + dist_mod

    E(m1 - m2) = A_m1 - A_m2
               = (ef_m1 - ef_m2) * Av
               = [(a1 + b1/Rv) - (a2 + b2/Rv)] * Av
               = [(a1 - a2) + (b1 - b2)/Rv] * Av
               = (a12 + b12/Rv) * Av
               = (a12 + b12/Rv) * R_V * E(B-V)
    (m1 - m2)_obs = (m1 - m2)_int + E(m1 - m2)
    (m1 - m2)_obs = (m1 - m2)_int + (a12 + b12/Rv) * R_V * E(B-V)
    """

    # Av = R_V * e
    mv = (ext_coefs[0] + ext_coefs[1] / R_V) * R_V * e + ext_coefs[2] * d
    iso_moved = isochrone + mv

    return iso_moved

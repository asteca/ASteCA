
import numpy as np


def main(
    isochrone, e, d, R_V, ext_coefs, N_fc, rand_unif, m_ini_idx,
        binar_flag, ext_diff=.0):
    """
    Receives an isochrone of a given age and metallicity and modifies
    its color and magnitude values according to given values for the extinction
    E(B-V) (e) and distance modulus (d).
    N_fc is the number of filters (N_fc[0]), and colors defined (N_fc[1]).

                 |------Nf-----|         |------Nc-----|
    isochrone = [f1, f2, .., fNf, M_ini, c1, c2, .., cNc, M_b]

                 |------Nf-----|  |------Nc-----|
    ext_coefs = [cf1, ...., cfNf, cc1, ...., ccNc]

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
    """
    # Copy to avoid overwriting
    isochrone = np.array(isochrone)

    # TODO: 'ext_diff' is i n place for #174

    Av = R_V * (e + rand_unif[:isochrone.shape[-1]] * ext_diff)
    Nf, Nc = N_fc

    def magmove(fi):
        Ax = (ext_coefs[fi][0] + ext_coefs[fi][1] / R_V) * Av
        return d + Ax

    def colmove(ci):
        Ex = ((ext_coefs[ci][0][0] + ext_coefs[ci][0][1] / R_V)
              - (ext_coefs[ci][1][0] + ext_coefs[ci][1][1] / R_V)) * Av
        return Ex

    # Move filters.
    for fi in range(Nf):
        Ax_d = magmove(fi)
        isochrone[fi] += Ax_d
        if binar_flag:
            # Move filters with binary data.
            isochrone[Nf + Nc + 1 + fi] += Ax_d

    # Move colors.
    for ci in range(Nf, Nf + Nc):
        Ex = colmove(ci)
        isochrone[ci] += Ex
        if binar_flag:
            # Move colors with binary data.
            isochrone[m_ini_idx + Nf + ci] += Ex

    # # import matplotlib.pyplot as plt
    # # # plt.subplot(121)
    # # plt.scatter(isochrone[1], isochrone[0], c='g')
    # # plt.scatter(isochrone[m_ini_idx + 2], isochrone[m_ini_idx + 1], c='r', alpha=.5)
    # # plt.gca().invert_yaxis()
    # # # Second color
    # # # plt.subplot(122)
    # # # plt.scatter(isochrone[2], isochrone[0], c='g')
    # # # plt.scatter(isochrone[6], isochrone[4], c='r', alpha=.5)
    # # # plt.gca().invert_yaxis()
    # # plt.show()

    return isochrone

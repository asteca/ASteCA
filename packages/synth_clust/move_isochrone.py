
import numpy as np


def isochMv(isochrone, e, d, R_V, ext_coefs, N_fc):
    """
    Obtain values used to move the magnitudes and colors.
    """
    Av = R_V * e
    mv = (ext_coefs[0] + ext_coefs[1] / R_V) * Av + ext_coefs[2] * d
    iso_moved = isochrone + np.hstack(mv[:sum(N_fc)])

    return iso_moved


def main(isochrone, e, d, R_V, ext_coefs):
    """
    Obtain values used to move the magnitudes and colors.
    """

    Av = R_V * e
    # Nf, Nc = N_fc

    # Method 1
    # # Magnitudes
    # # mx = Mx + dist_mod + cx * Av
    # # A_mv = dist_mod + cx * Av
    # #
    # A_mv = d + (ext_coefs[0][:, 0] + ext_coefs[0][:, 1] / R_V) * Av

    # # Colors
    # # (m1 - m2)o = (m1 - m2)i + E(m1 - m2)
    # # E(m1 - m2) = A_m1 - A_m2
    # # A_x = ef * Av ; ef = a + b/R_V (CCM model)
    # # E(m1 - m2) = (ef_m1 - ef_m2) * Av
    # # E(m1 - m2) = (ef_m1 - ef_m2) * R_V * E(B-V)
    # #
    # E_mv = ((ext_coefs[1][:, 0][:, 0] + ext_coefs[1][:, 0][:, 1] / R_V) -
    #         (ext_coefs[1][:, 1][:, 0] + ext_coefs[1][:, 1][:, 1] / R_V)) * Av

    # # Store two times: once for the original magnitudes and colors, and once
    # # for the binary values of these magnitudes and colors.
    # A_E_mv = np.vstack(np.hstack([A_mv, E_mv, A_mv, E_mv]))

    # iso_moved = np.append(
    #     isochrone[:len(A_E_mv)] + A_E_mv, isochrone[len(A_E_mv):], axis=0)

    # Method 2
    # a1 = np.vstack(d + (ext_coefs[0][:, 0] + ext_coefs[0][:, 1] / R_V) * Av)
    # a2 = np.vstack(
    #     ((ext_coefs[1][:, 0][:, 0] + ext_coefs[1][:, 0][:, 1] / R_V) -
    #      (ext_coefs[1][:, 1][:, 0] + ext_coefs[1][:, 1][:, 1] / R_V)) * Av)

    # isochrone[:Nf] = isochrone[:Nf] + a1
    # isochrone[Nf:(Nf + Nc)] = isochrone[Nf:(Nf + Nc)] + a2
    # isochrone[(Nf + Nc):(Nf + Nc + Nf)] =\
    #     isochrone[(Nf + Nc):(Nf + Nc + Nf)] + a1
    # isochrone[(Nf + Nc + Nf):(Nf + Nc + Nf + Nc)] =\
    #     isochrone[(Nf + Nc + Nf):(Nf + Nc + Nf + Nc)] + a2

    # # Method 3
    # a1 = d + (ext_coefs[0][:, 0] + ext_coefs[0][:, 1] / R_V) * Av
    # a2 = ((ext_coefs[1][:, 0][:, 0] + ext_coefs[1][:, 0][:, 1] / R_V) -
    #       (ext_coefs[1][:, 1][:, 0] + ext_coefs[1][:, 1][:, 1] / R_V)) * Av

    # a3 = np.concatenate((a1, a2, a1, a2))[:, None]

    # # This one works
    # # iso_moved = list(isochrone[:2 * sum(N_fc)] + a3)
    # # iso_moved = np.array(iso_moved + list(isochrone[(2 * sum(N_fc)):]))

    # # This one does not
    # isochrone[:2 * sum(N_fc)] = isochrone[:2 * sum(N_fc)] + a3
    # iso_moved = isochrone

    # Method 4
    mv = (ext_coefs[0] + ext_coefs[1] / R_V) * Av + ext_coefs[2] * d
    iso_moved = isochrone + mv

    return iso_moved


def main_old(isochrone, e, d, R_V, ext_coefs, N_fc):
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
    '''
    Av = R_V * e
    Nf, Nc = N_fc
    iso_moved = []
    # Move filters.
    for fi, mag in enumerate(isochrone[:Nf]):
        # mx = Mx + dist_mod + Ax
        # Ax = cx * Av
        #
        Ax = (ext_coefs[fi][0] + ext_coefs[fi][1] / R_V) * Av
        iso_moved.append(np.array(mag) + d + Ax)

    # Move colors.
    for ci, col in enumerate(isochrone[Nf:(Nf + Nc)]):
        # (m1 - m2)o = (m1 - m2)i + E(m1 - m2)
        # E(m1 - m2) = A_m1 - A_m2
        # A_x = ef * Av ; ef = a + b/R_V (CCM model)
        # E(m1 - m2) = (ef_m1 - ef_m2) * Av
        # E(m1 - m2) = (ef_m1 - ef_m2) * R_V * E(B-V)
        #
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
    iso_moved = iso_moved + list(isochrone[(2 * Nf + 2 * Nc):])

    return iso_moved


if __name__ == '__main__':
    import time as t
    import sys
    import pickle
    mpath = sys.path[0].replace('synth_clust', '').replace('packages/', '')
    with open(mpath + 'moved.pickle', 'rb') as f:
        isochrone, e, d, R_V, ext_coefs, N_fc = pickle.load(f)

    ext_coefs_new = [
        np.array([[9.99749022e-01], [-3.57235098e-02], [1.76371324e-04],
                  [9.99749022e-01], [-3.57235098e-02], [1.76371324e-04],
                  [0.00000000e+00], [0.00000000e+00], [0.00000000e+00],
                  [0.00000000e+00], [0.00000000e+00], [0.00000000e+00],
                  [0.00000000e+00], [0.00000000e+00]]),
        np.array([[-0.00462922], [0.83868144], [0.95016114], [-0.00462922],
                  [0.83868144], [0.95016114], [0.], [0.], [0.], [0.], [0.],
                  [0.], [0.], [0.]]),
        np.array([[1.], [0.], [0.], [1.], [0.], [0.], [0.], [0.], [0.], [0.],
                 [0.], [0.], [0.], [0.]])]

    N = 100000

    s = t.clock()
    for _ in range(N):
        iso_moved = main(isochrone, e, d, R_V, ext_coefs, N_fc)
    t1 = t.clock() - s

    s = t.clock()
    for _ in range(N):
        iso_moved2 = main_n(isochrone, e, d, R_V, ext_coefs_new)
    t2 = t.clock() - s

    print(t1, t2)
    print(np.all([np.allclose(
        x1, x2) for x1, x2 in zip(*[iso_moved, iso_moved2])]))

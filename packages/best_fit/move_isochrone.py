
import numpy as np
from ..inp import input_params as g


def main(isochrone, e, d):
    '''
    Receives an isochrone of a given age and metallicity and modifies
    its color and magnitude values according to given values for the extinction
    E(B-V) (e) and distance modulus (d).
    '''

    cmd_sel = g.ps_params[1]

    iso_moved = [[], []]

    if cmd_sel == 1:
        # For UBVI system.
        #
        # E(B-V) = (B-V) - (B-V)o
        # Av = 3.1*E(B-V)
        # (mv - Mv)o = -5 + 5*log(d) + Av
        #
        Av = 3.1 * e
        iso_moved = [np.array(isochrone[0]) + e,
                     np.array(isochrone[1]) + d + Av]

    elif cmd_sel == 2:
        # For UBVI system.
        #
        # E(V-I) = (V-I) - (V-I)o
        # E(V-I) = 1.244*E(B-V)
        # Av = 3.1*E(B-V)
        # (mv - Mv)o = -5 + 5*log(d) + Av
        #
        Av = 3.1 * e
        iso_moved = [np.array(isochrone[0]) + 1.244 * e,
                     np.array(isochrone[1]) + d + Av]
    elif cmd_sel == 3:
        # For UBVI system.
        #
        # E(U-B) = 0.72 * E(B-V) + 0.05 * E(B-V)^2
        # E(U-B) = (U-B) - (U-B)o
        # Av = 3.1*E(B-V)
        # (mv - Mv)o = -5 + 5*log(d) + Av
        #
        Av = 3.1 * e
        iso_moved = [np.array(isochrone[0]) + (0.72 * e + 0.05 * e ** 2),
                     np.array(isochrone[1]) + d + Av]

    elif cmd_sel == 4:
        # For Washington system.
        #
        # E(C-T1) = 1.97*E(B-V) = (C-T1) - (C-T)o
        # M_T1 = T1 + 0.58*E(B-V) - (m-M)o - 3.2*E(B-V)
        #
        # (C-T1) = (C-T1)o + 1.97*E(B-V)
        # T1 = M_T1 - 0.58*E(B-V) + (m-M)o + 3.2*E(B-V)
        #
        V_Mv = d + 3.2 * e
        iso_moved = [np.array(isochrone[0]) + 1.97 * e,
                     np.array(isochrone[1]) - 0.58 * e + V_Mv]

    elif cmd_sel == 5:
        # For 2MASS system.
        #
        # E(J-H) = 0.34*E(B-V) = (J-H) - (J-H)o
        # A_J = A_V - 2.28*E(B-V) = 0.82*E(B-V)
        #
        # (J-H) = (J-H)o + 0.34*E(B-V)
        # J = M_J + (m-M)o + A_J
        #
        A_J = 0.82 * e
        iso_moved = [np.array(isochrone[0]) + 0.34 * e,
                     np.array(isochrone[1]) + d + A_J]

    elif cmd_sel == 6:
        # For 2MASS system.
        #
        # E(J-H) = 0.34*E(B-V) = (J-H) - (J-H)o
        # A_H = A_J - 0.34*E(B-V) = 0.82*E(B-V) - 0.34*E(B-V) = 0.48*E(B-V)
        #
        # (J-H) = (J-H)o + 0.34*E(B-V)
        # H = M_H + (m-M)o + A_H
        #
        A_H = 0.48 * e
        iso_moved = [np.array(isochrone[0]) + 0.34 * e,
                     np.array(isochrone[1]) + d + A_H]

    elif cmd_sel == 7:
        # For 2MASS system.
        #
        # E(H-K) = 0.2*E(B-V) = (H-K) - (H-K)o
        # A_K = A_H - 0.2*E(B-V) = 0.48*E(B-V) - 0.2*E(B-V) = 0.28*E(B-V)
        #
        # (H-K) = (H-K)o + 0.2*E(B-V)
        # K = M_K + (m-M)o + A_K
        #
        A_K = 0.28 * e
        iso_moved = [np.array(isochrone[0]) + 0.2 * e,
                     np.array(isochrone[1]) + d + A_K]

    elif cmd_sel == 8:
        # For SDSS ugriz system.
        #
        # E(u-g) = 1.125951*E(B-V)
        # A_g = 1.20585*Av = 1.20585*3.1*E(B-V) = 3.738135*E(B-V)
        #
        A_g = 3.738135 * e
        iso_moved = [np.array(isochrone[0]) + 1.125951 * e,
                     np.array(isochrone[1]) + d + A_g]

    elif cmd_sel == 9:
        # For SDSS ugriz system.
        #
        # E(g-r) = 1.037353*E(B-V)
        # A_g = 1.20585*Av = 1.20585*3.1*E(B-V) = 3.738135*E(B-V)
        #
        A_g = 3.738135 * e
        iso_moved = [np.array(isochrone[0]) + 1.037353 * e,
                     np.array(isochrone[1]) + d + A_g]

    elif cmd_sel == 10:
        # For Stroemgren system.
        # From: Carreta et al. (2011) http://arxiv.org/abs/1109.3199
        #
        # E(b-y) = 0.772*E(B-V)
        # A_y = 3.277*E(B-V) (\sim 1.*A_V)
        #
        A_y = 3.277 * e
        iso_moved = [np.array(isochrone[0]) + 0.772*e,
                     np.array(isochrone[1]) + d + A_y]

    elif cmd_sel == 11:
        # For Stroemgren system.
        #
        # E(v-y) = 1.275*E(B-V)
        #
        A_y = 3.277 * e
        iso_moved = [np.array(isochrone[0]) + 1.275*e,
                     np.array(isochrone[1]) + d + A_y]

    elif cmd_sel == 12:
        # For Stroemgren system.
        #
        # E(u-y) = 1.954*E(B-V)
        #
        A_y = 3.277 * e
        iso_moved = [np.array(isochrone[0]) + 1.954*e,
                     np.array(isochrone[1]) + d + A_y]

    elif cmd_sel == 13:
        # For HST/ACS WFC system. Use average values for cool/hot stars from
        # http://adsabs.harvard.edu/abs/2005MNRAS.357.1038B
        #
        # E(F606W-F814W) = 1.072*E(B-V)
        #
        A_F606W = 2.6875 * e
        iso_moved = [np.array(isochrone[0]) + 0.921*e,
                     np.array(isochrone[1]) + d + A_F606W]

    return iso_moved

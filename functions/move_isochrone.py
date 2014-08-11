# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 16:35:24 2014

@author: gabriel
"""

import numpy as np


def move_isoch(cmd_sel, isochrone, e, d):
    '''
    Recieves an isochrone of a given age and metallicity and modifies
    its color and magnitude values according to given values for the extinction
    E(B-V) (e) and distance modulus (d).
    '''
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

    return iso_moved
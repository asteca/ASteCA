# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 16:35:24 2014

@author: gabriel
"""

import numpy as np
import get_in_params as g


def move_isoch(isochrone, e, d, ccm_coefs):
    '''
    Recieves an isochrone of a given age and metallicity and modifies
    its color and magnitude values according to given values for the extinction
    E(B-V) (e) and distance modulus (d).
    '''

    # Store initial and actual mass.
    iso_moved = [isochrone[0], isochrone[1], [], []]

    # Rv ratio is defined in the input params file.
    Av = g.ps_params[2] * e

    # Loop through magnitudes.
    for mi, mag in enumerate(isochrone[2]):

        # mx = Mx + dist_mod + Ax
        # Ax = ccm_m * Av
        #
        Ax = ccm_coefs[0][mi] * Av
        iso_moved[2].append(np.array(mag) + d + Ax)

    # Loop through colors.
    for ci, col in enumerate(isochrone[3]):

        # (C1 - C2) = (C1 - C2)o + Ax
        # Ax = ccm_c * Av
        #
        Ax = ccm_coefs[1][ci] * Av
        iso_moved[3].append(np.array(col) + Ax)

        ## For UBVI system.
        ##
        ## E(V-I) = (V-I) - (V-I)o
        ## E(V-I) = 1.244*E(B-V)
        ## Av = 3.1*E(B-V)
        ## (mv - Mv)o = -5 + 5*log(d) + Av
        ##
        #iso_moved = [np.array(isochrone[0]) + 1.244 * e,
                     #np.array(isochrone[1]) + d + Av]

        ## For UBVI system.
        ##
        ## E(U-B) = 0.72 * E(B-V) + 0.05 * E(B-V)^2
        ## E(U-B) = (U-B) - (U-B)o
        ## Av = 3.1*E(B-V)
        ## (mv - Mv)o = -5 + 5*log(d) + Av
        ##
        #iso_moved = [np.array(isochrone[0]) + (0.72 * e + 0.05 * e ** 2),
                     #np.array(isochrone[1]) + d + Av]

        ## For Washington system.
        ##
        ## E(C-T1) = 1.97*E(B-V) = (C-T1) - (C-T)o
        ## M_T1 = T1 + 0.58*E(B-V) - (m-M)o - 3.2*E(B-V)
        ##
        ## (C-T1) = (C-T1)o + 1.97*E(B-V)
        ## T1 = M_T1 - 0.58*E(B-V) + (m-M)o + 3.2*E(B-V)
        ##
        #V_Mv = d + 3.2 * e
        #iso_moved = [np.array(isochrone[0]) + 1.97 * e,
                     #np.array(isochrone[1]) - 0.58 * e + V_Mv]

        ## For 2MASS system.
        ##
        ## E(J-H) = 0.34*E(B-V) = (J-H) - (J-H)o
        ## A_J = A_V - 2.28*E(B-V) = 0.82*E(B-V)
        ##
        ## (J-H) = (J-H)o + 0.34*E(B-V)
        ## J = M_J + (m-M)o + A_J
        ##
        #A_J = 0.82 * e
        #iso_moved = [np.array(isochrone[0]) + 0.34 * e,
                     #np.array(isochrone[1]) + d + A_J]

        ## For 2MASS system.
        ##
        ## E(J-H) = 0.34*E(B-V) = (J-H) - (J-H)o
        ## A_H = A_J - 0.34*E(B-V) = 0.82*E(B-V) - 0.34*E(B-V) = 0.48*E(B-V)
        ##
        ## (J-H) = (J-H)o + 0.34*E(B-V)
        ## H = M_H + (m-M)o + A_H
        ##
        #A_H = 0.48 * e
        #iso_moved = [np.array(isochrone[0]) + 0.34 * e,
                     #np.array(isochrone[1]) + d + A_H]

        ## For 2MASS system.
        ##
        ## E(H-K) = 0.2*E(B-V) = (H-K) - (H-K)o
        ## A_K = A_H - 0.2*E(B-V) = 0.48*E(B-V) - 0.2*E(B-V) = 0.28*E(B-V)
        ##
        ## (H-K) = (H-K)o + 0.2*E(B-V)
        ## K = M_K + (m-M)o + A_K
        ##
        #A_K = 0.28 * e
        #iso_moved = [np.array(isochrone[0]) + 0.2 * e,
                     #np.array(isochrone[1]) + d + A_K]

    return iso_moved
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 16:35:24 2014

@author: gabriel
"""

def move_isoch(sys_select, isochrone, e, d):
    '''
    Recieves an isochrone of a given age and metallicity and modifies
    it according to given values for the extinction E(B-V) and distance
    modulus.
    '''
    iso_moved = [[], []]
    
    if sys_select == 'UBVI':
        # For UBVI system.
        #
        # E(B-V) = (B-V) - (B-V)o
        # Av = 3.1*E(B-V)
        # (mv - Mv)o = -5 + 5*log(d) + Av
        #
        Av = 3.1*e
        for item in isochrone[1]:
            # mv affected by extinction.
            iso_moved[1].append(item + d + Av)
        for item in isochrone[0]:
            # (B-V) affected by extinction.
            iso_moved[0].append(item + e)
    elif sys_select == 'WASH':
        # For Washington system.
        #
        # E(C-T1) = 1.97*E(B-V) = (C-T1) - (C-T)o
        # M_T1 = T1 + 0.58*E(B-V) - (m-M)o - 3.2*E(B-V)
        #
        # (C-T1) = (C-T1)o + 1.97*E(B-V)
        # T1 = M_T1 - 0.58*E(B-V) + (m-M)o + 3.2*E(B-V)
        #
        V_Mv = d + 3.2*e
        for item in isochrone[1]:
             # T1 magnitude affected by extinction.
            iso_moved[1].append(item - 0.58*e + V_Mv)
        for item in isochrone[0]:
             # C-T1 color affected by extinction.
            iso_moved[0].append(item + 1.97*e)
            
    # Store moved colors and magnitudes in final list.
    isoch_final = [iso_moved[0], iso_moved[1]]

    return isoch_final
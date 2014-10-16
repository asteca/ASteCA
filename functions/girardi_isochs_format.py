# -*- coding: utf-8 -*-
"""
Created on Tue Sep 23 13:03:39 2014

@author: gabriel
"""

import get_in_params as g
from cmd_phot_systs import phot_mags as pm


def age_f():
    '''
    Define reg expression to isolate the age of an isochrone.
    '''
    age_format = r"Age = \t(.+?) yr"
    return age_format


def i_format(syst):
    '''
    Read line start format and columns indexes for the selected set of
    Girardi isochrones and chosen CMD.
    '''

    iso_select = g.ps_params[0]

    # Assign values according to the system and set of isochrones selected.
    if iso_select == 'MAR':
        # String that identifies the beginning of a new isochrone.
        line_start = "#\tIsochrone\tZ ="
        # Mass columns and index where magnitudes begin.
        mass_i, mass_a, mags_i = 1, 2, 7
    elif iso_select == 'PAR':
        # String that identifies the beginning of a new isochrone.
        line_start = "#\tIsochrone  Z = "
        # Mass columns and index where magnitudes begin.
        mass_i, mass_a, mags_i = 2, 3, 8

    # Get photometric systems dictionary.
    all_systs = pm()

    # Identify indexes of magnitudes defined in the input paramrter file
    # as stored in its respective photometric system's metallicity files.
    mags = []
    # Name of photometric system as defined in the dictionary above.
    sys = syst[0]
    # For each magnitude defined in this photometric system.
    for mag in syst[1]:
        # Store the correct column index for this magnitude.
        mags.append(mags_i + all_systs[sys].index(mag))

    return line_start, mass_i, mass_a, mags
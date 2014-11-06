# -*- coding: utf-8 -*-
"""
Created on Tue Sep 23 13:03:39 2014

@author: gabriel
"""


def isoch_format(iso_select, syst):
    '''
    Read line start format and columns indexes for the selected set of
    Girardi isochrones and chosen CMD.
    '''

    # Define reg expression to isolate the age of an isochrone.
    age_format = r"Age = \t(.+?) yr"

    # Assign values according to the system and set of isochrones selected.
    if iso_select == 'MAR':
        # String that identifies the beginning of a new isochrone.
        line_start = "#\tIsochrone\tZ ="
        # Mass column.
        mass_i, mass_a = 1, 2
    elif iso_select == 'PAR':
        # String that identifies the beginning of a new isochrone.
        line_start = "#\tIsochrone  Z = "
        # Mass column.
        mass_i, mass_a = 2, 3

    # Dictionary that stores the names and column indexes for each
    # magnitude defined in each phoyometric system.
    all_systs = {
        'UBVRIJKH': (8, ['U', 'B', 'V', 'R', 'I', 'J', 'H', 'K']),
        'JHKs': (8, ['J', 'H', 'Ks']),
        'CMT1T2BVRI': (8, ['C', 'M', 'T1', 'T2', 'B', 'V', 'R', 'I'])}

    # Identify indexes of magnitudes defined in the input paramrter file
    # as stored in its respective photometric system's metallicity files.
    mags = []
    # Name of photometric system as defined in the dictionary above.
    sys = syst[0]
    # For each magnitude defined in this photometric system.
    for mag in syst[1]:
        # Store the correct column index for this magnitude.
        mags.append(all_systs[sys][0] + all_systs[sys][1].index(mag))

    return line_start, age_format, mass_i, mass_a, mags
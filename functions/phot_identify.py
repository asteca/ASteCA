# -*- coding: utf-8 -*-
"""
Created on Thu Oct 09 11:36:00 2014

@author: gabriel
"""

import get_in_params as g
from cmd_phot_systs import phot_mags as pm


def identify_phot_data():
    '''
    Reads the photometric data given in the input file by the user and
    separates the magnitude and colors indexes (and their errors)
    from their photometric system identities.
    '''

    # Unpack data.
    mags, cols, phot_diag = g.gd_params[1:]

    magnitudes, e_mags, colors, e_cols = [], [], [], []
    mag_names, col_names, all_mags = [], [], []

    for mag in mags:
        colum_indx, phot_name = mag.split(',')
        if not phot_name.startswith("e"):
            magnitudes.append(int(colum_indx))
            mag_names.append(phot_name)
            all_mags.append(phot_name)
        else:
            e_mags.append(int(colum_indx))

    for col in cols:
        colum_indx, phot_name = col.split(',')
        if not phot_name.startswith("e"):
            colors.append(int(colum_indx))
            col_names.append(phot_name)
            # Separate magnitudes from color name.
            m1, m2 = phot_name.split('-')
            all_mags.append(phot_name[0] + m1[1:])
            all_mags.append(phot_name[0] + m2)
        else:
            e_cols.append(int(colum_indx))

    # Remove duplicate magnitudes if they exist.
    all_mags = list(set(all_mags))

    # Initialize for up to 9 photometric systems.
    phot_mags = [[] for _ in range(9)]
    for mag in all_mags:
        # Extract full magnitude name.
        m = mag[1:] if len(mag) > 2 else mag[1]
        phot_mags[int(mag[0])].append(m)

    # Identifiers for different photometric systems.
    all_systs = pm()
    # Create list with names of photometric systems.
    phot_systs = []
    for i, phot_syst in enumerate(phot_mags):
        if phot_syst:
            # Store names of folders that hold the Girardi isochrones that
            # correspond to this photometric system.
            phot_systs.append([all_systs.keys()[i], phot_syst])

    # Identify photometric data to plot.
    diag_axis = []
    for axis in phot_diag:
        if axis[:3] == 'mag':
            diag_axis.append(0)  # 0 points to magnitudes list.
            diag_axis.append(int(axis[-1]) - 1)
        else:
            diag_axis.append(2)  # 2 points to colors list.
            diag_axis.append(int(axis[-1]) - 1)
    #
    # diag_axis = [index pointing to either mag or col, index pointing to
    # which mag or col, index pointing to mag or col, index pointing to which
    # mag or col]

    phot_indexes = [magnitudes, colors, e_mags, e_cols]
    phot_names = [mag_names, col_names]

    return [phot_indexes, phot_names, phot_systs, diag_axis]
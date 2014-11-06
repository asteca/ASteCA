# -*- coding: utf-8 -*-
"""
Created on Thu Oct 09 11:36:00 2014

@author: gabriel
"""

import get_in_params as g


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
            all_mags.append(phot_name[0] + phot_name[1])
            all_mags.append(phot_name[0] + phot_name[2])
        else:
            e_cols.append(int(colum_indx))

    # Remove duplicate magnitudes if they exist.
    all_mags = list(set(all_mags))

    # Initialize for up to 20 photometric systems.
    phot_systs = [[] for _ in range(20)]
    for mag in all_mags:
        phot_systs[int(mag[0]) - 1].append(mag[1])

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

    print magnitudes, colors, e_mags, e_cols
    print mag_names, col_names
    print phot_systs
    print diag_axis

    ## Fix isochrones location according to the CMD and set selected.
    #text1, text2 = 'none', 'none'
    #if cmd_select in {1, 2, 3, 8}:
        #text1 = 'ubvi'
    #elif cmd_select in {4}:
        #text1 = 'wash'
    #elif cmd_select in {5, 6, 7}:
        #text1 = '2mass'
    #text2 = 'marigo' if iso_select == 'MAR' else 'parsec'
    ## Set iso_path according to the above values.
    #iso_path = join(mypath + '/isochrones/' + text1 + '_' + text2)

    phot_indexes = [magnitudes, colors, e_mags, e_cols]
    phot_names = [mag_names, col_names]

    return [phot_indexes, phot_names, phot_systs, diag_axis]
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 09 11:36:00 2014

@author: gabriel
"""

import get_in_params as gip


def identify_phot_data():
    '''
    Reads the photometric data given in the input file by the user and
    separates the magnitude and colors indexes (and their errors)
    from their photometric system identities.
    '''

    # Unpack data.
    mags, cols, phot_diag = gip.gd_params[1:]

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
            all_mags.extend((phot_name[0], phot_name[1]))
        else:
            e_cols.append(int(colum_indx))

    # Remove duplicate magnitudes if they exist.
    all_mags = list(set(all_mags))

    print magnitudes, colors
    print e_mags, e_cols
    print mag_names, col_names
    print all_mags
    raw_input()

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

    ## Fix magnitude and color names for the CMD axis.
    ## m_1 is the y axis magnitude, m_2 is the magnitude used to obtain the
    ## color index and the third value in each key indicates how the color
    ## is to be formed, e.g: '12' means (m_1 - m_2)
    #cmds_dic = {1: ('V', 'B', 21), 2: ('V', 'I', 12), 3: ('V', 'U', 21),
        #4: ('{T_1}', 'C', 21), 5: ('J', 'H', 12), 6: ('H', 'J', 21),
        #7: ('K', 'H', 21), 8: ('(B-V)', '(U-B)', 3)}


    return
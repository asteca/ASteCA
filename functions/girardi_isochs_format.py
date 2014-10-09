# -*- coding: utf-8 -*-
"""
Created on Tue Sep 23 13:03:39 2014

@author: gabriel
"""


def isoch_format(iso_select, cmd_select):
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
        mass = 1
        if cmd_select == 1:
            # V, B
            mags = [9, 8]
        elif cmd_select == 2:
            # V, I
            mags = [9, 11]
        elif cmd_select == 3:
            # B, U
            mags = [8, 7]
        elif cmd_select == 4:
            # T1, C
            mags = [9, 7]
    elif iso_select == 'PAR':
        # String that identifies the beginning of a new isochrone.
        line_start = "#\tIsochrone  Z = "
        # Mass column.
        mass = 2
        if cmd_select == 1:
            # V, B
            mags = [10, 9]
        elif cmd_select == 2:
            # V, I
            mags = [10, 12]
        if cmd_select == 3:
            # V, U
            mags = [10, 8]
        elif cmd_select == 4:
            # T1, C
            mags = [10, 8]
        elif cmd_select == 5:
            # J, H
            mags = [8, 9]
        elif cmd_select == 6:
            # H, J
            mags = [9, 8]
        elif cmd_select == 7:
            # K_s, H
            mags = [10, 9]
        elif cmd_select == 8:
            # U, B, V
            mags = [8, 9, 10]

    return line_start, age_format, mass, mags
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
    if iso_select in ['PAR10', 'PAR11', 'PAR12']:
        # String that identifies the beginning of a new isochrone.
        line_start = "#\tIsochrone  Z = "
        # Mass column.
        mass = 2
        if cmd_select == 1:
            # V, B
            mag1, mag2 = 10, 9
        elif cmd_select == 2:
            # V, I
            mag1, mag2 = 10, 12
        if cmd_select == 3:
            # V, U
            mag1, mag2 = 10, 8
        elif cmd_select == 4:
            # T1, C
            mag1, mag2 = 10, 8
        elif cmd_select == 5:
            # J, H
            mag1, mag2 = 8, 9
        elif cmd_select == 6:
            # H, J
            mag1, mag2 = 9, 8
        elif cmd_select == 7:
            # K_s, H
            mag1, mag2 = 10, 9

    return line_start, age_format, mass, mag1, mag2
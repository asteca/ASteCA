# -*- coding: utf-8 -*-
"""
Created on Fri Dec 5 12:00:00 2014

@author: gabriel
"""

import numpy as np
import re
from .._in import get_in_params as g
from girardi_isochs_format import isoch_format as i_format


def read_met_file(met_f, age_values):
    '''
    Read a given metallicity file and return the isochrones for the ages
    within the age range.
    '''

    cmd_select = g.ps_params[1]

    # Read line start format and columns indexes for the selected set of
    # Girardi isochrones.
    line_start, age_format, imass_idx, mag1_idx, mag2_idx = i_format()

    # Initialize list that will hold all the isochrones for this
    # metallicity value.
    metal_isoch = []

    # Open the metallicity file.
    with open(met_f, mode="r") as f_iso:

        # Define empty lists.
        isoch_col, isoch_mag, isoch_mas = [], [], []

        # Initial value for age to avoid 'not defined' error.
        age = -99.

        # Iterate through each line in the file.
        for line in f_iso:

            # Identify beginning of a defined isochrone.
            if line.startswith(line_start):

                # Save stored values if these exist.
                # Skip first age for which the lists will be empty.
                if isoch_col:
                    # Store color, magnitudes and masses for this
                    # isochrone.
                    metal_isoch.append([isoch_col, isoch_mag,
                                        isoch_mas])
                    # Reset lists.
                    isoch_col, isoch_mag, isoch_mas = [], [], []

                # Read age value for this isochrone.
                age0 = re.findall(age_format, line)  # Find age in line.
                age = np.around(np.log10(float(age0[0])), 2)

            # If age value falls inside the given range, store the
            # isochrone's data.
            if age in age_values:

                # Save mag, color and mass values for each isochrone.
                if not line.startswith("#"):
                    reader = line.split()
                    # Color.
                    # Generate colors correctly <-- HARDCODED, FIX
                    if cmd_select in {2, 5, 9, 13}:
                        isoch_col.append(float(reader[mag1_idx]) -
                                         float(reader[mag2_idx]))
                    else:
                        isoch_col.append(float(reader[mag2_idx]) -
                                         float(reader[mag1_idx]))
                    # Magnitude.
                    isoch_mag.append(float(reader[mag1_idx]))
                    # Mass
                    isoch_mas.append(float(reader[imass_idx]))

        # Save the last isochrone when EOF is reached.
        else:
            # If list is not empty.
            if isoch_col:
                # Store colors, magnitudes and masses for this
                # isochrone.
                metal_isoch.append([isoch_col, isoch_mag, isoch_mas])

    return metal_isoch


def get_isochs(met_f_filter, age_values):
    '''
    Stores the available isochrones of different metallicities and
    ages, according to the ranges given to these parameters.
    '''

    # Lists that store the colors, magnitudes and masses of the isochrones.
    # isoch_list = [metal_1, ..., metal_M]
    # metal_i = [isoch_i1, ..., isoch_iN]
    # isoch_ij = [colors, magnitudes, mass]
    # isoch_list[i][j] --> i: metallicity index ; j: age index
    isoch_list = []

    # Iterate in order through all the metallicity files stored for the
    # selected set of isochrones.
    for met_f in met_f_filter:

        metal_isoch = read_met_file(met_f, age_values)

        # Store list holding all the isochrones with the same metallicity
        # in the final isochrone list.
        isoch_list.append(metal_isoch)

    return isoch_list

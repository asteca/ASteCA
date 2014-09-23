# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 15:08:39 2014

@author: gabriel
"""

import os
from os.path import join
import numpy as np


def get_isoch_format(iso_select, cmd_select):
    '''
    Read line start format and columns indexes for the selected set of
    Girardi isochrones and chosen CMD.
    '''

    # Assign values according to the system and set of isochrones selected.
    if iso_select == 'MAR':
        # String that identifies the beginning of a new isochrone.
        line_start = "#\tIsochrone\tZ ="
        # Mass column.
        mas_i = 1
        if cmd_select == 1:
            # V, B
            mag1_i, mag2_i = 9, 8
        elif cmd_select == 2:
            # V, I
            mag1_i, mag2_i = 9, 11
        elif cmd_select == 3:
            # B, U
            mag1_i, mag2_i = 8, 7
        elif cmd_select == 4:
            # T1, C
            mag1_i, mag2_i = 9, 7
    elif iso_select == 'PAR':
        # String that identifies the beginning of a new isochrone.
        line_start = "#\tIsochrone  Z = "
        # Mass column.
        mas_i = 2
        if cmd_select == 1:
            # V, B
            mag1_i, mag2_i = 10, 9
        elif cmd_select == 2:
            # V, I
            mag1_i, mag2_i = 10, 12
        if cmd_select == 3:
            # V, U
            mag1_i, mag2_i = 10, 8
        elif cmd_select == 4:
            # T1, C
            mag1_i, mag2_i = 10, 8
        elif cmd_select == 5:
            # J, H
            mag1_i, mag2_i = 8, 9
        elif cmd_select == 6:
            # H, J
            mag1_i, mag2_i = 9, 8
        elif cmd_select == 7:
            # K_s, H
            mag1_i, mag2_i = 10, 9

    return line_start, mas_i, mag1_i, mag2_i


def get_ranges(m_rs, a_rs, e_rs, d_rs):
    '''
    Calculate parameter ranges to be used by the selected best fit method.
    '''

    # Store ranges and steps for these parameters.
    z_min, z_max, z_step = m_rs
    age_min, age_max, age_step = a_rs
    e_bv_min, e_bv_max, e_bv_step = e_rs
    dis_mod_min, dis_mod_max, dis_mod_step = d_rs

    # UPDATE max values.
    # Add a small value to each max value to ensure that the range is a bit
    # larger than the one between the real min and max values. This simplifies
    # the input of data and ensures that the GA algorithm won't fail when
    # encoding/decoding the floats into their binary representations.
    z_max = z_max + min(z_max / 100., z_step / 2.)
    age_max = age_max + min(age_max / 100., age_step / 2.)
    e_bv_max = e_bv_max + min(e_bv_max / 100., e_bv_step / 2.)
    dis_mod_max = dis_mod_max + min(dis_mod_max / 100., dis_mod_step / 2.)

    # Store min, *UPDATED* max values and steps.
    ranges_steps = [[z_min, z_max, z_step], [age_min, age_max, age_step],
                    [e_bv_min, e_bv_max, e_bv_step],
                    [dis_mod_min, dis_mod_max, dis_mod_step]]

    # Create ranges for metallicity and age.
    z_range = np.arange(z_min, z_max, z_step)
    a_range = np.arange(age_min, age_max, age_step)

    # Store all possible extinction and distance modulus values in list.
    # isoch_ed = [extinction, dis_mod]
    # extinction = [e_1, e_2, ..., e_n]
    # dis_mod = [dm_1, dm_2, ..., dm_m]
    e_range = [round(i, 2) for i in np.arange(e_bv_min, e_bv_max, e_bv_step)]
    d_range = [round(i, 2) for i in np.arange(dis_mod_min, dis_mod_max,
        dis_mod_step)]
    isoch_ed = [e_range, d_range]

    return z_range, a_range, isoch_ed, ranges_steps


def read_met_file(iso_path, met_file, metal, a_range, cmd_select,
    isoch_format):
    '''
    Read a given metallicity file and return the isochrones for the ages
    within the age range.
    '''

    line_start, mini_indx, mag1_indx, mag2_indx = isoch_format

    # Initialize list that will hold all the isochrones for this
    # metallicity value.
    metal_isoch = []

    # Open the metallicity file.
    with open(join(iso_path, met_file), mode="r") as f_iso:

        # List that holds the each value of metallicity and age in a
        # given single metallicity file.
        met_params = []

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
                    # Save metallicity and age in list.
                    met_params.append([metal, age])
                    # Store colors, magnitudes and masses for this
                    # isochrone.
                    metal_isoch.append([isoch_col, isoch_mag,
                        isoch_mas])
                    # Reset lists.
                    isoch_col, isoch_mag, isoch_mas = [], [], []

                # Read age value for this isochrone.
                age_str = line.split("Age =")[1]
                age = round(np.log10(float(age_str[:-3])), 2)

            # If age value falls inside the given range, store the
            # isochrone's data (in this line).
            if np.isclose(a_range, age, atol=0.01).any():

                # Save mag, color and mass values for each isochrone
                # star.
                if not line.startswith("#"):
                    reader = line.split()
                    # Color.
                    # Generate colors correctty,
                    if cmd_select in {2, 5}:
                        isoch_col.append(float(reader[mag1_indx]) -
                        float(reader[mag2_indx]))
                    else:
                        isoch_col.append(float(reader[mag2_indx]) -
                        float(reader[mag1_indx]))
                    # Magnitude.
                    isoch_mag.append(float(reader[mag1_indx]))
                    # Mass
                    isoch_mas.append(float(reader[mini_indx]))

        # Save the last isochrone when EOF is reached.
        else:
            # If list is not empty.
            if isoch_col:
                # Save metallicity and age in list.
                met_params.append([metal, age])
                # Store colors, magnitudes and masses for this
                # isochrone.
                metal_isoch.append([isoch_col, isoch_mag, isoch_mas])

    return metal_isoch, met_params


def get_ip(ps_params, metal_files):
    '''
    Stores the available isochrones of different metallicities and
    ages, according to the ranges given to these parameters.
    '''

    # Unpack.
    iso_path, cmd_select, iso_select, m_rs, a_rs, e_rs, d_rs = ps_params

    # Read line start format and columns indexes for the selected set of
    # Girardi isochrones.
    isoch_format = get_isoch_format(iso_select, cmd_select)

    # Get parameters ranges.
    z_range, a_range, isoch_ed, ranges_steps = get_ranges(m_rs, a_rs, e_rs,
        d_rs)

    # Lists that store the colors, magnitudes and masses of the isochrones.
    # isoch_list = [metal_1, ..., metal_M]
    # metal_i = [isoch_i1, ..., isoch_iN]
    # isoch_ij = [colors, magnitudes, mass]
    # isoch_list[i][j] --> i: metallicity index ; j: age index
    isoch_list = []

    # List that will hold all the metallicities and ages associated with the
    # stored isochrones.
    # This way of storing does not need to assume that
    # each metallicity file contains the same number of isochrones with the
    # same ages.
    # isoch_ma = [met_1, ..., met_M]
    # met_i = [params_i1, ..., params_iN]
    # params_ij = [metallicity_i, age_j]
    # isoch_ma[i][*][0] --> metallicity i (float)
    # isoch_ma[i][j][1] --> age j (float)
    isoch_ma = []

    # Iterate in order through all the metallicity files stored for the
    # selected set of isochrones.
    for met_file in metal_files:

        # Extract metallicity value from the name of the file.
        metal = float(met_file[:-4])

        # Process metallicity file only if it's inside the given range.
        if np.isclose(z_range, metal, atol=0.0001).any():

            metal_isoch, met_params = read_met_file(iso_path, met_file, metal,
                a_range, cmd_select, isoch_format)

            # Store list holding all the isochrones with the same metallicity
            # in the final isochrone list.
            isoch_list.append(metal_isoch)
            # Store parameter values in list that holds all the metallicities
            # and ages.
            isoch_ma.append(met_params)

    return isoch_list, isoch_ma, isoch_ed, ranges_steps


def get_metal_files(iso_path):
    '''
    Read names of all metallicity files stored in isochrones path given.
    '''

    return sorted(os.listdir(iso_path))


def ip(ps_params, bf_flag):
    '''
    Read isochrones and parameters if best fit function is set to run.
    '''

    ip_list = []
    # Only read files of best fit method is set to run.
    if bf_flag is True:

        # Unpack.
        iso_path = ps_params[0]

        # Read names of all metallicity files stored in isochrones path given.
        metal_files = get_metal_files(iso_path)

        # Get isochonres, ranges and parameters values.
        ip_list = get_ip(ps_params, metal_files)

        total = len(ip_list[0]) * len(ip_list[0][0])
        print 'Theoretical isochrones read and stored ({}).'.format(total)

    return ip_list
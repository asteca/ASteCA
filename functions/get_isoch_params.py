# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 15:08:39 2014

@author: gabriel
"""

import os
from os.path import join
import numpy as np

    

def ip(ps_params):
    '''
    Stores all the available isochrones of different metallicities and
    ages according to the ranges given to this parameters.
    Also stores available parameter values for the stored isochrones
    between the specified ranges and with the given steps.
    '''
    
    # Paths to isochrone files.
    iso_path = ps_params[0]
    # String that identifies the beginning of a new isochrone.
    line_start = ps_params[1]
    # Read indexes for this Girardi output file.
    mini_indx, mag_indx, col_indx = ps_params[2]
    # Read ranges and steps for these parameters.
    z_min, z_max, z_step = ps_params[3]
    age_min, age_max, age_step = ps_params[4]
    e_bv_min, e_bv_max, e_bv_step = ps_params[5]
    dis_mod_min, dis_mod_max, dis_mod_step = ps_params[6]
    
    # Store ranges and steps.
    ranges_steps = [[z_min, z_max, z_step], [age_min, age_max, age_step],
                    [e_bv_min, e_bv_max, e_bv_step],\
                    [dis_mod_min, dis_mod_max, dis_mod_step]]
    
    # Add a small value to each max value to ensure that the range is a bit
    # larger than the one between the real min and max values. This simplifies
    # the input of data and ensures that the GA won't fail when encoding/
    # decoding the floats into their binary representations.
    z_max = z_max + z_step/2.
    age_max = age_max + age_step/2.
    e_bv_max = e_bv_max + e_bv_step/2.
    dis_mod_max = dis_mod_max + dis_mod_step/2.
   
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
    
    # Iterate through all metallicity files in order.
    for met_file in sorted(os.listdir(iso_path)):

        # Store the metallicity value.
        metal = float(met_file[:-4])
        
        # Process metallicity file only if it's inside the given range.
        if z_min<= metal <= z_max:
            
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
                            metal_isoch.append([isoch_col, isoch_mag, isoch_mas])
                            # Reset lists.
                            isoch_col, isoch_mag, isoch_mas = [], [], []
                            
                        # Read age value.
                        age_str = line.split("Age =")[1]
#                        age = float(age_str[:-3])/1.e09
                        age = round(np.log10(float(age_str[:-3])), 2)
                        
                    # Store age value in 'ages' list if it falls inside
                    # the given range.
                    if age_min<= age <=age_max:

                        # Save mag, color and mass values for each isochrone
                        # star.
                        if not line.startswith("#"):
                            reader = line.split()
                            # Color.
                            isoch_col.append(float(reader[col_indx]) -
                            float(reader[mag_indx]))
                            # Magnitude.
                            isoch_mag.append(float(reader[mag_indx]))
                            # Mass
                            isoch_mas.append(float(reader[mini_indx]))

            # Store list holding all the isochrones with the same metallicity
            # in the final isochrone list.
            isoch_list.append(metal_isoch)
            # Store parameter values in list that holds all the metallicities
            # and ages.
            isoch_ma.append(met_params)
                  
    
    # Store all possible extinction and distance modulus values in list.
    # isoch_ed = [extinction, dis_mod]
    # extinction = [e_1, e_2, ..., e_n]
    # dis_mod = [dm_1, dm_2, ..., dm_m]
    isoch_ed = [[], []]
    for e_bv in np.arange(e_bv_min, e_bv_max, e_bv_step):
        isoch_ed[0].append(round(e_bv, 2))
    for dis_mod in np.arange(dis_mod_min, dis_mod_max, dis_mod_step):
        # Store params for this isochrone.
        isoch_ed[1].append(round(dis_mod, 2))                  
            
    ip_list = [isoch_list, isoch_ma, isoch_ed, ranges_steps]
    print isoch_ma[17], '\n'
    print isoch_ed
            
    return ip_list

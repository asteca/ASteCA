# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 15:08:39 2014

@author: gabriel
"""

import os
from os.path import join

from get_mass_dist import mass_dist as md
from genetic_algorithm import gen_algor as g_a
import numpy as np



def get_ranges_paths(sys_select, iso_select):
    '''
    Reads parameters ranges and paths to stored isochrone files.
    '''
    # Girardi isochrones span a range of ~ 3.98+06 to 1.26e+10 yr.
    # MASSCLEAN clusters span a range of ~ 1.00e06 to 1.00e+10 yr.
    if sys_select == 'UBVI':
        # Range of values where the parameters will move.
        e_bv_min, e_bv_max, e_bv_step = 0., 1.01, 0.05
        dis_mod_min, dis_mod_max, dis_mod_step = 8., 13., 0.5
        z_min, z_max = 0.01, 0.03
        age_min, age_max = 0.003, 3.2
        
        # Path where isochrone files are stored.
        iso_path = '/media/rest/github/isochrones/iso_ubvi_marigo'
        
        # Data for files formatted for UBVI Marigo tracks.
        line_start = "#\tIsochrone\tZ = "
        # Index that points to the corresponding column in the file.
        mag_indx = 9
         
    elif sys_select == 'WASH':
        # Select cloud.
#        cloud = raw_input('SMC or LMC cloud?')
#        cloud = 'LMC'
#        
        # Range of values where the parameters will move.
        e_bv_min, e_bv_max, e_bv_step = 0., 0.41, 0.01
#        if cloud == 'SMC':
#            dis_mod_min, dis_mod_max, dis_mod_step = 18.9, 18.91, 1.
#        elif cloud == 'LMC':
#            dis_mod_min, dis_mod_max, dis_mod_step = 18.5, 18.51, 1.
        dis_mod_min, dis_mod_max, dis_mod_step = 18., 19.01, 0.05
        z_min, z_max = 0.0005, 0.021
        # age_val x10 yr
#        age_min, age_max = 0.003, 12.6
        age_min, age_max = 6.6, 10.11
    
        # Select Marigo or PARSEC tracks.        
        if iso_select == 'MAR':  # Marigo.
            line_start = "#\tIsochrone\tZ = "
            # Index that points to the corresponding column in the file.
            mini_indx, col_indx, mag_indx = 1, 7, 9
            # Path where isochrone files are stored.
            iso_path = '/media/rest/github/isochrones/iso_wash_marigo'
        elif iso_select == 'PAR':  # PARSEC.
            line_start = "#\tIsochrone  Z = "
            # Index that points to the corresponding column in the file.
            mini_indx, col_indx, mag_indx = 2, 8, 10
            # Path where isochrone files are stored.
            iso_path = '/media/rest/github/isochrones/iso_wash_parsec'
    
    # Store ranges and steps.
    ranges_steps = [[z_min, z_max], [age_min, age_max],
                    [e_bv_min, e_bv_max, e_bv_step],\
                    [dis_mod_min, dis_mod_max, dis_mod_step]]
    # Store indexes that point to columns in the file that stores isochrones.
    indexes = [mini_indx, col_indx, mag_indx]
    
    return ranges_steps, indexes, iso_path, line_start
    
    

def get_isoch_params(sys_select, iso_select):
    '''
    Reads and stores available parameter values for the stored isochrones
    between the specified ranges and with the given steps.
    Also stores all the available isochrones of different metallicities
    according to the ranges given to this parameters.
    '''

    # Call function to obtain the ranges and steps for the parameters along
    # with the path to the isochrone files and information about how they
    # are formatted (line_start).
    ranges_steps, indexes, iso_path, line_start = get_ranges_paths(sys_select,
                                                                   iso_select)
    z_min, z_max = ranges_steps[0]
    age_min, age_max = ranges_steps[1]
    e_bv_min, e_bv_max, e_bv_step = ranges_steps[2]
    dis_mod_min, dis_mod_max, dis_mod_step = ranges_steps[3]
    
    # Read indexes for this Girardi output file.
    mini_indx, col_indx, mag_indx = indexes
    
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
    
    # Iterate through all metallicity files.
    for met_file in os.listdir(iso_path):

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
            # Store in list that holds all the metallicities and ages.
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
                  
    return isoch_list, isoch_ma, isoch_ed, ranges_steps



def boostrap_resample(memb_prob_avrg_sort):
    '''
    Resamples the observed cluster to use in the bootstrap process.
    '''
#    http://www.astroml.org/_modules/astroML/resample.html#bootstrap
    obs_clust = memb_prob_avrg_sort
    return obs_clust
    


def gip(sys_select, iso_select, memb_prob_avrg_sort):
    '''
    Main function.
    
    Perform a best fitting process to find the cluster's parameters:
    E(B-V), distance (distance modulus), age and metallicity.
    '''

    # Select photometric system to be used.
    # sys_select: Select UBVI or Washington system as 1 or 2.
    
    # Store mass distribution used to produce a synthetic cluster based on
    # a given theoretic isochrone.
    # 1st param: 'chabrier_2001', 'kroupa_1993'
    # 2nd param: 'total_number', 'total_mass'
    # 3rd param: total cluster mass or number of stars in clusters, depending
    # on the chosen 2nd param.
    mass_dist = md('kroupa_1993', 'total_number', 500)
    
    # Genetic algorithm parameters.
    n_pop, n_gen, fdif, p_cross, p_mut = 100, 100, 3./5., 0.85, 0.01
    
    
    # Store all isochrones in all the metallicity files in isoch_list. We do
    # this so the files will only have to be accessed once.
    # Store metallicity values and isochrones ages between the allowed
    # ranges in isoch_ma; extinction and distance modulus values in isoch_ed.
    isoch_list, isoch_ma, isoch_ed, ranges_steps = get_isoch_params(sys_select, iso_select)
#    print isoch_list[3][17]
#    print isoch_ma[3]
#    print isoch_ed
#    raw_input()
   
    # Begin bootstrap block.
    # NUmber of times to run the bootstrap block.
    N_B = 500
    # List that holds the parameters values obtained by the bootstrap
    # process.
    params_boot = []
    for i in range(N_B+1):
        
        # The first pass is done with no resampling to calculate the final
        # values. After that we resample to get the uncertainty in each
        # parameter.
        if i == 0:
            obs_clust = memb_prob_avrg_sort
        else:
            obs_clust = boostrap_resample(memb_prob_avrg_sort)

        # Call algorithm to calculate the likelihoods for the set of
        # isochrones and return the best fitting parameters.
        
        if i==0:
            # Brute force algorithm.
#            isoch_fit_params = brute_force(sys_select, isoch_params, iso_path,
#                                           line_start, indexes, obs_clust,\
#                                           mass_dist)
            # Genetic algorithm.
            isoch_fit_params = g_a(sys_select, obs_clust, isoch_list, isoch_ma, isoch_ed,
                                   mass_dist, ranges_steps, n_pop, n_gen, fdif, p_cross, p_mut)
        else:
            # Brute force.
#            params_boot.append(brute_force(sys_select, isoch_params, iso_path,
#                                           line_start, indexes, obs_clust,\
#                                           mass_dist))
            # Genetic algorithm algorithm.
            params_boot.append(g_a(sys_select, obs_clust, isoch_list, isoch_ma, isoch_ed,
                                   mass_dist, ranges_steps, n_pop, n_gen, fdif, p_cross, p_mut))
        
    # Calculate errors for each parameter.
    isoch_fit_errors = np.mean(params_boot)
    
    return isoch_fit_params, isoch_fit_errors
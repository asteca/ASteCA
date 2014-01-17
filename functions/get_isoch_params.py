# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 15:08:39 2014

@author: gabriel
"""

import os
from os.path import join

from get_mass_dist import mass_dist as md
import numpy as np
import time


def get_ranges_paths(sys_select):
    '''
    Reads parameters ranges and paths to stored isochrone files.
    '''
    # Girardi isochrones span a range of ~ 3.98+06 to 1.26e+10 yr.
    # MASSCLEAN clusters span a range of ~ 1.00e06 to 1.00e+10 yr.
    if sys_select == '1':
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
         
    elif sys_select == '2':
        # Select cloud.
#        cloud = raw_input('SMC or LMC cloud?')
#        cloud = 'LMC'
#        
        # Range of values where the parameters will move.
        e_bv_min, e_bv_max, e_bv_step = 0., 0.21, 0.02
#        if cloud == 'SMC':
#            dis_mod_min, dis_mod_max, dis_mod_step = 18.9, 18.91, 1.
#        elif cloud == 'LMC':
#            dis_mod_min, dis_mod_max, dis_mod_step = 18.5, 18.51, 1.
        dis_mod_min, dis_mod_max, dis_mod_step = 18., 19., 0.2
#        z_min, z_max = 0.0005, 0.02
        z_min, z_max = 0.0005, 0.005
        # age_val x10 yr
#        age_min, age_max = 0.003, 12.6
        age_min, age_max = 0.003, 0.5
    
        # Select Marigo or PARSEC tracks.        
#        iso_select = raw_input('Select Marigo or PARSEC tracks as 1 or 2: ')
        iso_select = '1' 
        if iso_select == '1':
            # Marigo.
            line_start = "#\tIsochrone\tZ = "
            # Index that points to the corresponding column in the file.
            mini_indx, col_indx, mag_indx = 1, 7, 9
            # Path where isochrone files are stored.
            iso_path = '/media/rest/github/isochrones/iso_wash_marigo'
        elif iso_select == '2':
            # PARSEC.
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
    
    

def get_isoch_params(sys_select):
    '''
    Reads and stores available parameter values for the stored isochrones
    between the specified ranges and with the given steps.
    Also stores all the available isochrones of different metallicities
    according to the ranges given to this parameters.
    '''

    # Call function to obtain the ranges and steps for the parameters along
    # with the path to the isochrone files and information about how they
    # are formatted (line_start).
    ranges_steps, indexes, iso_path, line_start = get_ranges_paths(sys_select)
    z_min, z_max = ranges_steps[0]
    age_min, age_max = ranges_steps[1]
    e_bv_min, e_bv_max, e_bv_step = ranges_steps[2]
    dis_mod_min, dis_mod_max, dis_mod_step = ranges_steps[3]
    
    # Read indexes for this Girardi output file.
    mini_indx, mag_indx, col_indx = indexes
    
    # Lists that store the colors, magnitudes and masses of the isochrone.
    isoch_list = []
    
    # List that will hold all the ages, metallicities and extinctions associated
    # with the stored isochrones.
    # isoch_params = [[metal, age, E(B-V), dis_mod], [], ...]
    isoch_params = []
    
    # Iterate through all metallicity files.
    for met_file in os.listdir(iso_path):
        
        # Process metallicity file only if it's inside the given range.
        if z_min<= met_file[:-4] <= z_max:
            
            # Store the metallicity value.
            metal = met_file[:-4]
            
            # Open the metallicity file.
            with open(join(iso_path, met_file), mode="r") as f_iso:
                
                # List that holds all the isochrone ages in the file.
                ages = []
                # Iterate through each line in the file.
                for line in f_iso:
                    
                    # Identify begginning of a defined isochrone.
                    if line.startswith(line_start):
                        # Store age value in 'ages' list if it falls inside
                        # the given range.
                        age_str = line.split("Age =")[1]
                        age = float(age_str[:-3])/1.e09
                        if age_min<= age <=age_max:
                            
                            # Save age in list.
                            ages.append(age)

                            # Save mag, color and mass values for each
                            # isochrone star.
                            if not line.startswith("#"):
                                reader = line.split()
                                # Color.
                                isoch_col.append(float(reader[col_indx]) -
                                float(reader[mag_indx]))
                                # Magnitude.
                                isoch_mag.append(float(reader[mag_indx]))
                                # Mass
                                isoch_mas.append(float(reader[mini_indx]))

                # Store lists in isochrone list.
                isoch_age = [isoch_col, isoch_mag, isoch_mas]

            # Store all isochrones in this metallicity file.


            # Store in list all the available isochrones to be used, according
            # to the ranges and steps given to its parameters metallicity, age,
            # extinction and distance modulus.
            for age_val in ages:
                
                # Loop through all extinction values.
                for e_bv in np.arange(e_bv_min, e_bv_max, e_bv_step):
                    # Loop through all distance modulus values.    
                    for dis_mod in np.arange(dis_mod_min, dis_mod_max,
                                             dis_mod_step):
                        # Store params for this isochrone.
                        isoch_params.append([metal, age_val, round(e_bv, 2),
                                             round(dis_mod, 2)])
                    
    return isoch_list, isoch_params, iso_path, line_start, indexes


#def read_isoch(m, a, e, d, sys_select, iso_path, line_start, indexes):
#    '''
#    Reads and stores an isochrone given the values of the parameters: m, a, e, d.
#    '''
#
#    # Lists that store the color, magnitude and masses of the isochrone.
#    iso_list, masses = [[], []], []
#    
#    # Read indexes for this Girardi output file.
#    mini_indx, mag_indx, col_indx = indexes
#    
#    # Read the file corresponding to the value of 'm' and store the isochrone
#    # of age 'a'.
#    m_file = m+'.dat'
#    with open(join(iso_path, m_file), mode="r") as f_iso:
#
#        # Set initial age value.
#        age = -99.
#        # Iterate through each line in the file.
#        for line in f_iso:
#            
#            age_flag = False
#            if line.startswith(line_start) and not age_flag:
#                # Save age value.
#                for i, part in enumerate(line.split("Age =")):
#                    # i=0 indicates the first part of that line, before 'Age ='
#                    if i==1:
#                        age = float(part[:-3])/1.e09
#                        age_flag = True
#            elif line.startswith(line_start) and age_flag:
#                break
#
#            if age == a:
#                # Save mag, color and mass values for each isochrone star.
#                if not line.startswith("#"):
#                    reader = line.split()
#                    # Magnitude.
#                    iso_list[1].append(float(reader[mag_indx]))
#                    # Color.
#                    iso_list[0].append(float(reader[col_indx]) -
#                    float(reader[mag_indx]))
#                    # Mass
#                    masses.append(float(reader[mini_indx]))
#
#        # Move isochrone according to E(B-V) and dis_mod values.
#        iso_color, iso_magnitude = move_track(iso_list, sys_select, e, d)
#        
#    # Append final data to array.
#    isochrone = [iso_color, iso_magnitude, masses]
#                    
#    return isochrone






def likelihood(synth_clust, obs_clust):
    '''
    Takes a synthetic cluster, compares it to the observed cluster and
    returns the weighted (log) likelihood value.
    This function follows the recipe given in Monteiro, Dias & Caetano (2010).
    '''
    
    # Store cluster data in new arrays: color & magnitude, their errors and 
    # each star's membership probabilities (weights)
    # Store color and magnitude into array.
#    data0 = zip(*obs_clust)
#    col_mag = np.array([data0[5], data0[3]], dtype=float)
#    # Store color and magnitude errors into array.
#    err_col_mag = np.array([data0[6], data0[4]], dtype=float)
#    print 'synth_clust', synth_clust, '\n'
   
    clust_stars_probs = []
    for indx,star in enumerate(obs_clust):
        # The first term does not depend on the synth stars.
        sigma_c, sigma_m = star[6], star[4]
        A = 1/(sigma_c*sigma_m)
        
        # Get probability for this cluster star.
        synth_stars = []
        for synth_st in zip(*synth_clust):
            # synth_st[0] = color ; synth_st[1] = mag
            B = np.exp(-0.5*((star[5]-synth_st[0])/sigma_c)**2)
            C = np.exp(-0.5*((star[3]-synth_st[1])/sigma_m)**2)
            synth_stars.append(A*B*C)
            
        # The final prob for this cluster star is the sum over all synthetic
        # stars.
        sum_synth_stars = sum(synth_stars) if sum(synth_stars)>0. else 1e-06
        clust_stars_probs.append(sum_synth_stars)
        
    # Store weights data (membership probabilities) into array.
    weights = np.array([zip(*obs_clust)[7]], dtype=float)   
    # Weight probabilities for each cluster star.
    weighted_probs = clust_stars_probs*weights
    
    # Get weighted likelihood.
#    L_x = reduce(lambda x, y: x*y, weighted_probs)
    
    # Final score: sum log likelihoods for each star in cluster.
    isoch_score = -sum(np.log(np.asarray(weighted_probs[0])))
    
    return isoch_score



def synthetic_clust(isochrone, mass_dist):
    '''
    Takes an isochrone and returns a synthetic cluster created according to
    a certain mass distribution.
    '''
    
    # Interpolate masses from the mass distribution into the isochrone to
    # obtain its magnitude and color values. Reject masses that fall outside of
    # the isochrone mass range.
    
    # Randomly select a fraction of these stars to be binaries.
    
    # Calculate the secondary masses of these binary stars.
    
    # Add masses and update array.
    
    # Randomly move stars according to given error distributions.
    
    # Remove stars according to a completness limit function.

    synth_clust = isochrone
    
    return synth_clust
    
    
    
def isoch_likelihood(m, a, e, d, sys_select, obs_clust, mass_dist):
    '''
    Call with given values for metallicity, age, extinction and distance modulus
    to generate a synthetic cluster with those parameters and compare it wiht
    the observed cluster.
    
    m, a, e, d = metallicity, age, extinction, distance modulus.
    '''
    
    # Store isochrone of metallicity value 'm' and age 'a' moved by the
    # values 'e' and 'd'.
    isoch_final = move_track(isoch_list, sys_select, e, d)
#    isoch_final = read_isoch(m, a, e, d, sys_select, iso_path, line_start,
#                             indexes)
    
    # Generate synthetic cluster using this "moved" isochrone and a mass
    # distribution.
    synth_clust = synthetic_clust(isoch_final, mass_dist)
    
    # Call function to obtain the likelihood by comparing the synthetic cluster
    # with the observed cluster.
    isoch_lik = likelihood(synth_clust, obs_clust)
    
    return isoch_lik
    


def boostrap_resample(memb_prob_avrg_sort):
    '''
    Resamples the observed cluster to use in the bootstrap process.
    '''
#    http://www.astroml.org/_modules/astroML/resample.html#bootstrap
    obs_clust = memb_prob_avrg_sort
    return obs_clust
    


def gip(sys_select, memb_prob_avrg_sort):
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
    
    # Store all posible combinations of metallicity, age, extinction and
    # distance modulus in isoch_params and all isochrones stored in all the
    # metallicity files in isoch_list. We do this so the files will only have
    # to be accessed once, thus being a more efficient method.
    
    # isoch_params = [params_1, ..., params_N]
    # params_i = [met_i, age_i, ext_i, dist_mod_i]
    
    # isoch_list = [isoch_1, ..., isoch_2]
    # isoch_i = [[colors], [magnitudes], [masses]]
    isoch_list, isoch_params, iso_path, line_start, indexes = get_isoch_params(sys_select)
    print isoch_params[0]
   
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

        # Call  algorithm to calculate the likelihoods for the set of
        # isochrones and return the best fitting parameters.
        
        if i==0:
            # Brute force.
            isoch_fit_params = brute_force(sys_select, isoch_params, iso_path,
                                           line_start, indexes, obs_clust,\
                                           mass_dist)
            # Genetic algorithm.
#            isoch_fit_params = genetic_algor(obs_clust, mass_dist)
        else:
            # Brute force.
            params_boot.append(brute_force(sys_select, isoch_params, iso_path,
                                           line_start, indexes, obs_clust,\
                                           mass_dist))
            # Genetic algorithm.
#            params_boot.append(genetic_algor(obs_clust, mass_dist))
        
    # Calculate errors for each parameter.
    isoch_fit_errors = np.mean(params_boot)
    
    return isoch_fit_params, isoch_fit_errors
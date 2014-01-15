# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 15:08:39 2014

@author: gabriel
"""

import os
from os.path import join

from get_mass_dist import mass_dist as md
import numpy as np


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
        mag_index = 9
         
    elif sys_select == '2':
        # Select cloud.
        cloud = raw_input('SMC or LMC cloud?')
        
        # Range of values where the parameters will move.
        e_bv_min, e_bv_max, e_bv_step = 0., 0.21, 0.01
        if cloud == 'SMC':
            dis_mod_min, dis_mod_max, dis_mod_step = 18.9, 18.91, 1
        elif cloud == 'LMC':
            dis_mod_min, dis_mod_max, dis_mod_step = 18.5, 18.51, 1
        z_min, z_max = 0.0005, 0.02
        age_min, age_max = 0.003, 12.6
    
        # Select Marigo or PARSEC tracks.        
        iso_select = raw_input('Select Marigo or PARSEC tracks as 1 or 2: ')
        if iso_select == '1':
            # Marigo.
            line_start = "#\tIsochrone\tZ = "
            mag_index = 9
            # Path where isochrone files are stored.
            iso_path = '/media/rest/github/isochrones/iso_wash_marigo'
        elif iso_select == '2':
            # PARSEC.
            line_start = "#\tIsochrone  Z = "
            mag_index = 10
            # Path where isochrone files are stored.
            iso_path = '/media/rest/github/isochrones/iso_wash_parsec'
    
    return e_bv_min, e_bv_max, e_bv_step, dis_mod_min, dis_mod_max,
    dis_mod_step, iso_path, line_start, mag_index
    
    

  
    
    
def get_isoch_params(sys_select):
    '''
    Reads and stores all available parameter values for the stored isochrones.
    '''

    # Call function to obtain the ranges for the parameters.
    e_bv_min, e_bv_max, e_bv_step, dis_mod_min, dis_mod_max, dis_mod_step,\
    iso_path, line_start, mag_index = get_ranges_paths(sys_select)

    # This list will hold every isochrone of a given age, metallicity moved
    # according to given values of distance modulues and extinction.
    # Each isochrrone is stored as a list composed of two sub-lists, one for
    # the mag and one for the color:
    # funcs = [[[mag values],[color values]], [[],[]], ...]
    isochrones =[]
    # List that will hold all the ages, metallicities and extinctions associated
    # with the isochrone in the same indexed position in 'funcs'.
    # params = [[metal, age, E(B-V), dis_mod], [], ...]
    isoch_params = []
    # Iterate through all isochrone files.
    for iso_file in os.listdir(iso_path):
    
        # Read the entire file once to count the total number of ages in it.    
        with open(join(iso_path, iso_file), mode="r") as f_iso:
            ages_total = 0
            for line in f_iso:
                if line.startswith(line_start):
                    ages_total += 1
        
        # Open the isochrone for this cluster's metallicity.
        with open(join(iso_path, iso_file), mode="r") as f_iso:
            
            # Store metallicity value.
            if iso_file[:-4] == 'zams':
                metal = 0.019
            else:
                metal = float(iso_file[:-4])
        
            # List that will hold all the isochrones of the same metallicity.
            iso_list = [[[], []] for _ in \
                        range(len(os.listdir(iso_path))*ages_total)]
    
            # List that holds the ages in the isochrone file.
            ages = []
            # This var counts the number of ages within a given metallicity file.
            age_indx = -1
            # Iterate through each line in the file.
            for line in f_iso:
                
                if line.startswith(line_start):
                    age_indx += 1
                    # Store age value in 'ages' list.
                    for i, part in enumerate(line.split("Age =")):
                        # i=0 indicates the first part of that line, before 'Age ='
                        if i==1:
                            ages.append(float(part[:-3])/1.e09)
    
                # Save mag and color values for each isochrone.
                if not line.startswith("#"):
                    reader = line.split()
                    # Magnitude.
                    iso_list[age_indx][1].append(float(reader[mag_index]))
                    # Color.
                    iso_list[age_indx][0].append(float(reader[8]) -
                    float(reader[mag_index]))
    
            # Separate isochrones with different ages.
            for age_val in range(age_indx+1):
                
                # Move isochrone according to E(B-V) and dis_mod values thus
                # generating a new function. Store this function in
                # 'isochrones' list.
                for e_bv in np.arange(e_bv_min, e_bv_max, e_bv_step):
                    
                    for dis_mod in np.arange(dis_mod_min, dis_mod_max,
                                             dis_mod_step):
                        
                        # Store params for this isochrone.
                        isoch_params.append([metal, ages[age_val],
                                             round(e_bv, 2), round(dis_mod, 2)])
                        
                        # Call function to move track.
                        iso_moved_0, iso_moved_1 = move_track(sys_select, \
                        dis_mod, e_bv, iso_list, age_val) 
                        # Append final x and y CMD data to array.
                        isochrones.append([iso_moved_0, iso_moved_1])
                    
    return isochrones, isoch_params



def read_isoch(m, a, e, d):
    '''
    Reads and stores an isochrone given the values of the parameters: m, a, e, d.
    '''

    # Call function to obtain these values.
    e_bv_min, e_bv_max, e_bv_step, dis_mod_min, dis_mod_max, dis_mod_step,\
    iso_path, line_start, mag_index = get_ranges_paths(sys_select)

    # This list will hold every isochrone of a given age, metallicity moved
    # according to given values of distance modulues and extinction.
    # Each isochrrone is stored as a list composed of two sub-lists, one for
    # the mag and one for the color:
    # funcs = [[[mag values],[color values]], [[],[]], ...]
    isochrones =[]
    # List that will hold all the ages, metallicities and extinctions associated
    # with the isochrone in the same indexed position in 'funcs'.
    # params = [[metal, age, E(B-V), dis_mod], [], ...]
    isoch_params = []
    # Iterate through all isochrone files.
    for iso_file in os.listdir(iso_path):
    
        # Read the entire file once to count the total number of ages in it.    
        with open(join(iso_path, iso_file), mode="r") as f_iso:
            ages_total = 0
            for line in f_iso:
                if line.startswith(line_start):
                    ages_total += 1
        
        # Open the isochrone for this cluster's metallicity.
        with open(join(iso_path, iso_file), mode="r") as f_iso:
            
            # Store metallicity value.
            if iso_file[:-4] == 'zams':
                metal = 0.019
            else:
                metal = float(iso_file[:-4])
        
            # List that will hold all the isochrones of the same metallicity.
            iso_list = [[[], []] for _ in \
                        range(len(os.listdir(iso_path))*ages_total)]
    
            # List that holds the ages in the isochrone file.
            ages = []
            # This var counts the number of ages within a given metallicity file.
            age_indx = -1
            # Iterate through each line in the file.
            for line in f_iso:
                
                if line.startswith(line_start):
                    age_indx += 1
                    # Store age value in 'ages' list.
                    for i, part in enumerate(line.split("Age =")):
                        # i=0 indicates the first part of that line, before 'Age ='
                        if i==1:
                            ages.append(float(part[:-3])/1.e09)
    
                # Save mag and color values for each isochrone.
                if not line.startswith("#"):
                    reader = line.split()
                    # Magnitude.
                    iso_list[age_indx][1].append(float(reader[mag_index]))
                    # Color.
                    iso_list[age_indx][0].append(float(reader[8]) -
                    float(reader[mag_index]))
    
            # Separate isochrones with different ages.
            for age_val in range(age_indx+1):
                
                # Move isochrone according to E(B-V) and dis_mod values thus
                # generating a new function. Store this function in
                # 'isochrones' list.
                for e_bv in np.arange(e_bv_min, e_bv_max, e_bv_step):
                    
                    for dis_mod in np.arange(dis_mod_min, dis_mod_max,
                                             dis_mod_step):
                        
                        # Store params for this isochrone.
                        isoch_params.append([metal, ages[age_val],
                                             round(e_bv, 2), round(dis_mod, 2)])
                        
                        # Call function to move track.
                        iso_moved_0, iso_moved_1 = move_track(sys_select, \
                        dis_mod, e_bv, iso_list, age_val) 
                        # Append final x and y CMD data to array.
                        isochrones.append([iso_moved_0, iso_moved_1])
                    
    return isochrone


def likelihood(synth_clust, obs_clust):
    '''
    Takes a synthetic cluster, compares it to the observed cluster and
    returns the weighted (log) likelihood value.
    This function follows the recipe given in Monteiro, Dias & Caetano (2010).
    '''
    
    # Store cluster data in new arrays: color & magnitude, their errors and 
    # each star's membership probabilities (weights)
    data0 = obs_clust
    # Store color and magnitude into array.
    col_mag = np.array([zip(*[data0[5], data0[3]])], dtype=float)
    # Store color and magnitude errors into array.
    err_col_mag = np.array([zip(*[data0[6], data0[4]])], dtype=float)
    # Store weights data (membership probabilities) into array.
    weights = np.array([data0[7]], dtype=float)        
    
    clust_stars_probs = []
    for indx,star in enumerate(col_mag):
        # The first term does not depend on the synth stars.
        sigma_c, sigma_m = err_col_mag[indx][0], err_col_mag[indx][1]
        A = 1/(sigma_c*sigma_m)
        
        # Get probability for this cluster star.
        sum_synth_stars = []
        for synth_st in synth_clust:
            # synth_st[0] = color ; synth_st[0] = mag
            B = np.exp(-0.5*((star[0]-synth_st[0])/sigma_c)**2)
            C = np.exp(-0.5*((star[1]-synth_st[1])/sigma_m)**2)
            sum_synth_stars.append(A*B*C)
            
        # The final prob for this cluster star is the sum over all synthetic
        # stars.
        clust_stars_probs.append(sum(sum_synth_stars))
        
    # Weight probabilities for each cluster star.
    weighted_probs = clust_stars_probs*weights
    
    # Get weighted likelihood.
    L_x = reduce(lambda x, y: x*y, weighted_probs)
    
    # Final score: sum log likelihoods for each star in cluster.
    isoch_score = sum(-np.log(L_x))
    
    return isoch_score



def get_synthetic_clust(isochrone, mass_dist):
    '''
    Takes an isochrone and returns a synthetic cluster created according to
    a certain mass distribution.
    '''

    
    return synth_clust
    
    
    
def move_track(sys_select, dis_mod, e_bv, iso_list, age_val):
    '''
    Recieves an isochrone of a given age and metallicity and modifies
    it according to given values for the extinction E(B-V) and distance
    modulus.
    '''
    iso_moved = [[], []]
    
    if sys_select == '1':
        # For UBVI system.
        #
        # E(B-V) = (B-V) - (B-V)o
        # Av = 3.1*E(B-V)
        # (mv - Mv)o = -5 + 5*log(d) + Av
        #
        Av = 3.1*e_bv
        for item in iso_list[age_val][1]:
            # mv affected by extinction.
            iso_moved[1].append(item + dis_mod + Av)
        for item in iso_list[age_val][0]:
            # (B-V) affected by extinction.
            iso_moved[0].append(item + e_bv)
    else:
        # For Washington system.
        #
        # E(C-T1) = 1.97*E(B-V) = (C-T1) - (C-T)o
        # M_T1 = T1 + 0.58*E(B-V) - (m-M)o - 3.2*E(B-V)
        #
        # (C-T1) = (C-T1)o + 1.97*E(B-V)
        # T1 = M_T1 - 0.58*E(B-V) + (m-M)o + 3.2*E(B-V)
        #
        V_Mv = dis_mod + 3.2*e_bv
        for item in iso_list[age_val][1]:
             # T1 magnitude affected by extinction.
            iso_moved[1].append(item - 0.58*e_bv + V_Mv)
        for item in iso_list[age_val][0]:
             # C-T1 color affected by extinction.
            iso_moved[0].append(item + 1.97*e_bv)    

    return iso_moved[0], iso_moved[1]    
    
    

def isoch_likelihood(m, a, e, d, obs_clust, mass_dist):
    '''
    Call with given values for metallicity, age, extinction and distance modulus
    to generate a synthetic cluster with those parameters and compare it wiht
    the observed cluster.
    
    m, a, e, d = metallicity, age, extinction, distance modulus.
    '''
    
    # Open file, given by the metallicity value 'm'.
    # Call function that reads all isochrones from stored files and creates new
    # extra ones according to the E(B-V) and dist_mod ranges given
    isochrones, isoch_params = read_isochrones(sys_select)   
    
    # Read isochrone, given by the age value 'a'.
    
 
    
    
    # Move isochrone in both axis according to 'e' and 'm'.
    isoch_final = move_track(sys_select, dis_mod, e_bv, iso_list, age_val)
    
    # Generate synthetic cluster using this "moved" isochrone and a mass
    # distribution.
    synth_clust = get_synthetic_clust(isoch_final, mass_dist)
    
    # Call function to obtain the likelihood by comparing the synthetic cluster
    # with the observed cluster.
    isoch_lik = likelihood(synth_clust, obs_clust)
    
    return isoch_lik
    


def boostrap_resample(memb_prob_avrg_sort):
    '''
    Resamples the observed cluster to use in the bootstrap process.
    '''
    obs_clust = memb_prob_avrg_sort
    return obs_clust
    


def brute_force(isoch_params, obs_clust, mass_dist):
    '''
    Brute force algorithm that computes the likelihoods for *all* the defined
    isochrones.
    '''
    # Initiate list that will hold the values (scores) which defines how well
    # each isochrone fits the observed data.
    score = []
    # Iterate through all the tracks defined and stored.
    for isoch in isoch_params:

        # Get parameters value from this isochrone.
        m, a, e, d = isoch
        
        # Call function that returns the score for a given track.
        isoch_score = isoch_likelihood(m, a, e, d, obs_clust, mass_dist)
        
        # Store the scores for each function/track into list.
        score.append(isoch_score)    
        
    # Find index of function with smallest value of the likelihoods.
    # This index thus points to the isochrone that best fits the observed
    # group of stars.
    best_func= np.argmin(score)
#    min_score = score[best_func]

    z_met, age_gyr, e_bv, dis_mod = [i for i in isoch_params[best_func]]
    dist_kpc = round(10**(0.2*(dis_mod+5.))/1000., 2)
    
    isoch_fit_params = [z_met, age_gyr, e_bv, dis_mod, dist_kpc]         
        
    return isoch_fit_params
    
    

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
    # distance modulus in a list.
    # isoch_params = [isoch_1, ..., isoch_N]
    # icosh_i = [met_i, age_i, ext_i, dist_mod_i]
    isoch_params = get_isoch_params(sys_select)
   
    # Begin bootstrap block.
    # NUmber of times to run the bootstrap block.
    N_B = 500.
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
            isoch_fit_params = brute_force(isoch_params, obs_clust, mass_dist)
            # Genetic algorithm.
#            isoch_fit_params = genetic_algor(obs_clust, mass_dist)
        else:
            # Brute force.
            params_boot.append(brute_force(isoch_params, obs_clust, mass_dist))
            # Genetic algorithm.
#            params_boot.append(genetic_algor(obs_clust, mass_dist))
        
    # Calculate errors for each parameter.
    isoch_fit_errors = np.mean(params_boot)
    
    return isoch_fit_params, isoch_fit_errors
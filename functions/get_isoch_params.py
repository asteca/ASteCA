# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 15:08:39 2014

@author: gabriel
"""

def read_isochrones():
    '''
    Reads and stores all isochrones available in the indicated folder.
    '''

    # This list will hold every isochrone of a given age, metallicity moved
    # according to given values of distance modulues and extinction.
    # Each isochrrone is stored as a list composed of two sub-lists, one for
    # the mag and one for the color:
    # funcs = [[[mag values],[color values]], [[],[]], ...]
    funcs =[]
    # List that will hold all the ages, metallicities and extinctions associated
    # with the isochrone in the same indexed position in 'funcs'.
    # params = [[metal, age, E(B-V), dis_mod], [], ...]
    params = []
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
                # generating a new function. Store this function in 'funcs' list.
                for e_bv in np.arange(e_bv_min, e_bv_max, e_bv_step):
                    
                    for dis_mod in np.arange(dis_mod_min, dis_mod_max, dis_mod_step):
                        
                        # Store params for this isochrone.
                        params.append([metal, ages[age_val], round(e_bv, 2),
                                       round(dis_mod, 2)])
                        
                        # Call function to move track.
                        iso_moved_0, iso_moved_1 = move_track(sys_select, dis_mod,
                                                              e_bv) 
                        # Append final x and y CMD data to array.
                        funcs.append([iso_moved_0, iso_moved_1])
                    
                    


def gip(memb_prob_avrg_sort):
    '''
    Perform a best fitting process to find the cluster's parameters:
    E(B-V), distance (distance modulus), age and metallicity.
    '''
    
    
    return isoch_fit_params
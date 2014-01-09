# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 10:24:45 2013

@author: gabriel
"""

import numpy as np

def mpas(cluster_region, runs_fields_probs, n_c, center_cl, clust_rad):
    """
    Use the number of approximate cluster's members n_c to identify the n_c
    most likely members among the stars inside the calculated cluster's radius,
    based on the probabilities assigned to each of these stars by the field
    decontamination algorithm.
    """
    
    # cluster_region = [[id,x,y,mag,e_mag,color1,e_col1], [], [], ...]
    # clus_reg_decont = [0.,0.1,0.5,0.,0.2,0.6,0.4,0.,0.2,0.5,0.3,0.9,...]
    
    # Every star in 'clus_reg_decont' contains the probability between 0 & 1
    # that the corresponding star in the list 'cluster_region' (located in the
    # index that corresponds to the same index) belongs to the cluster ie,
    # can be regarded as a cluster member.
    
    
    # Initiallize final list that contains all the lists with the different
    # realizations of the most probable members.
    membership_prob_list = []
    
    for reg_indx, clus_reg_decont_2 in enumerate(clus_reg_decont):
    
        # Generate a temporary list where each item is a list containing the
        # star's data from 'cluster_region' and its decontamination index
        # from 'clus_reg_decont'.
        temp_prob_members = []        
        
        # Create the new list adding the decontamination index from
        # 'clus_reg_decont' to the star's data from 'cluster_region'.
        for st_indx, star in enumerate(cluster_region):
            
            # Reject stars located outside the cluster's limit.
            dist = np.sqrt((center_cl[0]-star[1])**2 + \
            (center_cl[1]-star[2])**2)
            
            if dist <= clust_rad:
                temp_prob_members.append([reg_indx] + star + \
                [round(clus_reg_decont_2[st_indx],3)])
        
        # Stars inside the cluster's radius are now saved in the list 
        # 'temp_prob_members' where each item contains the data for each
        # star: 
        # [reg_indx,id,x,y,mag,e_mag,color1,e_col1,decont_index].
        
        # Append this list to the main list that stores all the list
        # for each field region.
        membership_prob_list.append(temp_prob_members)
        

    # Find average of decont_index value for stars in each sublist in the
    # main membership_prob_list list.
    arr = np.asarray(membership_prob_list)
    arr2 = arr[:, :, 1:].mean(0)
    # Round decont_index values to 3 decimals.
    membership_prob_avrg = np.apply_along_axis(lambda arr: \
    np.round(arr, 3), -1, arr2).tolist()


    # We now sort this list first by the decontamination index from max
    # value (1) to min (0) and then by its error values and magnitude value,
    # (in that order) from min to max value.
    # item[8] is the star's decont index and item[4] its magnitude.
    membership_prob_avrg_sort = sorted(membership_prob_avrg,
                                       key=lambda item: (-item[7], item[4],
                                                         item[6], item[3]))

    return membership_prob_avrg_sort
        
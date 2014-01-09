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
        
# DEPRECATED   
        # We now trim this list up to 'N_c' items so it will contain only the
        # brightest stars with the largest decontamination index as the
        # approximate number of cluster members.
#            most_prob_memb_nc = most_prob_memb_sort[:n_c]
# DEPRECATED   
        

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

    
    
# DEPRECATED
    # We now re-scale this list so that the minimum value will be set
    # as 0 and the maximum as 1.
    # Put prob values in first sublist of temp_lst.
#        temp_lst = list(reversed(zip(*most_prob_memb_nc)))
#        min_p, max_p = np.min(temp_lst[0]), np.max(temp_lst[0])
#        for indx, star in enumerate(most_prob_memb_nc):
#            # Call function to get re-scaled probability value.
#            prob_r = rescale(star[8], min_p, max_p)
#            # Update value in list.
#            most_prob_memb_nc[indx][8] = prob_r
# DEPRECATED

    
# DEPRECATED        
#        # Store all stars in all lists into a single list.
#        temp_single_list = []
#        for KDE_reg in most_prob_memb_list:
#            for star in KDE_reg:
#                temp_single_list.append(star)
#
#        # Iterate through all stars in the list averaging those that are found
#        # more then once.
#        temp_single_avrg = []
#        nums = defaultdict(list)
#        for arr in temp_single_list:
#            key = tuple(arr[1:8]) # make the x,y,mag,e_mag,col,e_col the keys
#            nums[key].append(arr[8]) # append the probability for the given key
#        
#        temp_single_avrg = [list(k) + [round(sum(vals)/len(vals),3)] for k,
#                            vals in nums.items()]
# DEPRECATED


    
# DEPRECATED        
    # Sort this list again as above (decontamination index first from max
    # value (1) to min (0) and then by its error values and magnitude value,
    # (in that order) from min to max value).
    # item[7] is the star's decont index and item[3] its magnitude.
#        most_prob_memb_avrg_temp = sorted(temp_single_avrg, key=lambda item: \
#                                    (-item[7], item[4], item[6], item[3]))
    
# DEPRECATED

    
# DEPRECATED
    # We now trim this list up to 'N_c' items so it will contain only the
    # brightest stars with the largest decontamination index as the
    # approximate number of cluster members.
#        most_prob_memb_avrg = most_prob_memb_avrg_temp[:n_c]
# DEPRECATED


    return membership_prob_avrg_sort, membership_prob_list
        
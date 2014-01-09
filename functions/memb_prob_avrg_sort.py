# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 10:24:45 2013

@author: gabriel
"""

import numpy as np

def mpas(cluster_region, runs_fields_probs, n_c, center_cl, clust_rad):
    """
    Generate a list with the averaged membership probabilities obtained
    with the decontamination algorithm for each star in cluster_region.
    Genarate a second list appending these averaged probabilities to each
    stars inside the cluster radius, sorted by the probabilities values.
    """
    
    # cluster_region = [[id,x,y,mag,e_mag,color1,e_col1], [], [], ...]
    
    # Average all Bayesian membership probabilities into a single value for
    # each star inside 'cluster_region'.
    clust_reg_prob_avrg = np.asarray(runs_fields_probs).mean(1).mean(0)
    
    # Create new list appending the membership probability to each star inside
    # the cluster radius.
    temp_prob_members = []
    for st_indx, star in enumerate(cluster_region):
        # Reject stars located outside the cluster's limit.
        dist = np.sqrt((center_cl[0]-star[1])**2 + (center_cl[1]-star[2])**2)
        
        if dist <= clust_rad:
            temp_prob_members.append(star + \
            [round(clust_reg_prob_avrg[st_indx],3)])
    
    # Stars inside the cluster's radius are now saved in the list 
    # 'temp_prob_members' where each item contains the data for each
    # star: [id,x,y,mag,e_mag,color1,e_col1,memb_prob].

    # Sort this list first by the membership probability from max
    # value (1) to min (0) and then by its error values and magnitude value,
    # (in that order) from min to max value.
    # item[7] is the star's memb_prob and item[3] its magnitude.
    membership_prob_avrg_sort = sorted(temp_prob_members,
                                       key=lambda item: (-item[7], item[4],
                                                         item[6], item[3]))

    return membership_prob_avrg_sort, clust_reg_prob_avrg
        
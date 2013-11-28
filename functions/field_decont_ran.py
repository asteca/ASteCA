# -*- coding: utf-8 -*-
"""
Created on Fri Sep  6 14:22:53 2013

@author: gabriel
"""


import numpy as np
import random as rd

# Field decontamination algorithm.

def field_decont_ran(flag_area_stronger, cluster_region, field_region, center_cl,
                 clust_rad):
    '''
    Apply random assignation decontamination algorithm.
    '''
    
    # Apply algorithm if l is bigger than 2, else skip it.
    if not(flag_area_stronger):
        
        print 'Applying random decontamination algorithm.'

        # cluster_region = [[id,x,y,mag,emag,col1,ecol1], [], [], ...]
        # len(cluster_region) = number of stars inside the cluster region
        # len(field_region[i]) = number of stars inside this field region

        # We assign to each star in the cluster region a random probability
        # between 0 and 1.
        
        # Initialize the list that will contain a number associated with each
        # star in the cluster region.
        clus_reg_decont = []
        
        # Initialize the list that will store the coordinates in CMD space and
        # the size of the box for each star in each of the field regions.
        field_reg_box = [[] for _ in range(len(field_region))]
                
        # Iterate through all stars in cluster region.
        for star in cluster_region:
            
            # Only compute if star is inside the cluster estimated radius.
            dist = np.sqrt((center_cl[0]-star[1])**2 + \
            (center_cl[1]-star[2])**2)
            if dist <= clust_rad[0]:
                ran_prob = rd.uniform(0., 1.)
                clus_reg_decont.append(ran_prob)
            else:
                clus_reg_decont.append(0.)

    # Skipping decontamination algorithm
    else:
        clus_reg_decont, field_reg_box = [], []

    return clus_reg_decont, field_reg_box, flag_area_stronger


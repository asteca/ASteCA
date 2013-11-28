"""
@author: gabriel
"""

import numpy as np

# Field decontamination algorithm: Variable Box (VB).

def field_decont_VB(flag_area_stronger, cluster_region, field_region):
    '''
    Apply modified version of Piatti et al. (2012) field decontamination
    algorithm.
    '''
    
    # Apply algorithm if l is bigger than 2, else skip it.
    if not(flag_area_stronger):
        
        print 'Applying VB decontamination algorithm.'

        # cluster_region = [[id,x,y,mag,emag,col1,ecol1], [], [], ...]
        # len(cluster_region) = number of stars inside the cluster region
        # len(field_region[i]) = number of stars inside this field region
        
        # We have the stars in the cluster region in the list
        # 'cluster_region' and the rest of the field stars regions defined in
        # lists inside the list 'field_region[]'.
        # The next step is to iterate through each star inside each the cluster
        # region and count the number of times said star is "assigned" to a star
        # in a field region.
        # When we finish, we have for each star in the cluster region a new
        # number attached which indicates how many times it was "removed" by the
        # algorithm. The bigger this number, the less likely it is that this
        # star is a true cluster member.
        
        # Initialize the list that will contain a number associated with each
        # star in the cluster region. This number means how many times that star
        # was removed by the decontamination algorithm. The bigger the number,
        # the more likely it is that star belongs to the field and not to the
        # cluster.
        clus_reg_decont = []
        clus_reg_decont2 = [0]*len(cluster_region)
        
        # Select box width scale.
        box_scale = 1.0
        
        # Initialize the list that will store the coordinates in CMD space and
        # the size of the box for each star in each of the field regions. We use
        # 'areas-1' to avoid counting the cluster region.
        field_reg_box = [[] for _ in range(len(field_region))]
        
        # Iterate through all the 'field stars' regions that were populated.
        for reg_index in range(len(field_region)):
            
            reg = field_region[reg_index]
            
            # Iterate through all stars inside this field region.
            for star1 in reg:

                # Initialize distance between star1 and the rest of the stars in
                # the region.
                dist_12 = [1000., 0.]
                # Initialize box values (just in case there are no stars in
                # the region)
                box_mag, box_col1 = 1., 0.5
                # Iterate again to find the closest star (star2) in the CMD
                # diagram, to star1 in the region.
                for star2 in reg:
                    
                    # Do not compare star1 to itself; [0] represents the star's
                    # ID and 3 its mag value.
                    if star2[0] == star1[0] and star2[3] == star1[3]:
                        # Star is the same, skip.
                        pass
                    else:
                        # star2[3] is the mag coordinate and star2[5] the col1
                        # coordinate
                        delta_mag, delta_col1 = abs(star1[3]-star2[3]), \
                        abs(star1[5]-star2[5])
                        # Get distance in the CMD space.
                        dist_12[1] = np.sqrt((delta_mag**2 + delta_col1**2))
                        # If this star is closer to star1 than the previous one
                        if dist_12[1] <= dist_12[0]:
                            # Update the new minimum distance.
                            dist_12[0] = dist_12[1]
                            # Store the box values for star1. Each 'box' value
                            # is a side of the rectangle centered in
                            # star1.
                            box_mag, box_col1 = box_scale*delta_mag,\
                                                box_scale*delta_col1
                
                # If box values are too big, reduce them to a default value.
                box_mag = min(box_mag, 1.)
                box_col1 = min(box_col1, 0.5)
                # If box values are too small, increase them.
                box_mag = max(box_mag, 0.2)
                box_col1 = max(box_col1, 0.1)
                    
                # Store the mag & col1 coordinates and the box sizes of the
                # star in a list.
                field_reg_box[reg_index].append([star1[3], star1[5], box_mag,
                                                box_col1])

                # Initialize distance value between star1 and cluster stars.
                dist_1c = [1000., 0.]
                
                # Initialize the indexes that tells me which star to remove from
                # the cluster region, if any.
                rem_star, rem_star_index = False, 0
                
                # We have now the box for star 1 as defined by box_mag and
                # box_col1. We iterate through all stars in the cluster region
                # to find the one inside this box (if one exists) located the
                # closest to its center.
                for cl_st_ind, clust_star in enumerate(cluster_region):
                    
                    # Check if this star in the cluster region falls inside the
                    # box defined for field star 'star1'.
                    delta_mag, delta_col1 = abs(star1[3]-clust_star[3]), \
                    abs(star1[5]-clust_star[5])
                    if delta_mag <= box_mag and delta_col1 <= box_col1:
                        # If it does, check the distance between the two stars
                        # in the CMD space.
                        dist_1c[1] = np.sqrt(delta_mag**2 + delta_col1**2)
                        if dist_1c[1] <= dist_1c[0]:
                            # If this distance is smaller than the previous one,
                            # update the new minimum distance.
                            dist_1c[0] = dist_1c[1]
                            # Store index pointing to that star as the one
                            # closest to star1 and set 'rem_star' as 'True' to
                            # know we have to remove a star.
                            rem_star, rem_star_index = True, cl_st_ind
                     
                # We have now the star in the cluster region closest to star1
                # and inside its box pointed by the index 'rem_star_index'.
                                                
                # Check if a star in the cluster region associated to star1 has
                # to be "removed".
                if rem_star:
                    # Increase the number of times this star was "removed".
                    clus_reg_decont2[rem_star_index] += 1
                    
            print '   Field %d processed.' % (reg_index)
            
        # Normalize and invert the 'cluster_reg_decont' list so the stars with
        # the bigger numbers (closer to 1) are the ones with the higher
        # probability of being a cluster star.
        # Check if list is empty.
        if not clus_reg_decont2:
            # If list is empty
            max_num = 0
        else:
            max_num = float(max(clus_reg_decont2))
        if max_num == 0:
        # If the list is empty.
            flag_area_stronger = True
            clus_reg_decont2, field_reg_box = [], []           
        else:
            clus_reg_decont2 = [round((max_num-float(i)), 2) for i in\
            clus_reg_decont2]


        clus_reg_decont.append(clus_reg_decont2)
    # Skipping decontamination algorithm
    else:
        clus_reg_decont, field_reg_box = [], []

    return clus_reg_decont, field_reg_box, flag_area_stronger

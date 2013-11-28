# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 10:06:31 2013

@author: gabriel
"""

from functions.get_spiral import spiral as gs
import numpy as np


def get_regions(x_center_bin, y_center_bin, width_bins, histo, clust_rad,
                h_manual, stars_in, stars_out):
    '''
    Get cluster and field regions around the cluster's center, used by the
    decontamination algorithms(s).
    '''
    
    # Use the bin center obtained with the smallest bin width
    x_c_b, y_c_b = x_center_bin[0], y_center_bin[0]
    
    # Define 'cluster region' as a spiral centered around the cluster
    # and of area a bit bigger than that defined by the cluster's radius.
    # The 'cluster region' will have an area as close as possible to
    # (l*r)^2, composed of bins of the minimum width used, where 'l' represents
    # the number of times the side of the 'square' that contains the cluster
    # region is bigger than the cluster's radius. The
    # 'field_regions' will have an area equal to that of the 'cluster_region'.

    # Get area as total number of bins in 2D hist times the area of each bin.
    area = len(histo[0][0])*len(histo[0])*(width_bins[0]**2)

    # We set the number of areas used: one is for the cluster region and the
    # rest for the field regions.
    areas = 21

    # Calculate 'length' = constant that multiplied by the cluster's radius gives
    # the length of the side of the square, called 'cluster_region', that
    # contains the cluster.
    # We set 'length' so that all the regions of area
    # (l*r)^2 sum up to less than 50% of area: 50*area/100 > 5*(l*r)^2. This
    # means that length = sqrt(0.1*area)/clust_rad[0]. First we see if length
    # obtained this way is equal or greater than 3. A value less than 3 would
    # mean that the area around the cluster is comparable to the cluster's area.
    # If this is the case we increase the 50% value by 10% and try again; we
    # repeat this process until a minimum value of 3 is obtained or until we
    # reach 100% for the value of the frame's area. If we encounter this
    # last scenario, we set length = 2 and raise a flag.
    
    for percent in range(50, 101, 10):
        length = int(np.sqrt(percent*area/(areas*100.))/clust_rad[0])
        if length < 3:
            # Increase the percentage and try again.
            pass
        else:
            # Exit the for loop.
            break

    # If length is still smaller than 3 we try with less field areas and
    # increase the threshold to accept a value of 2.
    if length < 3.:
        for areas_2 in range(areas, 1, -1):           
            length = int(np.sqrt(percent*area/(areas_2*100.))/clust_rad[0])
            if length >= 2:
                break
        
    # If after this length is still smaller than 2, this means that the area
    # taken by the cluster is bigger than half of the total frame. In this case
    # the algorithm can not be applied as is since we have no equal sized area
    # with field stars to compare it with the cluster region. We skip the
    # decontamination algorithm alltoghether and raise a flag.
    flag_area_stronger = False
    if length < 2:
        flag_area_stronger = True

    # Apply algorithm if length is bigger than 2, else skip it.
    if not(flag_area_stronger):

        # Calculate number of bins such that their combined area is
        # approximately (l*r)^2. See that: num_bins_area * width_bins[0]^2 = 
        # (l*r)^2.
        num_bins_area = int((length*clust_rad[0]/width_bins[0])**2)
    
        # This list holds all the lists that will hold the bins that
        # compose the areas for the cluster (first list) and the remaining
        # field regions.
        regions = [[[],[]] for _ in range(areas)]
    
        # We'll populate all regions with 'num_bins_area' bins each, starting
        # from the central bin for the 'cluster_region' and when that area is
        # filled then we move on to the next area (ie: first 'field region') and
        # so on. The bins are located in a spiral ring-like fashion around the
        # center bin.
        
        # Get list that contains the spiral as a list of x,y coordinates (also
        # stored as lists) starting from the initial bin [0, 0].
        spiral = gs()
        # Initialize the bin counter that tells me how many bins are already
        # added to a given region.
        bin_count = [0 for _ in range(areas)]
    
        # We add the bin to the corresponging region until 'num_bins_area'
        # are added to that region, then we move on to the next region.
        # Since the initial bin in the spiral corresponds to the center
        # of the cluster, the first region to be populated will be the
        # 'cluster_region'.
        reg_index = 0
    
        # Iterate through all the bins in the spiral, starting from
        # the central bin which corresponds to the center of the cluster
        for sp_item in spiral:
            
            
            # THIS IS AN ADDITION MADE SO THAT THE KDE ALGORITHM USES CMD's
            # OF AREAS EQUAL TO THE CLUSTER AREA FOR THE FIELDS SINCE ONLY STARS
            # INSIDE THE CLUSTER RADIUS ARE USED TO OBTAIN THE CLUSTER'S CMD.
            if reg_index == 0:
                # We're inside the cluster region: use big square area
               num_bins_area_2 = num_bins_area
            else:
                # We're inside a field region: use cluster's area.
                num_bins_area_2 = int(np.pi*((clust_rad[0]/width_bins[0])**2))
    
    
            # If the region is not filled yet, we add this bin to this region.
            if bin_count[reg_index] <= num_bins_area_2:
                # Check if the bin exists in the 2D histogram.
                try:
                    histo[0][x_c_b+sp_item[0]][y_c_b+sp_item[1]]
                except IndexError:
                    # Item out of histogram range. Do nothing.
                    pass
                else:
                    # Check that the index is not negative because python
                    # will assign items in lists even if they are pointed
                    # as negative values, ie: list1[1][-2]; which in this
                    # case makes no sense because it would fall out of the
                    # 2D histogram.
                    if (x_c_b+sp_item[0])>=0 and (y_c_b+sp_item[1])>=0:
                        # If the item exists, we store the coordinates of
                        # that bin in this region in both coordinates.
                        regions[reg_index][0].append(x_c_b+sp_item[0])
                        regions[reg_index][1].append(y_c_b+sp_item[1])
                        # And we increase the bin count.
                        bin_count[reg_index] += 1
            # If this region is filled, then we move on to the next one.
            else:
                reg_index += 1
    
            # If all regions have been filled, we exit the spiral for loop.
            if reg_index == areas:
                break


        # At this point we have the list 'regions' containing the cluster region
        # and a given number of field stars region. Each list
        # inside 'regions' is composed of two lists, the first one containing
        # the x bin coordinates for every bin that corresponds to that region
        # and the second list the y bin coordinates. We need to obtain the stars
        # located within each of those areas so as to be able to construct the
        # CMD diagrams of the cluster and the field stars regions. To do this
        # we use 'h_manual' wich is a list that already contains the stars that
        # fall inside each of the bins in the 2D histogram along with their
        # relevant data (ID, x, y, T1, etc..)
        
        
        # Initialize the lists that will hold the bins that compose the
        # areas for the cluster and all the field regions.
        cluster_region, field_region = [], [[] for _ in range(reg_index-1)]
        
        # Iterate through all the regions that were populated (not necessarily
        # all the regions defined if the cluster area is big or the frame too
        # small).
        for index in range(reg_index):
            # Iterate through all the x values of the bins in this region.
            for xbin_index, xbin in enumerate(regions[index][0]):
                ybin = regions[index][1][xbin_index]
                
                if index == 0:
                    # We are inside the 'cluster_region'.
                    # 'h_manual[xbin][ybin][0]' is the number of stars inside
                    # this bin.
                    for star in range(h_manual[xbin][ybin][0]):
                        # Add all the stars inside this bin to the cluster area.
                        # We use 'star+1' to skip the first item which holds the
                        # number of stars in the bin.
                        # Check to see if star belongs to the group of stars
                        # that were NOT removed because of its large error.
                        # h_manual[xbin][ybin][star+1][0] is the ID of the star.
                        st_id = h_manual[xbin][ybin][star+1][0]
                        st_x = h_manual[xbin][ybin][star+1][1]
                        if any(st_id == i[0] for i in stars_in) and \
                        any(st_x == i[1] for i in stars_in) or \
                        any(st_id == i[0] for i in stars_out) and \
                        any(st_x == i[1] for i in stars_out):
                            cluster_region.append(h_manual[xbin][ybin][star+1])
                else:
                    # We are inside one of the field stars regions.
                    for star in range(h_manual[xbin][ybin][0]):
                        st_id = h_manual[xbin][ybin][star+1][0]
                        st_x = h_manual[xbin][ybin][star+1][1]
                        if any(st_id == i[0] for i in stars_in) and \
                        any(st_x == i[1] for i in stars_in) or \
                        any(st_id == i[0] for i in stars_out) and \
                        any(st_x == i[1] for i in stars_out):
                            field_region[index-1].append(h_manual[xbin]\
                            [ybin][star+1])

                            
        # Now we have the stars inside the cluster region in the list
        # 'cluster_region' and the rest of the field stars regions defined in
        # lists inside the list 'field_region[]'.
        
        # If any of the field regions has less than 4 stars then we remove it
        # from the list otherwise the KDE or the p-value algorithms will fail.
        # This lists stores the indexes of the empty regions.
        field_regs_del = []
        for indx, s_lst in enumerate(field_region):
            if len(s_lst) < 4:
                field_regs_del.append(indx)
        # Delete empty regions this way to avoid messing with the indexes.
        for index in sorted(field_regs_del, reverse=True):
            del field_region[index]

        # If after removing the empty regions no regions are left, raise the
        # flag.
        if not(field_region):
            flag_area_stronger = True
            
    # Skipping regions separation.
    else:
        cluster_region, field_region = [], []
        
    return flag_area_stronger, cluster_region, field_region
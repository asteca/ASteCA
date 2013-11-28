# -*- coding: utf-8 -*-
"""
Created on Sun Sep  8 09:25:39 2013

@author: gabriel
"""

import numpy as np
from scipy.ndimage.filters import gaussian_filter

# Field decontamination algorithm.

def field_decont_dias(flag_area_stronger, cluster_region, center_cl, clust_rad):
    
    # Apply algorithm of l is bigger than 2, else skip it.
    if not(flag_area_stronger):
        
        print 'Applying Dias decontamination algorithm.'
        
        # cluster_region = [[id,x,y,T1,eT1,CT1,eCT1], [], [], ...]
        # len(cluster_region) = number of stars inside the cluster region
        # len(field_region[i]) = number of stars inside this field region
        
        # We have the stars in the cluster region in the list
        # 'cluster_region' and the rest of the field stars regions defined in
        # lists inside the list 'field_region[]'.
        
        # We iterate through each star inside the cluster
        # region and obtain the probability 'prob_f' of that star of belonging
        # to it by applying Dias et al 2012 membership likelihood.
        
        # Initialize 'clus_reg_decont' empty list. The numbers inside this
        # list are the probabilities obtained for that star in the cluster
        # region. The bigger this number, the more likely it is that star is a
        # true cluster member.
        clus_reg_decont = []
        

        # Obtain width of Gaussian kernel to use.
#        dist_list = []
#        for star in cluster_region:
#            for star2 in cluster_region:
#                if star != star2:
#                    dist = np.sqrt((star2[1]-star[1])**2 + (star2[2]-star[2])**2)
#                    dist_list.append(dist)
#        st_dv = np.std(dist_list)
        st_dv = 140.
        
        # Calculate density values for the cluster region.
        x_data = list(zip(*cluster_region)[1])
        y_data = list(zip(*cluster_region)[2])        
        
        xmin, xmax = min(x_data), max(x_data)
        ymin, ymax = min(y_data), max(y_data)
        rang = [[xmin, xmax], [ymin, ymax]]
        
        # Set bin width to be used in the histogram.
        d_b = 1.
        
        binsxy = [int((xmax - xmin) / d_b), int((ymax - ymin) / d_b)]
        
        # hist is the 2D histogran, xedges & yedges store the edges of the bins
        hist, xedges, yedges = np.histogram2d(x_data, y_data, range=rang, 
                                              bins=binsxy)

        # h_g is the 2D histogram with a gaussian filter applied
        h_g = gaussian_filter(hist, st_dv, mode='constant')

        # x_cent_g & y_cent_g store the x,y bin coordinates of the bin with the
        # maximum value in the 2D filtered histogram, ie: the center of the 
        # putative cluster.
        x_cent_g, y_cent_g = np.unravel_index(h_g.argmax(), h_g.shape)
        
        # Obtain maximum value of density in cluster region.
        max_dens = h_g[x_cent_g][y_cent_g]

#        # Calculate centers in pixel coordinates.
#        center_x_g = np.average(xedges[x_cent_g:x_cent_g + 2])
#        center_y_g = np.average(yedges[y_cent_g:y_cent_g + 2])
 
        # Get probability of belonging to the cluster for every star in the
        # cluster region.
        for star in cluster_region:

#            dist = np.sqrt((center_cl[0]-star[1])**2 + \
#            (center_cl[1]-star[2])**2)
            
            # Convert pixel to bin coordinates.
            x_bc = int((star[1]-min(x_data))/d_b)
            y_bc = int((star[2]-min(y_data))/d_b)
            try: # Check to see if the bin coords exist.
                h_g[x_bc][y_bc]
            except IndexError:
            # Bin doesn't exist so the star is located at the border.
                # Assign density value of bin 2 values below.
                dens = h_g[x_bc-2][y_bc-2]
            else:
            # Bins exist so star is not at the border.
                # Get density value for the position.
                dens = h_g[x_bc][y_bc]
            
            # Store all likelihoods in array.
            likhd_sum = []
            for star2 in cluster_region:
                
                # Skip same star.
                if star2 != star:
                
                    dist = np.sqrt((star2[1]-star[1])**2 + (star2[2]-star[2])**2)
                    
                    # Calculate likelihood.
                    
                    # C filter value
                    c_col = star[5]+ star[3]
                    c_err = abs(star[6] - star[4])
                    c_col2 = star2[5]+ star2[3]
                    
#                    np.exp(-0.5*((star[5]-star2[5])/star[6])**2)*\

                    likhd = np.exp(-0.5*((star[3]-star2[3])/star[4])**2)*\
                    np.exp(-0.5*((c_col-c_col2)/c_err)**2)*\
                    np.exp(-0.5*((clust_rad[0]-dist)/13.)**2)*\
                    np.exp(-0.5*((dens-max_dens)/st_dv)**2)*\
                    (1./(star[4]*star[6]*13.*st_dv))
                    
                    # Store value in list.
                    likhd_sum.append(likhd)
                
            # Get likelihood value for 'star' in cluster region.
            prob = sum(likhd_sum)
            clus_reg_decont.append(round(prob, 2))

        print '   Cluster region processed.'
        
    # Skipping decontamination algorithm
    else:
        clus_reg_decont = []
        
    return clus_reg_decont

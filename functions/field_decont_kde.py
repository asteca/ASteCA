"""
@author: gabriel
"""

import numpy as np
import random as rd
from scipy import stats

'''
KDE field decontamination algorithm.
'''


def gauss_error(col_lst, e_col_lst, mag_lst, e_mag_lst):
    # Randomly move mag and color through a Gaussian function.
    col_gauss = rd.gauss(np.array(col_lst), np.array(e_col_lst))
    mag_gauss = rd.gauss(np.array(mag_lst), np.array(e_mag_lst))
    
    return col_gauss, mag_gauss


def bw_val(runs, run_num, num_stars):
    # Calculate bandwidth value based on the number of stars in the region
    # and the run_num value.
    bw_min = 0.1
    bw_max = num_stars**(-1./6.) # Scott's rule
    bw_step = (bw_max-bw_min)/(runs-1)
    bw_range = np.arange(bw_min, bw_max+(bw_step/2.), bw_step)
    # Value to return    
    bw_go = bw_range[run_num]
    return bw_go
    


def field_decont_kde(flag_area_stronger, cluster_region, field_region,
                     col1_data, mag_data, center_cl, clust_rad):
    
    # Apply algorithm of l is bigger than 2, else skip it.
    if not(flag_area_stronger):
        
        # Set number of runs for the KDE algorithm.
        runs = 5
        
        # Set the number of samples used in the Monte Carlo itegration.
        mc_sample = 10
        
        print 'Applying KDE decontamination algorithm.'
        
        # cluster_region = [[id,x,y,T1,eT1,CT1,eCT1], [], [], ...]
        # len(cluster_region) = number of stars inside the cluster region
        # len(field_region[i]) = number of stars inside this field region
        
        # We have the stars in the cluster region in the list
        # 'cluster_region' and the rest of the field stars regions defined in
        # lists inside the list 'field_region[]'.
        
        # We iterate through each star inside the cluster
        # region and obtain the probability 'prob_f' of that star of belonging
        # to a given field region. The value '1-prob_f' will be the probability
        # of that star of being a cluster star. We repeat this for each field
        # region and average the probabilities obtained. Then we get the
        # probability for that star in the cluster region of belonging to
        # the cluster region sequence using the same method and average this
        # value with the one above. This last step prevents stars in crowded
        # regions from being assigned very low probabilities and also helps
        # to lower the probabilities of stars in the cluster region located
        # distant from other stars in the cluster CMD.
        
        # Initialize 'clus_reg_decont' empty list. The numbers inside this
        # list are the average of the probabilities obtained for each star in
        # the cluster region in each field region (inverted by taking the 1-*
        # substraction) and the probability in the cluster region. The
        # bigger this number, the more likely it is that that star is a true
        # cluster member.
        clus_reg_decont = []
            
        # Obtain max and min values for the axis. We do it this way so that
        # the plots of the KDEs will be aligned with the scatters plots
        # of the field region later on. Note that the limits are independent
        # of mag_lst and col_lst (used below).
        xmin, xmax = max(-0.9, min(col1_data)-0.2),\
                             min(3.9, max(col1_data)+0.2)
        ymin, ymax = max(mag_data)+0.5, min(mag_data)-0.5        
   
        
        # Run 'runs' times.
        flag_25, flag_50, flag_75 = False, False, False
        for run_num in range(runs):

            # Format data used to obtain kernel estimate in cluster region. Only
            # use stars inside cluster's radius.
            col_lst, e_col_lst, mag_lst, e_mag_lst = [], [], [], []
            for star in cluster_region:
                
                dist = np.sqrt((center_cl[0]-star[1])**2 + \
                (center_cl[1]-star[2])**2)
                
                if dist <= clust_rad[0]: 
                    # Color data.
                    col_lst.append(star[5])
                    # Color error.
                    e_col_lst.append(star[6])
                    # Magnitude data.
                    mag_lst.append(star[3])
                    # Magnitude error.
                    e_mag_lst.append(star[4]) 

            # Move magnitude and colors for cluster stars randomly using a
            # Gaussian function.
            col_gauss, mag_gauss = gauss_error(col_lst, e_col_lst, mag_lst,
                                               e_mag_lst)
            
            # Obtain the KDE for the cluster region (stars inside cluster's
            # radius).
            x, y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
            positions = np.vstack([x.ravel(), y.ravel()])
            values = np.vstack([col_gauss, mag_gauss])
            # Call function to calculate bandwidth to use.
            bw_choice = bw_val(runs, run_num, len(col_gauss))
            # The results are HEAVILY dependant on the bandwidth used here.
            # See: http://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.gaussian_kde.html
            kernel_cl = stats.gaussian_kde(values, bw_method = bw_choice)
            # Store the KDE for plotting it later on.
            if run_num == 4:
                kde_cl = np.reshape(kernel_cl(positions).T, x.shape)
            
            # Initialize 'clus_reg_decont_cl' empty list which will hold the
            # probs obtained for each star with the cluster's KDE.
            clus_reg_decont_cl = []
            
            # Get probability of belonging to the cluster for every star in the
            # cluster region.
            for star in cluster_region:
                # Only compute if star is inside the cluster estimated radius.
                # We do this to save processing time.
                dist = np.sqrt((center_cl[0]-star[1])**2 + \
                (center_cl[1]-star[2])**2)
                
                if dist <= clust_rad[0]:            
                    # Compute the point below which to integrate.
                    iso = kernel_cl((star[5], star[3]))
                    
                    # Random sample from the KDE distribution.
                    sample = kernel_cl.resample(size= mc_sample)
                    
                    # Filter the sample.
                    insample = kernel_cl(sample) < iso
                    
                    # As per Monte Carlo, the integral is equivalent to the
                    # probability of drawing a point that gets through the
                    # filter.
                    integral = insample.sum() / float(insample.shape[0])
                    # This is the probability associated with that point of
                    # belonging to the cluster sequence.
                    clus_reg_decont_cl.append(integral)
                else:
                    clus_reg_decont_cl.append(0.)
                
                # Cluster KDE obtained.
            
            
            # This list will hold the KDEs used for plotting.
            if run_num == 4:
                kde = [[] for _ in range(len(field_region))]
            # This list will hold the probabilities for each field region.
            clus_reg_decont_fl = [[] for _ in range(len(field_region))]
            
            # Iterate through all the 'field stars' regions that were populated.
            for reg_index in range(len(field_region)):
                
                reg = field_region[reg_index]
                
                # Format field region data.
                col_lst, e_col_lst, mag_lst, e_mag_lst = [], [], [], []
                for star in reg:
                    # Color data.
                    col_lst.append(star[5])
                    # Color error.
                    e_col_lst.append(star[6])
                    # Magnitude data.
                    mag_lst.append(star[3])
                    # Magnitude error.
                    e_mag_lst.append(star[4])
                    
                # Move magnitude and colors randomly using a Gaussian function.
                col_gauss, mag_gauss = gauss_error(col_lst, e_col_lst, mag_lst,
                                                   e_mag_lst)
                
                # Obtain the KDE for this field region.            
                x, y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
                positions = np.vstack([x.ravel(), y.ravel()])
                values = np.vstack([col_gauss, mag_gauss])
                # Call function to calculate bandwidth to use.
                bw_choice = bw_val(runs, run_num, len(col_gauss))
                # The results are HEAVILY dependant on the bandwidth used here.
                # See: http://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.gaussian_kde.html
                kernel = stats.gaussian_kde(values, bw_method = bw_choice)
                # Store the KDE for plotting it later on.
                if run_num == 4 and reg_index == 0:
                    kde_f = np.reshape(kernel(positions).T, x.shape)
                
                # We iterate through all stars in the cluster region to obtain
                # the probability for each one of belonging to this field region.
                for star in cluster_region:
                    # Only compute if star is inside the cluster estimated radius.
                    # We do this to save sime processing time.
                    dist = np.sqrt((center_cl[0]-star[1])**2 + \
                    (center_cl[1]-star[2])**2)
                    if dist <= clust_rad[0]:
    
                        # Compute the point below which to integrate.
                        iso = kernel((star[5], star[3]))
                        
                        # Sample from the KDE distribution
                        sample = kernel.resample(size= mc_sample)
                        
                        # Filter the sample
                        insample = kernel(sample) < iso
                        
                        # The integral is equivalent to the probability of
                        # drawing a point that gets through the filter
                        integral = insample.sum() / float(insample.shape[0])
                        # Get the probability associated with that star of
                        # belonging to the cluster as 1 minus the probability of
                        # belonging to the field.
                        prob = 1 - integral

                        # Save probab value for this star of belonging to the
                        # cluster.
                        clus_reg_decont_fl[reg_index].append(prob)
                    else:
                        clus_reg_decont_fl[reg_index].append(0.)
                    
            # Field processed.
    
            # Now we have the probabilities of each star of belonging to the
            # cluster sequence obtained with the cluster KDE in
            # 'clus_reg_decont_cl' and the probs obtained with the field stars
            # regions in 'clus_reg_decont_fl'
            
            # First we average all the values in the field stars regions lists.
            # Change to array because its faster.
            arr_clus_reg_decont_fl = np.array(clus_reg_decont_fl)
            # Create new list with averaged values.
            avrg_field_reg = np.mean(arr_clus_reg_decont_fl, 0)
            
            # Now we average that list with the one obtained using the cluster
            # field.
            # Create new array of both array and list.
            data = np.array([avrg_field_reg, clus_reg_decont_cl])
            # Store average in temporary array.
            temp = np.average(data, axis=0)
            # Convert to final list.
            clus_reg_decont.append(temp.tolist())
            
            if run_num+1 >= runs/4 and flag_25 == False:
                print '  25% done'
                flag_25 = True
            elif run_num+1 >= runs/2 and flag_50 == False:
                print '  50% done'
                flag_50 = True
            elif run_num+1 >= (runs/2 + runs/4) and flag_75 == False:
                print '  75% done'
                flag_75 = True
            elif run_num+1 == runs:
                print '  100% done'
                
        
    # Skipping decontamination algorithm
    else:
        clus_reg_decont, kde_cl, kde_f = [], [], []
    
    return clus_reg_decont, kde_cl, kde_f

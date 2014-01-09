"""
@author: gabriel
"""

import numpy as np
import random as rd
from scipy import stats

'''
Bayesian KDE field decontamination algorithm.
'''


def gauss_error(col_lst, e_col_lst, mag_lst, e_mag_lst):
    # Randomly move mag and color through a Gaussian function.
    col_gauss = rd.gauss(np.array(col_lst), np.array(e_col_lst))
    mag_gauss = rd.gauss(np.array(mag_lst), np.array(e_mag_lst))
    
    return col_gauss, mag_gauss


def bw_val(runs, run_num, num_stars):
    '''
    Calculate bandwidth value based on the number of stars in the region
    and the run_num value.
    '''
    bw_min = 0.1
    bw_max = num_stars**(-1./6.) # Scott's rule
    bw_step = (bw_max-bw_min)/(runs-1)
    bw_range = np.arange(bw_min, bw_max+(bw_step/2.), bw_step)
    # Value to return    
    bw_go = bw_range[run_num]
    return bw_go
    
    
def mc_probability(reg, xmin, xmax, ymin, ymax, runs, run_num, cluster_region,
                   center_cl, clust_rad, mc_sample):
    '''
    Calculate probability/likelihood for each cluster region star through
    Monte Carlo integration.
    '''   
    # Format region data.
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
    
    # Obtain the KDE for this region.            
    x, y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    positions = np.vstack([x.ravel(), y.ravel()])
    values = np.vstack([col_gauss, mag_gauss])
    # Call function to calculate bandwidth to use.
    bw_choice = bw_val(runs, run_num, len(col_gauss))
    # The results are HEAVILY dependant on the bandwidth used here.
    # See: http://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.gaussian_kde.html
    kernel = stats.gaussian_kde(values, bw_method = bw_choice)
    
    reg_decont = []
    # We iterate through all stars in the cluster region to obtain
    # the probability/likelihood for each one of belonging to this region.
    for star in cluster_region:
        # Only compute if star is inside the cluster estimated radius.
        dist = np.sqrt((center_cl[0]-star[1])**2 + (center_cl[1]-star[2])**2)
        if dist <= clust_rad[0]:

            # Compute the value below which to integrate.
            iso = kernel((star[5], star[3]))
            
            # Sample from the KDE distribution
            sample = kernel.resample(size=mc_sample)
            
            # Filter the sample
            insample = kernel(sample) < iso
            
            # The integral is equivalent to the probability of
            # drawing a point that gets through the filter
            integral = insample.sum() / float(insample.shape[0])
            # Avoid 'nan' and/or 'infinite' solutions.
            integral = integral if integral > 0. else 0.000001
            
            # Save probability value for this star of belonging to this region.
            reg_decont.append(integral)
        else:
            reg_decont.append(0.000001)
    
    return reg_decont, kernel, positions, x



def field_decont_kde(cluster_region, field_region, col1_data, mag_data,
                     center_cl, clust_rad):
    '''
    Main function.
    '''
    
    # Set total number of runs for the KDE algorithm to 100.
    runs = int(100/len(field_region))
    
    # Set the number of samples used by the Monte Carlo integration.
    mc_sample = 10
    
    print 'Applying KDE decontamination algorithm.'
    
    # cluster_region = [[id,x,y,T1,eT1,CT1,eCT1], [], [], ...]
    # len(cluster_region) = number of stars inside the cluster region
    # len(field_region[i]) = number of stars inside this field region
    # Stars in the cluster region in the list 'cluster_region' and the
    # rest of the field stars regions defined in lists inside the list
    # 'field_region[]'.
    
    # Obtain max and min values for the axis. We do it this way so that
    # the plots of the KDEs will be aligned with the scatters plots
    # of the field region later on. Note that the limits are independent
    # of mag_lst and col_lst (used below) and depend only on the full
    # ranges of magnitude and color for the data.
    xmin, xmax = max(-0.9, min(col1_data)-0.2),\
                         min(3.9, max(col1_data)+0.2)
    ymin, ymax = max(mag_data)+0.5, min(mag_data)-0.5  
    
    # Initialize 'clus_reg_decont' empty list. The numbers inside this
    # list are the average of the probabilities obtained for each star in
    # the cluster region in each field region (inverted by taking the 1-*
    # substraction) and the probability in the cluster region. The
    # bigger this number, the more likely it is that that star is a true
    # cluster member.
    clus_reg_decont = []  
    
    # Run 'runs' times.
    flag_25, flag_50, flag_75 = False, False, False
    for run_num in range(runs):
        
        # This list will hold the probabilities for each field region.
        clus_reg_decont_fl = [[] for _ in field_region]
        # Iterate through all the 'field stars' regions that were populated.
        for indx, fl_region in enumerate(field_region):
            
            # Obtain likelihoods for each star in the cluster region
            # using this field region, ie: P(A)
            reg_decont_fl, kernel, positions, x = \
            mc_probability(fl_region, xmin, xmax, ymin, \
            ymax, runs, run_num, cluster_region, center_cl, clust_rad, \
            mc_sample)
            # Store number of stars in field region.
            n_fl = len(fl_region)
            # Store the KDE for plotting it later on.
            if run_num == 4 and indx == 0:
                kde_f = np.reshape(kernel(positions).T, x.shape)

            # Randomly shuffle the stars in the cluster region.
            clust_reg_shuffle = np.random.permutation(cluster_region)
            # Remove n_fl random stars from the cluster region and pass
            # it to the function that obtains the likelihoods for each
            # star in the "cleaned" cluster region, ie: P(B)
            if n_fl < len(cluster_region):
                clust_reg_clean = clust_reg_shuffle[n_fl:]
            else:
                # If field region has more stars than the cluster region,
                # don't remove any star. This should not happen though.
                clust_reg_clean = clust_reg_shuffle
            n_cl = len(clust_reg_clean)
            reg_decont_cl, kernel, positions, x = \
            mc_probability(clust_reg_clean, xmin, xmax, \
            ymin, ymax, runs, run_num, cluster_region, center_cl, \
            clust_rad, mc_sample)
            # Cluster KDE obtained. 
            # Store the KDE for plotting it later on.
            if run_num == 4:
                kde_cl = np.reshape(kernel(positions).T, x.shape)

            # Obtain Bayesian probability for each star in the cluster
            # region.
            p_a, p_b = np.array(reg_decont_fl), np.array(reg_decont_cl)
            bayes_prob = 1./(1. + (n_fl*p_a)/(n_cl*p_b))
            # Store probabilities obtained with this field region.
            clus_reg_decont_fl[indx] = bayes_prob

        # Now we have the probabilities of each star of belonging to the
        # cluster sequence in 'clus_reg_decont_fl'
        
        # Average all the probability values for each star.
        avrg_field_reg = np.mean(np.array(clus_reg_decont_fl), 0)
        # Store in final list as list.
        clus_reg_decont.append(avrg_field_reg.tolist())
        
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
                
    return clus_reg_decont, kde_cl, kde_f

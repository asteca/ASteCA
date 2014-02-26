"""
@author: gabriel
"""

import os.path
import numpy as np
from scipy import stats


def gauss_error(col, e_col, mag, e_mag):
    # Randomly move mag and color through a Gaussian function.
    col_gauss = col + np.random.normal(0, 1, len(col))*e_col
    mag_gauss = mag + np.random.normal(0, 1, len(col))*e_mag
    
    return col_gauss, mag_gauss


def bw_val(runs, run_num, num_stars):
    '''
    Calculate bandwidth value based on the number of stars in the region
    and the run_num value.
    '''
    bw_min = 0.1
    bw_max = num_stars**(-1./6.) # Scott's rule
    if runs > 1:
        bw_step = (bw_max-bw_min)/(runs-1)
        bw_range = np.arange(bw_min, bw_max+(bw_step/2.), bw_step)
        # Value to return    
        bw_go = bw_range[run_num]
    else:
        bw_go = bw_max
    return bw_go
    
    
def mc_probability(reg_call, reg, xmin, xmax, ymin, ymax, runs, run_num,
                   cluster_region, center_cl, clust_rad, mc_sample):
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
        if dist <= clust_rad:

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
            if reg_call == 'field':
                reg_decont.append(1.)
            elif reg_call == 'clust':
                reg_decont.append(0.000001)
    
    return reg_decont, kernel, positions, x



def field_decont_kde(flag_area_stronger, cluster_region, field_region,
                     col1_data, mag_data, center_cl, clust_rad, clust_name,
                     sub_dir, da_params):
    '''
    Bayesian KDE field decontamination algorithm.
    '''

    mode, run_n, mc_samp, mypath2 = da_params

    # Check if at least one field region was obtained.
    if flag_area_stronger:
        print "WARNING: no field regions found. Using 'skip'."
        mode = 'skip'
    
    # Check if 'mode' was correctly set, else use 'skip'.
    if mode not in ['auto', 'manual', 'read', 'skip']:
        print "WARNING: Wrong name for 'mode' in input file. Using 'skip'."
        mode = 'skip'
        
    if mode == 'read':
        # Check if file exists.
        memb_file = mypath2+'/'+sub_dir+'/'+clust_name+'_memb.dat'
        if not os.path.isfile(memb_file):
            # File does not exist.
            print "WARNING: members file does not exist. Using 'skip'."
            mode = 'skip'
        
    # Run algorithm for any of these selections.
    if mode == 'auto' or mode == 'manual':

        if mode == 'auto':
            # Set total number of runs for the KDE algorithm to 100.
            runs = int(100/len(field_region))
            # Set the number of samples used by the Monte Carlo integration.
            mc_sample = 1000
        elif mode == 'manual':
            # Take values from input data file.
            runs, mc_sample = run_n, mc_samp        
        
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
        
        # This list holds one sub-list per run. Each of those sub-lists holds N
        # sub-sub-lists (where N = len(field_region)) with the membership
        # probabilities assigned to each star in the 'cluster_region'.
        #
        # runs_fields_probs = [A_1, A_2, ..., A_runs]
        # A_i = [B_1, B_2, ..., B_N] ; N = len(field_regions)
        # B_j = [p_1, p_2, ..., p_n] ; n = len(cluster_region)
        # A,B --> arrays ; p --> floats (Bayesian probabilities)
        runs_fields_probs = []
        
        # Run 'runs' times.
        flag_25, flag_50, flag_75 = False, False, False
        for run_num in range(runs):
            
            # This list will hold the probabilities for each field region.
            field_reg_probs = [[] for _ in field_region]
            # Iterate through all the 'field stars' regions that were populated.
            for indx, fl_region in enumerate(field_region):
                
                # Obtain likelihoods for each star in the cluster region
                # using this field region, ie: P(A)
                reg_call = 'field'
                reg_decont_fl, kernel, positions, x = \
                mc_probability(reg_call, fl_region, xmin, xmax, ymin, \
                ymax, runs, run_num, cluster_region, center_cl, clust_rad, \
                mc_sample)
                # Store number of stars in field region.
                n_fl = len(fl_region)
                # For plotting purposes.
#                if run_num == 4 and indx == 0:
#                    kde_f = np.reshape(kernel(positions).T, x.shape)
    
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
                reg_call = 'clust'
                reg_decont_cl, kernel, positions, x = \
                mc_probability(reg_call, clust_reg_clean, xmin, xmax, \
                ymin, ymax, runs, run_num, cluster_region, center_cl, \
                clust_rad, mc_sample)
                # Cluster KDE obtained. 
                # For plotting purposes.
#                if run_num == 4:
#                    kde_cl = np.reshape(kernel(positions).T, x.shape)
    
                # Obtain Bayesian probability for each star in cluster_region.
                p_a, p_b = np.array(reg_decont_fl), np.array(reg_decont_cl)
                bayes_prob = 1./(1. + (n_fl*p_a)/(n_cl*p_b))
                # Store probabilities obtained with this field region.
                field_reg_probs[indx] = bayes_prob
    
            # Now we have the probabilities of each star in 'cluster_region' of
            # being an actual cluster member (membership) in 'field_reg_probs',
            # one sub-list per field region.
            
            # Append this list to the list that holds all the runs.
            runs_fields_probs.append(field_reg_probs)
            
            if run_num+1 >= 0.25 * runs and flag_25 is False:
                print '  25% done'
                flag_25 = True
            elif run_num+1 >= 0.5 * runs and flag_50 is False:
                print '  50% done'
                flag_50 = True
            elif run_num+1 >= 0.75 * runs and flag_75 is False:
                print '  75% done'
                flag_75 = True
            elif run_num+1 == runs:
                print '  100% done'            

    elif mode == 'read':
        print 'Reading membership probabilities from file.'
        # File where membership probabilities are stored.
        memb_file = mypath2+'/'+sub_dir+'/'+clust_name+'_memb.dat'
        # Read probabilities from file.           
        data = np.loadtxt(memb_file, unpack=True)
        id_list = data[0].tolist()
        
        probs = []
        # Assign probabilities read from file according to the star's IDs.
        # Those stars outside the cluster's radius are assigned a very low value.
        for indx,star in enumerate(cluster_region):
            if star[0] in id_list:
                # Index of star in file.
                i = id_list.index(star[0])
                # Assign the probability stored in file for this star.
                probs.append(data[7][i])
            else:
                probs.append(0.000001)
            
        # Store probabilities in list.
        runs_fields_probs = [[probs]]
        
    elif mode == 'skip':
        print 'Assigning equal probabilities to all stars inside cluster radius.'
        # Assign equal probabilities to all stars.
        runs_fields_probs = [[[1.]*len(cluster_region)]]
              
#    return runs_fields_probs, kde_cl, kde_f
    return runs_fields_probs

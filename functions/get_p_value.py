# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 17:58:28 2013

@author: gabriel
"""

import numpy as np
import random as rd
from scipy.stats import ks_2samp
from scipy import stats

import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
ks = importr('ks')
kde_test = ks.kde_test
hpi_kfe = ks.Hpi_kfe

'''
Compare the cluster region KDE with all the field region KDEs using Duong's
ks package (developed in R) to obtain a p-value. This value will be close
to 1 if the cluster region is very similar to the field regions and
closer to 0 as it differentiates from it.
As a rule of thumb, a p-value > 0.05 (ie: 5%) indicates that the overdensity
is most likely not a real cluster.
'''


def gauss_error(col_lst, e_col_lst, mag_lst, e_mag_lst):
    # Randomly move mag and color through a Gaussian function.
    col_gauss = rd.gauss(np.array(col_lst), np.array(e_col_lst))
    mag_gauss = rd.gauss(np.array(mag_lst), np.array(e_mag_lst))
    
    return col_gauss, mag_gauss



def get_CMD(region):
    
    # Obtain CMD for given region.
    # Format field region data.
    col_lst, e_col_lst, mag_lst, e_mag_lst = [], [], [], []
    for star in region:
        # Color data.
        col_lst.append(star[5])
        # Color error.
        e_col_lst.append(star[6])
        # Magnitude data.
        mag_lst.append(star[3])
        # Magnitude error.
        e_mag_lst.append(star[4])
        
    # Move magnitude and colors randomly according to their errors,
    # using a Gaussian function.
    col_gauss, mag_gauss = gauss_error(col_lst, e_col_lst, mag_lst,
                                       e_mag_lst)

    matrix_0 = [col_gauss, mag_gauss]
    # Format values so the list will be composed of [star1, star2, ...] where
    # star1 = [color, magnitude]
    matrix_1 = map(list, zip(*matrix_0))
    # Put all stars into a single list.
    matrix = [star for axis in matrix_1 for star in axis]
    return matrix



def get_pval(flag_area_stronger, cluster_region, field_region,
                     col1_data, mag_data, center_cl, clust_rad):
    
    # Apply algorithm if enough field regions around the cluster region
    # were selected.
    if not(flag_area_stronger):
        
        print 'Obtaining p_value for cluster region vs field regions.'
        
        # Set number of runs for the p_value algorithm with a maximum of
        # 100 if only one field region was used.
#        runs = int(100/len(field_region))
        print 'runs:', int(100/len(field_region))
        runs = 5
        
        # Only use stars inside cluster's radius.
        cluster_region_r = []
        for star in cluster_region:
            dist = np.sqrt((center_cl[0]-star[1])**2 + \
            (center_cl[1]-star[2])**2)
            if dist <= clust_rad[0]: 
                cluster_region_r.append(star)
                
        # The first list holds all the p_values obtained comparing the cluster
        # region with the field regions, the second one holds p_values for field
        # vs field comparisions.
        p_vals_cl, p_vals_f = [], []
        # Iterate a given number of times.
        for run_num in range(runs):
            # Loop through all the field regions.
            for indx, f_region in enumerate(field_region):
                
                # CMD for cluster region.
                matrix_cl = get_CMD(cluster_region_r)
                rows_cl = int(len(matrix_cl)/2)
                # CMD for 1st field region.
                matrix_f1 = get_CMD(f_region)
                rows_f1 = int(len(matrix_f1)/2)
                
                # Create matrices for these CMDs.
                m_cl = robjects.r.matrix(robjects.FloatVector(matrix_cl),
                                       nrow=rows_cl, byrow=True)
                m_f1 = robjects.r.matrix(robjects.FloatVector(matrix_f1),
                                         nrow=rows_f1, byrow=True)
                                         
                # Bandwith matrices.
                hpic = hpi_kfe(x=m_cl, binned=True)
                hpif1 = hpi_kfe(x=m_f1, binned=True)
                
                # Call 'ks' function to obtain p_value.
                # Cluster vs field p_value.
                res_cl = kde_test(x1=m_cl, x2=m_f1, H1=hpic, H2=hpif1)
                p_val_cl = res_cl.rx2('pvalue')
                # Store cluster vs field p-value.
                p_vals_cl.append(float(str(p_val_cl)[4:]))

                # Compare the field region used above with all the remaining
                # field regions.
                for f_region2 in field_region[indx:]:
                
                    # CMD for 2nd field region.
                    matrix_f2 = get_CMD(f_region2)
                    rows_f2 = int(len(matrix_f2)/2)
                    # Matrix.
                    m_f2 = robjects.r.matrix(robjects.FloatVector(matrix_f2),
                                             nrow=rows_f2, byrow=True)
                    # Bandwith.
                    hpif2 = hpi_kfe(x=m_f2, binned=True)
            
                    # Field vs field p_value.
                    res_f = kde_test(x1=m_f1, x2=m_f2, H1=hpif1, H2=hpif2)
                    p_val_f = res_f.rx2('pvalue')
                    # Store field vs field p-value.
                    p_vals_f.append(float(str(p_val_f)[4:]))


            if run_num == runs/4:
                print '  25% done'
            elif run_num == runs/2:
                print '  50% done'
            elif run_num == (runs/2 + runs/4):
                print '  75% done'
            
        # Obtain average for all p-values.
        p_val_cl_avrg = round(np.average(p_vals_cl), 2)
        p_val_f_avrg = round(np.average(p_vals_f), 2)
        
        # Obtain final p_value for the *distributions* of p_values obtained
        # above.
        m_cl = robjects.FloatVector(p_vals_cl)
        m_f = robjects.FloatVector(p_vals_f)
        res_cl_f = kde_test(x1=m_cl, x2=m_f)
        p_val_cl_f = res_cl_f.rx2('pvalue')
        
        kol_smir_d_p = ks_2samp(p_vals_cl, p_vals_f)
        kol_smir_d, kol_smir_p = kol_smir_d_p[0], kol_smir_d_p[1]
        print '  K-S statistic', kol_smir_d
        print '  K-S p-value', kol_smir_p
        
        p_value = round(float(str(p_val_cl_f)[4:]), 2)
        print '  KDE Test p_value', p_value
        
#        print '\n', p_vals_cl, '\n'
#        print p_vals_f, '\n'
#        print p_value
#        raw_input()

    # Skipping decontamination algorithm
    else:
        p_value, p_vals_cl, p_vals_f, kde_cl_norm, kde_f_norm, p_val_cl_avrg,\
        p_val_f_avrg = -1., [], [], [], [], -1., -1.
        

    # For plotting purposes only.        
        
    # Obtain the 1D KDE for the cluster region (stars inside cluster's
    # radius) vs field regions.
    xmin, xmax = min(p_vals_cl), max(p_vals_cl)
    x_cl = np.mgrid[xmin:xmax:100j]
    kernel_cl = stats.gaussian_kde(p_vals_cl)
    kde_cl = np.reshape(kernel_cl(x_cl).T, x_cl.shape)
    # 1-D KDE for cluster vs field region p_values.
    kde_cl_norm = (np.array(kde_cl) - min(kde_cl))/max(kde_cl)

    # Obtain the 1D KDE for the field regions vs field regions.
    xmin, xmax = min(p_vals_f), max(p_vals_f)
    x_f = np.mgrid[xmin:xmax:100j]
    kernel_f = stats.gaussian_kde(p_vals_f)
    kde_f = np.reshape(kernel_f(x_f).T, x_f.shape)
    # 1-D KDE for field vs field region p_values.
    kde_f_norm = (np.array(kde_f) - min(kde_f))/max(kde_f)
 
    return p_value, p_vals_cl, p_vals_f, kde_cl_norm, kde_f_norm, p_val_cl_avrg,\
    p_val_f_avrg

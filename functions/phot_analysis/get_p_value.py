# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 17:58:28 2013

@author: gabriel
"""

import numpy as np
from scipy import stats
from scipy.integrate import quad
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr

# Define variables to communicate with package 'R'.
ks = importr('ks')
kde_test = ks.kde_test
hpi_kfe = ks.Hpi_kfe


def gauss_error(col, e_col, mag, e_mag):
    # Randomly move mag and color through a Gaussian function.
    col_gauss = col + np.random.normal(0, 1, len(col)) * e_col
    mag_gauss = mag + np.random.normal(0, 1, len(col)) * e_mag

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


def get_pval(cl_region, field_region, col1_data, mag_data, pv_params,
        flag_area_stronger):
    '''
    Compare the cluster region KDE with all the field region KDEs using Duong's
    ks package (developed in R) to obtain a p-value. This value will be close
    to 1 if the cluster region is very similar to the field regions and
    closer to 0 as it differentiates from it.
    As a rule of thumb, a p-value > 0.05 (ie: 5%) indicates that one should
    reject the null hypothesis that the KDEs arose from ther same distribution.
    We assign a probability of the overdensity being a real cluster as 1
    minus the overlap between the KDEs of the distributions of p-values for
    the cluster vs field and field vs field comparisons.
    '''

    mode_pv, num_runs = pv_params
    flag_pval_test = True if mode_pv in {'auto', 'manual'} else False

    # skip test if < 10 members are found within the cluster's radius.
    flag_few_members = False if len(cl_region) > 10 else True
    if flag_pval_test and flag_few_members:
        print '  WARNING: < 10 stars in cluster region, skipping test.'
        flag_pval_test = False

    # Check if test is to be applied or skipped. Check if field regions
    # where found.
    if flag_pval_test and not flag_area_stronger:

        print 'Obtaining p_value for cluster region vs field regions.'

        # Set number of runs for the p_value algorithm with a maximum of
        # 100 if only one field region was used.
        runs = int(100 / len(field_region)) if mode_pv == 'auto' else num_runs

        # The first list holds all the p_values obtained comparing the cluster
        # region with the field regions, the second one holds p_values for
        # field vs field comparisions.
        p_vals_cl, p_vals_f = [], []
        # Iterate a given number of times.
        milestones = [25, 50, 75, 100]
        for run_num in range(runs):
            # Loop through all the field regions.
            for indx, f_region in enumerate(field_region):

                # CMD for cluster region.
                matrix_cl = get_CMD(cl_region)
                rows_cl = int(len(matrix_cl) / 2)
                # CMD for 1st field region.
                matrix_f1 = get_CMD(f_region)
                rows_f1 = int(len(matrix_f1) / 2)

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
                p_vals_cl.append(float(str(p_val_cl)[4:].replace(',', '.')))

                # Compare the field region used above with all the remaining
                # field regions. This results in [N*(N+1)/2] combinations of
                # field vs field comparisions.
                for f_region2 in field_region[indx:]:

                    # CMD for 2nd field region.
                    matrix_f2 = get_CMD(f_region2)
                    rows_f2 = int(len(matrix_f2) / 2)
                    # Matrix.
                    m_f2 = robjects.r.matrix(robjects.FloatVector(matrix_f2),
                                             nrow=rows_f2, byrow=True)
                    # Bandwith.
                    hpif2 = hpi_kfe(x=m_f2, binned=True)

                    # Field vs field p_value.
                    res_f = kde_test(x1=m_f1, x2=m_f2, H1=hpif1, H2=hpif2)
                    p_val_f = res_f.rx2('pvalue')
                    # Store field vs field p-value.
                    p_vals_f.append(float(str(p_val_f)[4:].replace(',', '.')))

            percentage_complete = (100.0 * (run_num + 1) / runs)
            while len(milestones) > 0 and percentage_complete >= milestones[0]:
                print "  {}% done".format(milestones[0])
                # Remove that milestone from the list.
                milestones = milestones[1:]

        # For plotting purposes.

        # Define KDE limits.
        xmin, xmax = -1., 2.
        x_kde = np.mgrid[xmin:xmax:1000j]

        # Obtain the 1D KDE for the cluster region (stars inside cluster's
        # radius) vs field regions.
        kernel_cl = stats.gaussian_kde(p_vals_cl)
        # KDE for plotting.
        kde_cl_1d = np.reshape(kernel_cl(x_kde).T, x_kde.shape)

        # Obtain the 1D KDE for the field regions vs field regions.
        kernel_f = stats.gaussian_kde(p_vals_f)
        # KDE for plotting.
        kde_f_1d = np.reshape(kernel_f(x_kde).T, x_kde.shape)

        # Calculate overlap between the two KDEs.
        def y_pts(pt):
            y_pt = min(kernel_cl(pt), kernel_f(pt))
            return y_pt

        overlap = quad(y_pts, -1., 2.)
        # Store y values for plotting the overlap filled.
        y_over = [float(y_pts(x_pt)) for x_pt in x_kde]

        # Probability value for the cluster.
        prob_cl_kde = 1 - overlap[0]

        # Store all return params in a single list.
        pval_test_params = [prob_cl_kde, p_vals_cl, p_vals_f, kde_cl_1d,
            kde_f_1d, x_kde, y_over]

        print 'Probability of physical cluster obtained ({:.2f}).'.format(
            prob_cl_kde)

    # Skip process.
    else:
        print 'Skipping p-value test for cluster.'
        # Pass empty lists to make_plots.
        pval_test_params = [-1., [], [], [], [], [], []]

    return pval_test_params, flag_pval_test

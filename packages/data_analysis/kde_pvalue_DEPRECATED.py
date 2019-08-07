
from .. import update_progress
import numpy as np
from scipy import stats
from scipy.integrate import quad


def gauss_error(col, e_col, mag, e_mag):
    # Randomly move mag and color through a Gaussian function.
    col_gauss = col + np.random.normal(0, 1, len(col)) * e_col
    mag_gauss = mag + np.random.normal(0, 1, len(col)) * e_mag

    return col_gauss, mag_gauss


def get_CMD(region):
    # TODO generalize to N dimensions

    # Obtain CMD for given region.
    # Format field region data.
    col_lst, e_col_lst, mag_lst, e_mag_lst = [], [], [], []
    for star in region:
        # Color data.
        col_lst.append(star[5][0])
        # Color error.
        e_col_lst.append(star[6][0])
        # Magnitude data.
        mag_lst.append(star[3][0])
        # Magnitude error.
        e_mag_lst.append(star[4][0])

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


def KDE_test(clp, pvalue_runs):
    """
    Compare the cluster region KDE with all the field region KDEs using Duong's
    R-ks package to obtain a p-value. This value will be close to 1 if the
    cluster region is very similar to the field regions and closer to 0 as it
    differentiates from it.

    As a rule of thumb, a p-value > 0.05 (ie: 5%) indicates that one should
    reject the null hypothesis that the KDEs arose from the same distribution.

    We assign a probability of the overdensity being a real cluster as 1
    minus the overlap between the KDEs of the distributions of p-values for
    the cluster vs field and field vs field comparisons.
    """
    # mags, cols = cld['mags'], cld['cols']
    flag_pval_test = True if pvalue_runs > 0 else False

    # Skip test if < 10 members are found within the cluster's radius.
    flag_few_members = False if len(clp['cl_region_c']) > 10 else True

    # Check if test is to be applied or skipped. Check if field regions
    # where found.
    if not flag_pval_test:
        print('Skipping KDE p-value test for cluster.')
        # Pass empty lists and re-write flag.
        flag_pval_test, pval_test_params = False, [-1.]

    elif flag_pval_test and clp['flag_no_fl_regs_c']:
        print('No field regions. Skipping KDE p-value test for cluster.')
        # Pass empty lists and re-write flag.
        flag_pval_test, pval_test_params = False, [-1.]

    elif flag_pval_test and flag_few_members:
        print('  WARNING: < 10 stars in cluster region.'
              ' Skipping KDE p-value test.')
        # Pass empty lists and re-write flag.
        flag_pval_test, pval_test_params = False, [-1.]

    # Run process.
    elif flag_pval_test:

        # Define variables to communicate with package 'R'.
        import rpy2.robjects as robjects
        from rpy2.robjects.packages import importr
        ks = importr('ks')
        kde_test = ks.kde_test
        hpi_kfe = ks.Hpi_kfe

        print('Obtaining KDE p-value for cluster vs field regions.')

        if len(clp['field_regions_c']) == 1 and pvalue_runs == 1:
            pvalue_runs = 2
            print("  WARNING: only 1 field region defined. Function\n"
                  "  will run 2 times.")

        # The first list holds all the p_values obtained comparing the cluster
        # region with the field regions, the second one holds p_values for
        # field vs field comparisons.
        p_vals_cl, p_vals_f = [], []
        # Iterate a given number of times.
        run_total, runs = int(pvalue_runs * len(clp['field_regions_c'])), 0
        for run_num in range(pvalue_runs):
            # Loop through all the field regions.
            for indx, f_region in enumerate(clp['field_regions_c']):

                # CMD for cluster region.
                matrix_cl = get_CMD(clp['cl_region_c'])
                rows_cl = int(len(matrix_cl) / 2)
                # CMD for 1st field region.
                matrix_f1 = get_CMD(f_region)
                rows_f1 = int(len(matrix_f1) / 2)

                # Create matrices for these CMDs.
                m_cl = robjects.r.matrix(robjects.FloatVector(matrix_cl),
                                         nrow=rows_cl, byrow=True)
                m_f1 = robjects.r.matrix(robjects.FloatVector(matrix_f1),
                                         nrow=rows_f1, byrow=True)
                # Bandwidth matrices.
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
                # field vs field comparisons.
                for f_region2 in clp['field_regions_c'][(indx + 1):]:

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

                runs += 1
                update_progress.updt(run_total, runs)

        # For plotting purposes.

        # Define KDE limits.
        xmin, xmax = -1., 2.
        x_kde = np.mgrid[xmin:xmax:1000j]

        # Obtain the 1D KDE for the cluster region (stars inside cluster's
        # radius) vs field regions.
        kernel_cl = stats.gaussian_kde(p_vals_cl)
        # KDE for plotting.
        kde_cl_1d = np.reshape(kernel_cl(x_kde).T, x_kde.shape)

        # Check if field regions were compared among each other.
        if p_vals_f:
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
        else:
            # If not, assign probability as 1 minus the average of the cluster
            # vs the single field region used.
            prob_cl_kde = 1. - np.mean(p_vals_cl)
            # Pass empty lists for plotting.
            kde_f_1d, y_over = np.asarray([]), []

        # Store all parameters in a single list.
        pval_test_params = [prob_cl_kde, kde_cl_1d, kde_f_1d, x_kde, y_over]

        print('Probability of physical cluster obtained ({:.2f}).'.format(
            prob_cl_kde))

    return pval_test_params, flag_pval_test


def main(clp, pvalue_runs, R_in_place, **kwargs):
    """
    Only run function if the necessary application and packages are in place.
    """
    if R_in_place:
        pval_test_params, flag_pval_test = KDE_test(clp, pvalue_runs)
    else:
        print("  WARNING: missing package, skipping KDE p-value test.")
        pval_test_params, flag_pval_test = [np.nan], False

    clp['pval_test_params'], clp['flag_pval_test'] =\
        pval_test_params, flag_pval_test
    return clp

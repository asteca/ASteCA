# -*- coding: utf-8 -*-

from os.path import join, realpath, dirname, exists, isfile
from os import getcwd, mkdir, rmdir, makedirs
import time
import shutil
import gc  # Garbage collector.

# Import files with defined functions
from functions.read_paths import read_paths as rp
from functions.get_in_params import get_in_params as gip
from functions.create_out_data_file import create_out_data_file as c_o_d_f
from functions.get_data_semi import get_semi as g_s
from functions.get_phot_data import get_data as gd
from functions.trim_frame import trim_frame as t_f
from functions.get_center import get_center as g_c
from functions.manual_histo import manual_histo as mh
from functions.get_background import get_background as gbg
from functions.get_dens_prof import get_dens_prof as gdp
from functions.get_radius import get_clust_rad as gcr
from functions.get_king_prof import get_king_profile as gkp
from functions.err_accpt_rejct import err_accpt_rejct as ear
from functions.get_in_out import get_in_out as gio
from functions.get_integ_mag import integ_mag as g_i_m
from functions.get_members_number import get_memb_num as g_m_n
from functions.get_cont_index import cont_indx as g_c_i
from functions.get_regions import get_regions as g_r
from functions.field_decont_bys import field_decont_bys as fdb
from functions.get_p_value import get_pval as g_pv
from functions.get_qqplot import qqplot as g_qq
from functions.get_completeness import mag_completeness as m_c
from functions.get_isoch_params import ip
from functions.reduce_membership import red_memb as rm
from functions.best_fit_synth_cl import best_fit as bfsc
from functions.make_plots import make_plots as mp
from functions.add_data_output import add_data_output as a_d_o
from functions.cl_members_file import cluster_members_file as c_m_f


def ocaat_main(f_indx, sub_dir, out_file_name, gip_params):
    '''
    This is the main container which holds the calls to all functions for
    a given photometric data file.
    '''

    # Start timing this loop.
    start = time.time()

    mode, in_dirs, gd_params, gc_params, cr_params, er_params,\
    gr_params, pv_params, da_params, ps_params, bf_params, \
    sc_params, ga_params, rm_params, pl_params, flag_move_file, axes_params =\
    gip_params
    input_dir, output_dir, done_dir = in_dirs

    # Generate output subdir.
    output_subdir = join(output_dir, sub_dir)
    # Check if subdir already exists, if not create it.
    if not exists(output_subdir):
        mkdir(output_subdir)

    # Store name of file in 'myfile'.
    myfile = dir_files[1][f_indx]
    # Store cluster's name
    clust_name = myfile[:-4]
    print 'Analizing cluster %s.' % (clust_name)

    # Get data from semi-data input file.
    mode, semi_return = g_s(input_dir, clust_name, mode)

    # Get cluster's photometric data from file.
    phot_data = gd(input_dir, sub_dir, myfile, gd_params)
    # If Manual mode is set, display frame and ask if it should be trimmed.
    phot_data = t_f(phot_data, mode)
    # Unpack coordinates, magnitude and color.
    x_data, y_data, mag_data, col1_data = phot_data[1], phot_data[2], \
    phot_data[3], phot_data[5]
    print 'Data obtained from input file (N stars: %d).' % len(phot_data[0])

    # Get cluster's center values and errors, filtered 2D hist, non-filtered
    # 2D hist, x,y bin centers and width of each bin
    # used
    center_params = g_c(x_data, y_data, gc_params, mode, semi_return)
    # Unpack values from list.
    bin_width, h_not_filt, hist_xyedges, bin_center = center_params[0][0], \
    center_params[1], center_params[2], center_params[4]
    center_cl = [center_params[5][0][0], center_params[5][0][1]]

    # Get density profile
    rdp_params = gdp(h_not_filt, bin_center, bin_width)
    radii, ring_density = rdp_params[:2]
    print 'Density profile calculated.'

    # Get background value in stars/px^2.
    backg_value = gbg(ring_density)
    print 'Background calculated (%0.5f stars/px^2).' % backg_value

    # Get cluster radius
    radius_params = gcr(phot_data, backg_value, cr_params, center_params,
        rdp_params, semi_return, mode)
    clust_rad = radius_params[0]

    # Get approximate number of cluster's members.
    n_c, flag_num_memb_low, a_clust, n_clust = g_m_n(backg_value, clust_rad,
        rdp_params, bin_width)
    print 'Approximate number of members in cluster obtained (%d).' % (n_c)

    # Get contamination index.
    cont_index = g_c_i(backg_value, a_clust, n_clust)
    print 'Contamination index obtained (%0.2f).' % cont_index

    # Get King profiles based on the density profiles.
    delta_xy = max((max(x_data) - min(x_data)), (max(y_data) - min(y_data)))
    k_prof, k_pr_err, d_b_k, n_c_k, flag_king_no_conver = \
    gkp(clust_rad, backg_value, radii, ring_density, delta_xy, x_data, y_data,
        bin_width)

    # Accept and reject stars based on their errors.
    popt_mag, popt_col1, acpt_stars, rjct_stars, err_plot, rjct_errors_fit = \
    ear(phot_data, axes_params, er_params, mode, semi_return)

    # Get stars in and out of cluster's radius.
    stars_in, stars_out, stars_in_rjct, stars_out_rjct = gio(center_cl,
        clust_rad, acpt_stars, rjct_stars)
    print "Stars separated in/out of cluster's boundaries."

    # Obtain manual 2D histogram for the field with star's values attached
    # to each bin.
    H_manual = mh(phot_data, hist_xyedges)
    print 'Manual 2D histogram obtained.'

    # Get cluster + field regions around the cluster's center.
    flag_area_stronger, cluster_region, field_region = g_r(bin_center,
        bin_width, h_not_filt, clust_rad, H_manual, stars_in, stars_out,
        gr_params)
    print 'Cluster + field stars regions obtained (%d).' % len(field_region)

    # Calculate integrated magnitude.
    integr_return = g_i_m(center_cl, clust_rad, cluster_region, field_region,
        flag_area_stronger)
    print 'Integrated color magnitude distribution obtained (%0.2f).' % \
        (integr_return[5] - integr_return[2])

    # Check if test is to be applied or skipped.
    flag_pval_test = pv_params[0]
    # Get physical cluster probability based on p_values distribution.
    if flag_pval_test:
        # Check if field regions where found.
        if not flag_area_stronger:
            # pval_test_params = prob_cl_kde, p_vals_cl, p_vals_f, kde_cl_1d,
            #                    kde_f_1d, x_kde, y_over
            pval_test_params = g_pv(cluster_region, field_region, col1_data,
                                    mag_data, center_cl, clust_rad, pv_params)
            # Add flag to list.
            pval_test_params = pval_test_params + [flag_pval_test]
            print 'Probability of physical cluster obtained (%0.2f).' % \
            pval_test_params[0]

            # Get QQ plot for p-values distributions.
            # qq_params = ccc, quantiles, r_squared, slope, intercept
            qq_params = g_qq(pval_test_params[1], pval_test_params[2])
            print 'QQ-plot obtained (CCC = %0.2f).' % qq_params[0]

    # Skip process.
    if not flag_pval_test or flag_area_stronger:
        print 'Skipping p-value test for cluster.'
        # Pass empty lists to make_plots.
        pval_test_params, qq_params = [-1., False], [-1.]

    # Apply decontamination algorithm if at least one equal-sized field region
    # was found around the cluster.
    print 'Applying decontamination algorithm.'
    decont_algor_return = fdb(flag_area_stronger, cluster_region, field_region,
                            col1_data, mag_data, center_cl, clust_rad,
                            clust_name, sub_dir, da_params)
    memb_prob_avrg_sort = decont_algor_return[0]

    # Create data file with membership probabilities.
    c_m_f(output_dir, sub_dir, clust_name, memb_prob_avrg_sort)
    print 'Membership probabilities saved to file.'

    # Get the completeness level for each magnitude bin. This will be used by
    # the isochrone/synthetic cluster fitting algorithm.
    completeness = m_c(mag_data)
    print 'Completeness magnitude levels obtained.'

    # Store all isochrones in all the metallicity files in isoch_list.
    # Store metallicity values and isochrones ages between the allowed
    # ranges in isoch_ma; extinction and distance modulus values in isoch_ed.
    # isoch_list, isoch_ma, isoch_ed = ip_list
    # Only read files if best fit process is set to run.
    # bf_flag = bf_params[0]
    if bf_params[0]:
        ip_list = ip(ps_params)
        print 'Theoretical isochrones read and stored (%d).' % \
            (len(ip_list[0]) * len(ip_list[0][0]))
    else:
        ip_list = []

    # Reduce number of stars in cluster according to a lower membership
    # probability limit.
    red_return = rm(decont_algor_return, bf_params, rm_params)
    red_memb_prob = red_return[0]

    # Obtain best fitting parameters for cluster.
    err_lst = [popt_mag, popt_col1, er_params[2]]
    bf_return = bfsc(err_lst, red_memb_prob, completeness, ip_list, bf_params,
        sc_params, ga_params, ps_params)

    # New name for cluster? Useful when there's a single photometric file
    # with multiple clusters in it.
    if mode == 'm':
        wrong_answer = True
        while wrong_answer:
            answer_rad = raw_input('New name for cluster? (y/n) ')
            if answer_rad == 'n':
                wrong_answer = False
            elif answer_rad == 'y':
                new_name = str(raw_input('Input new name: '))
                clust_name = new_name
                wrong_answer = False
            else:
                print 'Wrong input. Try again.\n'

    # Add cluster data and flags to output file
    a_d_o(out_file_name, sub_dir, output_dir, clust_name, center_params,
        radius_params, k_prof, k_pr_err, n_c_k, flag_king_no_conver, cont_index,
        n_c, pval_test_params[0], qq_params[0], integr_return,
        rjct_errors_fit, flag_num_memb_low, bf_return)
    print 'Data added to output file.'

    # Make plots
    if pl_params[0]:
        mp(output_subdir, clust_name, x_data, y_data, center_params, rdp_params,
            backg_value, radius_params[0:3],
            cont_index, mag_data, col1_data, popt_mag, popt_col1,
            err_plot, rjct_errors_fit, k_prof, k_pr_err, d_b_k,
            flag_king_no_conver, stars_in, stars_out, stars_in_rjct,
            stars_out_rjct, integr_return, n_c, flag_area_stronger,
            cluster_region, field_region, pval_test_params, qq_params,
            memb_prob_avrg_sort, completeness, bf_params, red_return,
            bf_return, ga_params, er_params, axes_params, ps_params, pl_params)
        print 'Plots created.'

    # Move file to 'done' dir.
    if flag_move_file:
        dst_dir = join(done_dir, sub_dir)
        # If the sub-dir doesn't exist, create it before moving the file.
        if not exists(dst_dir):
            makedirs(dst_dir)
        shutil.move(join(input_dir, sub_dir, myfile), dst_dir)
        # Also move *memb_data.dat file if it exists.
        if isfile(join(input_dir, sub_dir, clust_name + '_memb.dat')):
            shutil.move(join(input_dir, sub_dir, clust_name + '_memb.dat'),
                dst_dir)
        # If sub-dir left behind is empty, remove it.
        try:
            rmdir(join(input_dir, sub_dir))
        except OSError:
            # Sub-dir not empty, skip.
            pass
        print 'Photometric data file moved.'

    elapsed = time.time() - start
    m, s = divmod(elapsed, 60)
    print 'End of analysis for %s in %dm %02ds.\n' % (clust_name, m, s)

    # Force the Garbage Collector to release unreferenced memory.
    gc.collect()


# Begin code.
print '            OCAAT v1.0\n'
print '-------------------------------------------\n'

# Path where the code is running
mypath = realpath(join(getcwd(), dirname(__file__)))

# Read input parameters from ocaat_input.dat file.
gip_params = gip(mypath)
# Read input/output paths.
input_dir, output_dir = gip_params[1][:2]
# Read paths and names of all clusters stored in input_dir.
dir_files = rp(input_dir)
# Create output data file in output_dir (append if file already exists)
out_file_name = c_o_d_f(output_dir)

# Iterate through all cluster files.
for f_indx, sub_dir in enumerate(dir_files[0]):

    try:
        # Call main function for this cluster.
        ocaat_main(f_indx, sub_dir, out_file_name, gip_params)
    except Exception, err:
        import traceback
        print 'FATAL: %s could not be processed. ' % dir_files[1][f_indx][:-4]
        print 'Error:', str(err), '\n'
        print traceback.format_exc()

print 'Full iteration completed.'
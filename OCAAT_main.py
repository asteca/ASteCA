# -*- coding: utf-8 -*-


from os.path import join, realpath, dirname, exists
import matplotlib.pyplot as plt
import time

# Import files with defined functions
import functions.get_in_params as gip

from functions.create_out_data_file import create_out_data_file as c_o_d_f
from functions.get_data_semi import get_semi as g_s
import functions.get_phot_data as gd
from functions.display_frame import disp_frame as d_f
from functions.trim_frame import trim_frame as t_f
import functions.get_center as g_c
from functions.display_cent import disp_cent as d_c
from functions.manual_histo import manual_histo as mh
import functions.get_background as gbg
import functions.get_dens_prof as gdp
import functions.get_radius as gr
import functions.get_king_prof as gkp
from functions.display_rad import disp_rad as d_r
import functions.err_accpt_rejct as ear
from functions.display_errors import disp_errors as d_e
from functions.err_accpt_rejct_max import err_a_r_m as e_a_r_m
import functions.get_in_out as gio
from functions.get_integ_mag import integ_mag as g_i_m
import functions.get_members_number as g_m_n
import functions.get_cont_index as g_c_i
from functions.get_regions import get_regions as g_r

from functions.field_decont_kde import field_decont_kde as fdk
from functions.get_p_value import get_pval as g_pv
from functions.get_qqplot import qqplot as g_qq
from functions.memb_prob_avrg_sort import mpas as m_p_a_s
from functions.get_completeness import mag_completeness as m_c
from functions.get_isoch_params import ip
from functions.best_fit_synth_cl import best_fit as bfsc

from functions.make_plots import make_plots as mp
from functions.add_data_output import add_data_output as a_d_o
from functions.cl_members_file import cluster_members_file as c_m_f

from os import makedirs
from os import listdir, getcwd, walk, mkdir, rmdir
import shutil
# Garbage collector.
import gc


print '            OCAAT v1.0\n'
print '-------------------------------------------\n'

# Path where the code is running
mypath = realpath(join(getcwd(), dirname(__file__)))

# Read input parameters for code from file.
mode, in_dirs, gd_params, gc_params, br_params, cr_params, er_params,\
gr_params, pv_params, da_params, ps_params, bf_params, sc_params, ga_params,\
flag_make_plot, flag_move_file, axes_params = gip.get_in_params(mypath)

# Read paths.
mypath2, mypath3, output_dir = in_dirs

# Create output data file (append if file already existis)
c_o_d_f(output_dir)

# Store subdir names [0] and file names [1] inside each subdir.
dir_files = [[], []]
for root, dirs, files in walk(mypath2):
    if dirs:
        for subdir in dirs:
            for name in listdir(join(mypath2, subdir)):
                # Check to see if it's a valid data file.
                if name.endswith(('.DAT', '.MAG', '.OUT', '.TEX')):
                    dir_files[0].append(subdir)
                    dir_files[1].append(name)


# Iterate through all cluster files.
for f_indx, sub_dir in enumerate(dir_files[0]):

    # Store name of file in 'myfile'.
    myfile = dir_files[1][f_indx]

    # Start timing this loop.
    start = time.time()

    # Store cluster's name
    clust_name = myfile[:-4]
    print 'Analizing cluster %s.' % (clust_name)

    # If Semi mode set, get data from data iput file.
    if mode == 's':
        cent_cl_semi, cl_rad_semi, cent_flag_semi, rad_flag_semi, \
        err_flag_semi = g_s(mypath, clust_name)

    # Get cluster's photometric data from file.
    phot_data = gd.get_data(mypath2, sub_dir, myfile, gd_params)
    x_data, y_data, mag_data, col1_data = phot_data[1], phot_data[2], \
    phot_data[3], phot_data[5]
    print 'Data correctly obtained from input file (N stars: %d).'\
    % len(phot_data[0])

    # If Manual mode is set, display frame and ask if it should be trimmed.
    if mode == 'm':
        # Show plot with center obtained.
        d_f(x_data, y_data, mag_data)
        plt.show()

        wrong_answer = True
        while wrong_answer:
            temp_cent, temp_side = [], []
            answer_fra = raw_input('Trim frame? (y/n) ')
            if answer_fra == 'n':
                wrong_answer = False
            elif answer_fra == 'y':
                print 'Input center of new frame (in px).'
                temp_cent.append(float(raw_input('x: ')))
                temp_cent.append(float(raw_input('y: ')))
                print 'Input side lenght for new frame (in px).'
                temp_side.append(float(raw_input('x_side: ')))
                temp_side.append(float(raw_input('y_side: ')))
                # Trim frame.
                phot_data = t_f(temp_cent, temp_side, phot_data)
                wrong_answer = False
            else:
                print 'Wrong input. Try again.\n'
    elif mode == 's':
        # If there are too many stars in the frame, trim it.
        if len(phot_data[0]) > 25000:
            temp_cent, temp_side = [], []
            # Set center.
            temp_cent.append(cent_cl_semi[0])
            temp_cent.append(cent_cl_semi[1])
            # Set side length.
            temp_side.append(2000.)
            temp_side.append(2000.)
            # Trim frame.
            phot_data = t_f(temp_cent, temp_side, phot_data)


    # Get cluster's center values and errors, set of 4 filtered 2D hist,
    # set of 4 non-filtered 2D hist, x,y bin centers and width of each bin
    # used
    center_cl, cent_cl_err, h_filter, h_not_filt, xedges_min_db, \
    yedges_min_db, x_center_bin, y_center_bin, width_bins, flag_center, \
    flag_std_dev = g_c.get_center(x_data, y_data, gc_params)
    print 'Auto center found: (%d, %d) px.' % (center_cl[0], center_cl[1])


    # If Manual mode is set, display center and ask the user to accept it or
    # input new one.
    flag_center_manual = False
    if mode == 'm':

        # Show plot with center obtained.
        d_c(x_data, y_data, mag_data, center_cl, cent_cl_err, x_center_bin, \
        y_center_bin, h_filter)
        plt.show()

        wrong_answer = True
        while wrong_answer:
            answer_cen = raw_input('Input new center values? (y/n) ')
            if answer_cen == 'n':
                print 'Value accepted.'
                wrong_answer = False
            elif answer_cen == 'y':
                print 'Input new center values (in px).'
                center_cl[0] = float(raw_input('x: '))
                center_cl[1] = float(raw_input('y: '))
                # Update values.
                cent_cl_err[0], cent_cl_err[1] = 13., 13.
                # Store center values in bin coordinates. We substract
                # the min (x,y) coordinate values otherwise the bin
                # coordinates won't be aligned.
                x_center_bin[0], y_center_bin[0] = int(round((center_cl[0]-\
                min(x_data))/width_bins[0])), \
                int(round((center_cl[1]-min(y_data))/width_bins[0]))
                wrong_answer = False
                flag_center_manual = True
            else:
                print 'Wrong input. Try again.\n'
    elif mode == 's':
        if cent_flag_semi == 1:
            # Update center values.
            center_cl[0] = cent_cl_semi[0]
            center_cl[1] = cent_cl_semi[1]
            # Update error values.
            cent_cl_err[0], cent_cl_err[1] = 13., 13.
            print 'Semi-auto center set: (%d, %d) px.' % (center_cl[0], center_cl[1])
            # Store center values in bin coordinates. We substract
            # the min (x,y) coordinate values otherwise the bin
            # coordinates won't be aligned.
            x_center_bin[0], y_center_bin[0] = int(round((center_cl[0]-\
            min(x_data))/width_bins[0])), \
            int(round((center_cl[1]-min(y_data))/width_bins[0]))
            flag_center_manual = True


    # Obtain manual 2D histogram for the field with star's values attached
    # to each bin.
    H_manual = mh(phot_data, xedges_min_db, yedges_min_db)
    print 'Manual 2D histogram obtained.'


    # Get background value in stars/area
    # Inner and outer radii for obtaining the background values. Both are
    # calculated as a given fraction of the minimum width between the x and y
    # axis spans.
    inn_fr, out_fr = br_params
    inner_ring = (min((max(x_data)-min(x_data)),
                      (max(y_data)-min(y_data))))/inn_fr
    outer_ring = (min((max(x_data)-min(x_data)),
                      (max(y_data)-min(y_data))))/out_fr

    if mode == 'm':
        wrong_answer = True
        while wrong_answer:
            answer_bkg = raw_input('Input new inner and outer radius for \
background (%d, %d) px? (y/n) ' % (inner_ring, outer_ring))
            if answer_bkg == 'n':
                print 'Values accepted.'
                wrong_answer = False
            elif answer_bkg == 'y':
                print 'Input new values (in px).'
                inner_ring = float(raw_input('inner_ring: '))
                outer_ring = float(raw_input('outer_ring: '))
                wrong_answer = False
            else:
                print 'Wrong input. Try again.\n'

    backg_value, flag_bin_count = gbg.get_background(x_data, y_data, \
                                    x_center_bin, y_center_bin, h_not_filt,\
                                    width_bins, inner_ring, outer_ring)
    print 'Background calculated (%0.5f stars/px^2).' % backg_value


    # Get density profile
    radii, ring_density, poisson_error = gdp.get_dens_prof(h_not_filt, \
    x_center_bin, y_center_bin, width_bins, inner_ring)
    print 'Density profile calculated.'


    # Get cluster radius
    radius_params = gr.get_clust_rad(backg_value, radii, ring_density,
                                     width_bins, cr_params)
    clust_rad = radius_params[0]
    print 'Radius calculated: %d px.' % clust_rad


    # If Manual mode is set, display radius and ask the user to accept it or
    # input new one.
    flag_radius_manual = False
    if mode == 'm':
        print 'Radius found: ', clust_rad
        d_r(x_data, y_data, mag_data, center_cl, cent_cl_err,
            radius_params[0:3], x_center_bin, y_center_bin, h_filter,\
            backg_value, radii, ring_density, clust_name, poisson_error,\
            width_bins, inner_ring, outer_ring)
        plt.show()

        wrong_answer = True
        while wrong_answer:
            answer_rad = raw_input('Accept radius (otherwise input new one \
manually)? (y/n) ')
            if answer_rad == 'y':
                print 'Value accepted.'
                wrong_answer = False
            elif answer_rad == 'n':
                clust_rad = float(raw_input('Input new radius value (in \
px): '))
                wrong_answer = False
                flag_radius_manual = True
            else:
                print 'Wrong input. Try again.\n'
    elif mode == 's':
        if rad_flag_semi == 1:
            # Update value.
            radius_params[0] = cl_rad_semi
            clust_rad = radius_params[0]
            flag_radius_manual = True


    # Get King profiles based on the density profiles.
    k_prof, k_pr_err, n_c_k, flag_king_no_conver = \
    gkp.get_king_profile(clust_rad, backg_value, radii, ring_density)
    if flag_king_no_conver == False:
        print '3-P King profile obtained.'
    else:
        print 'King profile fitting did not converge.'


    # Apply auto rejecting of errors if flag is True.
    e_max = er_params[2]
    # Accept and reject stars based on their errors.
    popt_mag, popt_col1, acpt_stars, rjct_stars, err_plot = \
    ear.err_accpt_rejct(phot_data, er_params)
    # This indicates if we are to use the output of the 'err_accpt_rejct'
    # function or all stars with errors < e_max.
    rjct_errors_fit = False

    # If list of accepted stars is empty, halt the code.
    if not acpt_stars:
        print '  No stars accepted based on their errors.'
        print '  This probably means the exponential error function\n\
did not converge.'
        print 'Using all stars with errors < %0.2f.' % e_max
        # Call function to reject stars with errors > e_max.
        popt_mag, popt_col1, acpt_stars, rjct_stars = e_a_r_m(phot_data,
                                                              er_params)
        rjct_errors_fit = True
    else:
        print 'Stars accepted/rejected based on their errors.'

    # Is 'semi' is set, check for the flag that indicates whether to use
    # auto errors rejecting or all stars with errors < e_max.
    if mode == 's':
        if err_flag_semi == 1:
            # Reject error fit.
            print 'Semi: using all stars with errors < %0.2f.' % e_max
            # Call function to reject stars w errors > e_max.
            popt_mag, popt_col1, acpt_stars, rjct_stars = e_a_r_m(phot_data,
                                                                  er_params)
            rjct_errors_fit = True

    # If 'manual' mode is set, display errors distributions and ask the user
    # to accept it or else use all stars except those with errors > e_max in
    # either the magnitude or the color.
    elif mode == 'm':
        print 'Plot error distributions.'
        # Display automatic errors rejection.
        d_e(mag_data, popt_mag,  popt_col1, acpt_stars, rjct_stars, err_plot,
            er_params, axes_params)
        plt.show()
        # Ask if keep or reject.
        wrong_answer = True
        while wrong_answer:
            answer_rad = raw_input('Accept fit for errors (otherwise use \
all stars with photom errors < %0.2f)? (y/n) ' % e_max)
            if answer_rad == 'y':
                print 'Fit accepted.'
                wrong_answer = False
            elif answer_rad == 'n':
                print 'Using stars with errors < %0.2f.' % e_max
                # Call function to reject stars w errors > e_max.
                popt_mag, popt_col1, acpt_stars, rjct_stars = e_a_r_m(phot_data,
                                                                      er_params)
                rjct_errors_fit = True
                wrong_answer = False
            else:
                print 'Wrong input. Try again.\n'


    # Get stars in and out of cluster's radius.
    stars_in, stars_out, stars_in_rjct, stars_out_rjct =  \
    gio.get_in_out(center_cl, clust_rad, acpt_stars, rjct_stars)
    print "Stars separated in/out of cluster's boundaries."


    # Calculate integrated magnitude.
    stars_in_mag, stars_in_all_mag = g_i_m(stars_in, stars_in_rjct)
    print 'Integrated magnitude distribution obtained.'


    # Get approximate number of cluster's members.
    n_c, flag_num_memb_low = g_m_n.get_memb_num(backg_value, clust_rad,
                                                stars_in, stars_in_rjct)
    print 'Approximate number of members in cluster obtained (%d).' % (n_c)


    # Get contamination index.
    cont_index = g_c_i.cont_indx(backg_value, clust_rad, stars_in,
                                 stars_in_rjct)
    print 'Contamination index obtained (%0.2f).' % cont_index


    # Get cluster + field regions around the cluster's center.
    flag_area_stronger, cluster_region, field_region = \
    g_r(x_center_bin, y_center_bin, width_bins, h_not_filt, clust_rad,
        H_manual, stars_in, stars_out, gr_params)
    print 'Cluster + field stars regions obtained (%d).' % len(field_region)


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
            print 'QQ-plot obtained (CCC = %0.2f)' % qq_params[0]

    # Skip process.
    if not flag_pval_test or flag_area_stronger:
        print 'Skipping p-value test for cluster.'
        # Pass empty lists to make_plots.
        pval_test_params, qq_params = [-1., False], [-1.]



    # Apply decontamination algorithm if at least one equal-sized field region
    # was found around the cluster.
    print 'Applying decontamination algorithm.'
    runs_fields_probs = fdk(flag_area_stronger, cluster_region, field_region,
                            col1_data, mag_data, center_cl, clust_rad,
                            clust_name, sub_dir, da_params)


    # Average and sort all membership probabilities for each star.
    memb_prob_avrg_sort, clust_reg_prob_avrg = m_p_a_s(cluster_region, \
    runs_fields_probs, center_cl, clust_rad)
    print 'Averaged probabilities for all runs.'


    # Get the completeness level for each magnitude bin. This will be used by
    # the isochrone/synthetic cluster fitting algorithm.
    completeness = m_c(mag_data)
    print 'Completeness magnitude levels obtained.'


    # Store all isochrones in all the metallicity files in isoch_list.
    # Store metallicity values and isochrones ages between the allowed
    # ranges in isoch_ma; extinction and distance modulus values in isoch_ed.
    # isoch_list, isoch_ma, isoch_ed = ip_list
    ip_list = ip(ps_params)
    print 'Theoretical isochrones read and stored.'


    # Obtain best fitting parameters for cluster.
    print 'Searching for optimal parameters.'
    err_lst = [popt_mag, popt_col1, e_max]
    bf_return = bfsc(err_lst, memb_prob_avrg_sort, completeness, ip_list,
                     bf_params, sc_params, ga_params, ps_params)
    print 'Best fit parameters obtained.'


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


    # Create data file with membership probabilities.
    c_m_f(output_dir, sub_dir, clust_name, memb_prob_avrg_sort)
    print 'Membership probabilities for stars in cluster region saved to file.'


    # Generate output subdir.
    output_subdir = join(output_dir, sub_dir)
    # Check if subdir already exists, if not create it
    if not exists(output_subdir):
        makedirs(output_subdir)

    # Add cluster data and flags to output file
    a_d_o(sub_dir, output_dir, clust_name, center_cl, clust_rad, k_prof,
          n_c_k, flag_king_no_conver, cont_index, n_c, pval_test_params[0],
          qq_params[0], stars_in_mag, flag_center, flag_std_dev,
          flag_center_manual, flag_radius_manual, rjct_errors_fit,
          flag_bin_count, radius_params[3:], flag_num_memb_low, bf_return)
    print 'Data added to output file.'

    # Make plots
    if flag_make_plot:
        mp(output_subdir, clust_name, x_data, y_data, center_cl, cent_cl_err,
           x_center_bin, y_center_bin, h_filter, radii, backg_value, inner_ring,
           outer_ring, radius_params[0:3], ring_density, poisson_error,
           cont_index, width_bins, mag_data, col1_data, popt_mag, popt_col1,
           err_plot, rjct_errors_fit, k_prof, k_pr_err, flag_king_no_conver,
           stars_in, stars_out, stars_in_rjct, stars_out_rjct, stars_in_mag,
           stars_in_all_mag, n_c, flag_area_stronger,
           cluster_region, field_region, pval_test_params, qq_params,
           clust_reg_prob_avrg, memb_prob_avrg_sort, bf_params, bf_return,
           ga_params, er_params, axes_params, ps_params)
        print 'Plots created.'
    
    # Move file to 'done' dir.
    if flag_move_file:
        dst_dir = join(mypath3, sub_dir)
        # If the sub-dir doesn't exist, create it before moving the file.
        if not exists(dst_dir):
            mkdir(dst_dir)
        shutil.move(join(mypath2, sub_dir, myfile), dst_dir)
        # If sub-dir left behind is empty, remove it.
        try:
            rmdir(join(mypath2, sub_dir))
        except OSError as ex:
            # Sub-dir not empty, skip.
            pass
        print 'Data file moved to folder: %s' % dst_dir


    elapsed = time.time() - start
    m, s = divmod(elapsed, 60)
    print 'End of analysis for %s in %dm %02ds.\n'% (clust_name, m, s)

    # Force the Garbage Collector to release unreferenced memory.
    gc.collect()

    # Store memory information.
#    h = hpy()
#    with open('/home/gabriel/clusters/mem_out', "a") as f:
#        f.write('\n')
#        f.write(strftime("%Y-%m-%d %H:%M:%S")+'\n')
#        f.write(str(h.heap()))
#        f.write('\n')

print '\nFull iteration completed.'
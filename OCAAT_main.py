# -*- coding: utf-8 -*-


# Import these to be able to grab the dir from where the code is running and
# iterate for all data files in that dir
from os import listdir, getcwd, walk, mkdir, rmdir
from os.path import join, realpath, dirname, exists
from os import makedirs
from os.path import expanduser

import matplotlib.pyplot as plt

import shutil
import time
import gc

# Import files with defined functions
from functions.create_out_data_file import create_out_data_file as c_o_d_f
from functions.get_isochrones import get_isochrones as g_i
from functions.get_data_semi import get_semi as g_s
import functions.get_data as gd
from functions.display_frame import disp_frame as d_f
from functions.trim_frame import trim_frame as t_f
import functions.get_center as g_c
from functions.display_cent import disp_cent as d_c
from functions.manual_histo import manual_histo as mh
import functions.get_background as gbg
import functions.get_dens_prof as gdp
import functions.get_clust_rad as gcr
import functions.get_king_prof as gkp
from functions.display_rad import disp_rad as d_r
import functions.err_accpt_rejct as ear
from functions.display_errors import disp_errors as d_e
from functions.err_accpt_rejct_03 import err_a_r_03 as e_a_r_03
import functions.get_in_out as gio
from functions.get_integ_mag import integ_mag as g_i_m
import functions.get_members_number as g_m_n
import functions.get_cont_index as g_c_i
from functions.get_regions import get_regions as g_r

#from functions.field_decont_ran import field_decont_ran as fdr
#from functions.field_decont_VB import field_decont_VB as fdvb
from functions.field_decont_kde import field_decont_kde as fdk
#from functions.field_decont_dias import field_decont_dias as fdd

from functions.get_p_value import get_pval as g_pv
from functions.get_qqplot import qqplot as g_qq
from functions.memb_prob_avrg_sort import mpas as m_p_a_s
from functions.make_plots import make_plots as mp
from functions.add_data_output import add_data_output as a_d_o
from functions.cl_members_file import cluster_members_file as c_m_f
        


print '            OCAAT\n'
print '-------------------------------------------\n'


# Select automatic or manual mode. Manual mode will ask for some values (center,
# radius, background) but will give you the chance to have them automatically
# calculated. Automatic mode does not require user intervention of any kind.
wrong_mode = True
while wrong_mode:
    mode = raw_input('\nSelect Auto / Semi / Manual mode (a/s/m): ')
    if mode == 'a':
        print 'Automatic mode selected.\n'
        wrong_mode = False
    elif mode == 's':
        print 'Semi automatic mode selected.\n'
        wrong_mode = False
    elif mode == 'm':
        print 'Manual mode selected.\n'
        wrong_mode = False
    else:
        print 'Wrong input. Try again.\n'



# Define paths to folders/files.

# Set 'home' dir.
home = expanduser("~")

# Path where the code is running
mypath = realpath(join(getcwd(), dirname(__file__)))
# Path where the cluster files to be processed are stored.
mypath2 = home+'/clusters/clusters_in'
# Path where the sub-folders are moved once processed.
mypath3 = home+'/clusters/clusters_done'
# Set output directory.        
output_dir = home+'/clusters/clusters_out/'


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


    # If Semi mode set, get data from 'data_input' file.
    if mode == 's':
        cent_cl_semi, cl_rad_semi, cent_flag_semi, rad_flag_semi, \
        err_flag_semi = g_s(mypath, clust_name) 


    # Get cluster's photometric data from file.
    id_star, x_data, y_data, mag_data, e_mag, col1_data, e_col1 = \
    gd.get_data(mypath2, sub_dir, myfile)
    print 'Data correctly obtained from input file  (N stars: %d).' % len(x_data)


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
                id_star, x_data, y_data, mag_data, e_mag, col1_data, e_col1 =\
                t_f(temp_cent, temp_side, id_star, x_data, y_data, mag_data, 
                    e_mag, col1_data, e_col1)
                wrong_answer = False
            else:
                print 'Wrong input. Try again.\n'
    elif mode == 's':
        # If there are too many stars in the frame, trim it.
        if len(x_data) > 25000:
            temp_cent, temp_side = [], []
            # Set center.
            temp_cent.append(cent_cl_semi[0])
            temp_cent.append(cent_cl_semi[1])
            # Set side length.
            temp_side.append(2000.)
            temp_side.append(2000.)
            # Trim frame.
            id_star, x_data, y_data, mag_data, e_mag, col1_data, e_col1 =\
            t_f(temp_cent, temp_side, id_star, x_data, y_data, mag_data, 
                e_mag, col1_data, e_col1)


    # Get cluster's center values and errors, set of 4 filtered 2D hist,
    # set of 4 non-filtered 2D hist, x,y bin centers and width of each bin
    # used
    center_cl, cent_cl_err, h_filter, h_not_filt, xedges_min_db, \
    yedges_min_db, x_center_bin, y_center_bin, width_bins, flag_center, \
    flag_std_dev = g_c.get_center(x_data, y_data)
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
    H_manual = mh(id_star, x_data, y_data, mag_data, e_mag, col1_data, e_col1, \
    xedges_min_db, yedges_min_db)
    print 'Manual 2D histogram obtained.'
    
    
    # Get background value in stars/area
    
    # Inner and outer radii for obtaining the background values. Both are 
    # calculated as a given fraction of the minimum widht between the x and y 
    # axis spans.
    inner_ring = (min((max(x_data)-min(x_data)), (max(y_data)-min(y_data))))/4.
    outer_ring = (min((max(x_data)-min(x_data)), (max(y_data)-min(y_data))))/3.  
    
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
    print 'Background calculated (%0.5f stars/px^2).' % backg_value[0]


    # Get density profile
    radii, ring_density, poisson_error = gdp.get_dens_prof(h_not_filt, \
    x_center_bin, y_center_bin, width_bins, inner_ring)
    print 'Density profile calculated.'
    
    
    # Get cluster radius
    clust_rad, delta_backg, delta_percentage, \
    flag_delta_total, flag_not_stable, flag_rad_500, flag_delta, \
    flag_delta_points = gcr.get_clust_rad(backg_value, radii, 
                                          ring_density, width_bins)
    print 'Radius calculated: %d px.' % clust_rad[0]


    # If Manual mode is set, display radius and ask the user to accept it or
    # input new one.
    flag_radius_manual = False
    if mode == 'm':
        print 'Radius found: ', clust_rad[0]
        d_r(x_data, y_data, mag_data, center_cl, cent_cl_err, clust_rad, \
        x_center_bin, y_center_bin, h_filter, backg_value, radii,
        delta_backg, ring_density, clust_name, poisson_error, width_bins,
        delta_percentage, inner_ring, outer_ring)
        plt.show()
                
        wrong_answer = True
        while wrong_answer:
            answer_rad = raw_input('Accept radius (otherwise input new one \
manually)? (y/n) ')
            if answer_rad == 'y':
                print 'Value accepted.'
                wrong_answer = False
            elif answer_rad == 'n':
                clust_rad[0] = float(raw_input('Input new radius value (in \
px): '))
                wrong_answer = False
                flag_radius_manual = True
            else:
                print 'Wrong input. Try again.\n'
    elif mode == 's':
        if rad_flag_semi == 1:
            # Update value.
            clust_rad[0] = cl_rad_semi
            flag_radius_manual = True

                
    # Accept and reject stars based on their errors.
    bright_end, popt_mag, popt_umag, pol_mag, popt_col1, popt_ucol1, \
    pol_col1, mag_val_left, mag_val_right, col1_val_left, col1_val_right, \
    acpt_stars, rjct_stars = ear.err_accpt_rejct(id_star, x_data, y_data,
                                                 mag_data, e_mag, col1_data,
                                                 e_col1)
    print 'Stars accepted/rejected based on their errors.'


    # Get King profiles based on the density profiles.
    k_prof, k_pr_err, n_c_k, flag_king_no_conver = \
    gkp.get_king_profile(clust_rad, backg_value, radii, ring_density)
    if flag_king_no_conver == False:
        print '3-P King profile obtained.'
    else:
        print 'King profile fitting did not converge.'

    
    # This indicates if we are to use the output of the 'err_accpt_rejct'
    # function or all stars with errors < 0.3.
    use_errors_fit = True
    # If Manual mode is set, display errors distributions and ask the user
    # to accept it or else use all stars except those with errors >0.3 in
    # either the magnitude or the color.
    flag_errors_manual = False
    if mode == 'm':
        print 'Plot error distributions.'
        d_e(mag_data, bright_end, popt_mag, popt_umag, pol_mag, popt_col1,
            popt_ucol1, pol_col1, mag_val_left, mag_val_right,
            col1_val_left, col1_val_right, acpt_stars, rjct_stars)
        plt.show()
                
        wrong_answer = True
        while wrong_answer:
            answer_rad = raw_input('Accept fit for errors (otherwise use \
all stars with photom errors < 0.3)? (y/n) ')
            if answer_rad == 'y':
                print 'Fit accepted.'
                wrong_answer = False
            elif answer_rad == 'n':
                print 'Using stars with errors < 0.3.'
                # Call function to reject stars w errors > 0.3.
                acpt_stars, rjct_stars = e_a_r_03(id_star, x_data, y_data,
                                                  mag_data, e_mag, col1_data,
                                                  e_col1)
                flag_errors_manual = True
                use_errors_fit = False
                wrong_answer = False
            else:
                print 'Wrong input. Try again.\n'
    elif mode == 's':
        if err_flag_semi == 1:
            # Reject error fit.
            print 'Using stars with errors < 0.3.'
            # Call function to reject stars w errors > 0.3.
            acpt_stars, rjct_stars = e_a_r_03(id_star, x_data, y_data,
                                              mag_data, e_mag, col1_data,
                                              e_col1)
            flag_errors_manual = True
            use_errors_fit = False


    # Get stars in and out of cluster's radius.
    stars_in, stars_out, stars_in_rjct, stars_out_rjct =  \
    gio.get_in_out(center_cl, clust_rad[0], acpt_stars, rjct_stars)
    print "Stars separated in/out of cluster's boundaries."
    
    
    # Calculate integrated magnitude.
    stars_in_mag, stars_in_all_mag = g_i_m(stars_in, stars_in_rjct)
    print 'Integrated magnitude distribution obtained.'
    
    
    # Get approximate number of cluster's members.
    n_c, flag_num_memb_low, flag_no_memb = g_m_n.get_memb_num(backg_value[0],\
    clust_rad[0], stars_in, stars_in_rjct)
    print 'Approximate number of members in cluster obtained (%d).' % (n_c)
            
    
    # Get contamination index.
    cont_index = g_c_i.cont_indx(backg_value[0], clust_rad[0], stars_in,
                                 stars_in_rjct)
    print 'Contamination index obtained (%0.2f).' % cont_index
    
    
    # Get cluster + field regions around the cluster's center.        
    flag_area_stronger, cluster_region, field_region = \
    g_r(x_center_bin, y_center_bin, width_bins, h_not_filt, clust_rad,
        H_manual, stars_in, stars_out)
    print 'Cluster + field stars regions obtained (%d).' % len(field_region)


    # Get physical cluster probability based on p_values distribution.
#    prob_cl_kde, p_vals_cl, p_vals_f, kde_cl_1d, kde_f_1d, x_kde, y_over = \
#    g_pv(flag_area_stronger, cluster_region, field_region, col1_data, \
#    mag_data, center_cl, clust_rad)
#    print 'Probability of physical cluster obtained (%0.2f).' % prob_cl_kde
#
#
#    # Get QQ plot for p-values distributions.
#    quantiles, r_squared, slope, intercept, ccc = g_qq(p_vals_cl, p_vals_f)
#    print 'QQ-plot obtained (R^2 = %0.2f)' % r_squared
    
    prob_cl_kde, p_vals_cl, p_vals_f, kde_cl_1d, kde_f_1d, x_kde, y_over, \
    quantiles, r_squared, slope, intercept, ccc = 0., [], [], [], [], [], [],\
    [], 0., 0., 0., 0.
    
################## Decontamination Algorithm Selection #######################
    # Apply decontamination algorithm only to stars with accepted photom
    # errors.
    
    # Uncomment to use KDE decontamination algorithm.
    field_reg_box = [] # Generated by the var box alg.
    # Apply algorithm if at least one equal-sized field region was found
    # around the cluster.
    if not(flag_area_stronger):
        print 'Applying KDE decontamination algorithm.'
        runs_fields_probs, kde_cl, kde_f = fdk(cluster_region, field_region,
                                             col1_data, mag_data,\
                                             center_cl, clust_rad)
    else:
        print 'WARNING: Decontamination algorithm was skipped.'
        # Skipping decontamination algorithm
        runs_fields_probs, kde_cl, kde_f = [], [], []
                                       
    # Uncomment to use variable box decontamination algorithm.
#    kde_cl, kde = [], [] # Generated by the KDE alg.
#    clus_reg_decont, field_reg_box, flag_area_stronger = \
#    fdvb(flag_area_stronger, cluster_region, field_region)

    # Uncomment to use Dias decontamination algorithm.
#    kde_cl, kde = [], [] # Generated by the KDE alg.
#    field_reg_box = [] # Generated by the var box alg.
#    clus_reg_decont = fdd(flag_area_stronger, cluster_region, center_cl,
#                          clust_rad)     
    
    # Uncomment to use random decontamination algorithm.
#    kde_cl, kde = [], [] # Generated by the KDE alg.
#    clus_reg_decont, field_reg_box, flag_area_stronger = \
#    fdr(flag_area_stronger, cluster_region, field_region, center_cl, clust_rad)

###############################################################################
    
    # Check if decont alg was applied.
    if flag_area_stronger:
        memb_prob_avrg_sort, clust_reg_prob_avrg = [], []
    else:
        # Average and sort all membership probabilities for each star and
        # store in list.
        memb_prob_avrg_sort, clust_reg_prob_avrg = m_p_a_s(cluster_region,\
        runs_fields_probs, n_c, center_cl, clust_rad[0])
        print 'Averaged probabilities for all runs.'
    
    
    # New name for cluster used in output data file.
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
        
        
    # Get manually fitted parameters for cluster, if these exist.
    cl_e_bv, cl_age, cl_feh, cl_dmod, iso_moved, zams_iso = g_i(mypath,
                                                                clust_name)
    print 'Isochrone obtained.'
    
        
    # Generate output subdir.   
    output_subdir = join(output_dir, sub_dir)
    # Check if subdir already exists, if not create it
    if not exists(output_subdir):
        makedirs(output_subdir)        
        
        
    # Make plots
    mp(output_subdir, clust_name, x_data, y_data, center_cl, cent_cl_err,
       x_center_bin, y_center_bin, h_filter, radii, backg_value, inner_ring,
       outer_ring, delta_backg, ring_density, poisson_error, clust_rad,
       cont_index, width_bins, delta_percentage, mag_data, col1_data,
       bright_end, popt_mag, popt_umag, pol_mag, popt_col1, popt_ucol1,
       pol_col1, mag_val_left, mag_val_right, col1_val_left, col1_val_right,
       use_errors_fit, k_prof, k_pr_err, flag_king_no_conver, stars_in,
       stars_out, stars_in_rjct,
       stars_out_rjct, stars_in_mag, stars_in_all_mag, n_c, flag_area_stronger,
       cluster_region, field_region,
       prob_cl_kde, p_vals_cl, p_vals_f, kde_cl_1d, kde_f_1d, x_kde, y_over,
       quantiles, r_squared, slope, intercept, ccc,
       clust_reg_prob_avrg, field_reg_box,
       kde_cl, kde_f, memb_prob_avrg_sort, iso_moved, zams_iso, cl_e_bv,
       cl_age, cl_feh, cl_dmod)
    print 'Plots created.'

    # Create data file with most probable members.
    c_m_f(output_dir, sub_dir, clust_name, memb_prob_avrg_sort)
    print 'Most probable members saved to file.'
 
   
    # Add cluster data and flags to output file
    a_d_o(sub_dir, output_dir, clust_name, center_cl, clust_rad, k_prof, 
          n_c_k, flag_king_no_conver, cont_index, n_c, prob_cl_kde, ccc,
          stars_in_mag,
          flag_center, flag_std_dev, flag_center_manual,
          flag_radius_manual, flag_errors_manual, flag_bin_count,
          flag_delta_total, flag_not_stable, flag_rad_500, flag_delta,
          flag_delta_points, flag_num_memb_low, flag_no_memb)
    print 'Data added to output file.'
    

    # Move file to 'done' dir.
#    dst_dir = join(mypath3, sub_dir)
#    # If the sub-dir doesn't exist, create it before moving the file.
#    if not exists(dst_dir):
#        mkdir(dst_dir)
#    shutil.move(join(mypath2, sub_dir, myfile), dst_dir)
#    # If sub-dir is empty, remove it.
#    try:
#        rmdir(join(mypath2, sub_dir))
#    except OSError as ex:
#        # Sub-dir not empty, skip.
#        pass
#    print "Data file moved to 'done' folder."
   

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
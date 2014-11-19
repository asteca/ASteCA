# -*- coding: utf-8 -*-

import time
import gc  # Garbage collector.
# Import files with defined functions.
from functions.get_in_params import get_in_params as gip
from functions.get_names_paths import names_paths as n_p
from functions.get_data_semi import get_semi as g_s
from functions.get_data import get_data as gd
from functions.trim_frame import trim_frame as t_f
from functions.get_2d_histo import get_2d_histo as g2dh
from functions.get_center import get_center as g_c
from functions.get_field_dens import field_dens as gfd
from functions.get_dens_prof import get_dens_prof as gdp
from functions.get_radius import get_clust_rad as gcr
from functions.get_king_prof import get_king_profile as gkp
from functions.err_accpt_rejct import err_accpt_rejct as ear
from functions.get_in_out import get_in_out as gio
from functions.get_integ_mag import integ_mag as g_i_m
from functions.get_members_number import get_memb_num as g_m_n
from functions.get_cont_index import cont_indx as g_c_i
from functions.get_regions import get_regions as g_r
from functions.decont_algor_bys import bys_da as dab
from functions.get_lf import lf
from functions.reduce_membership import red_memb as rm
from functions.synth_cl_err import synth_clust_err as sce
from functions.best_fit_synth_cl import best_fit as bfsc
from functions.make_plots import make_plots as mp
from functions.create_out_data_file import create_out_data_file as c_o_d_f
from functions.add_data_output import add_data_output as a_d_o
from functions.cl_members_file import cluster_members_file as c_m_f
from functions.done_move import done_move as dm


def asteca_funcs(mypath, cl_file, ip_list, R_in_place):
    '''
    Container that holds the calls to all the functions.
    '''

    # Start timing this loop.
    start = time.time()

    # Read input parameters from params_input.dat file.
    mode, done_dir, gd_params, gh_params, cr_params, kp_flag, \
    im_flag, er_params, fr_number, pv_params, da_params, ps_params, bf_params,\
    sc_params, ga_params, rm_params, pl_params, flag_move_file, \
    axes_params = gip(mypath)

    # Define system of coordinates used.
    px_deg = gd_params[-1]
    coord_lst = ['px', 'x', 'y'] if px_deg == 'px' else ['deg', 'ra', 'dec']
    coord, x_name, y_name = coord_lst

    # Get file names and paths.
    clust_name, data_file, memb_file, output_dir, output_subdir, \
    memb_file_out, write_name = n_p(mypath, cl_file)
    print 'Analizing cluster {}/{} ({}).'.format(cl_file[0], clust_name, mode)

    # Get data from semi-data input file.
    mode, semi_return = g_s(clust_name, mode)

    # Get cluster's photometric data from file.
    phot_data = gd(data_file, gd_params)
    # If Manual mode is set, display frame and ask if it should be trimmed.
    phot_data = t_f(phot_data, mode)
    # Unpack coordinates, magnitude and color.
    x_data, y_data, mag_data, col1_data = phot_data[1], phot_data[2], \
    phot_data[3], phot_data[5]

    # Obtain 2D histograms for the observed frame using several bin widths.
    hist_lst = g2dh(x_data, y_data, gh_params)
    bin_width = hist_lst[-1]
    print "Frame's 2D histogram obtained"

    # Get cluster's center coordinates and errors.
    center_params = g_c(x_data, y_data, mag_data, hist_lst, mode,
        semi_return, coord_lst)
    # Unpack values from list.
    cent_bin, kde_center = center_params[0], center_params[1]

    # Get density profile
    rdp_params = gdp(hist_lst, cent_bin)
    radii, ring_density = rdp_params[:2]
    print 'Radial density profile (RDP) calculated.'

    # Get field density value in stars/px^2.
    field_dens = gfd(ring_density)
    print 'Field density calculated ({:.1E} stars/{c}^2).'.format(field_dens,
        c=coord)

    # Get cluster radius
    radius_params = gcr(phot_data, field_dens, cr_params, center_params,
        rdp_params, semi_return, mode, bin_width, coord_lst)
    clust_rad = radius_params[0]

    # Get King profiles based on the density profiles.
    kp_params = gkp(kp_flag, clust_rad, field_dens, radii, ring_density,
        coord_lst)

    # Get approximate number of cluster's members.
    n_memb, flag_num_memb_low, a_clust, n_clust = g_m_n(field_dens, clust_rad,
        rdp_params, bin_width)
    print ("Approximate number of members in cluster obtained "
        "({:.0f}).".format(n_memb))

    # Get contamination index.
    cont_index = g_c_i(field_dens, a_clust, n_clust)
    print 'Contamination index obtained ({:.2f}).'.format(cont_index)

    # Accept and reject stars based on their errors.
    acpt_stars, rjct_stars, err_plot, err_flags, err_pck, er_params = \
    ear(phot_data, axes_params, er_params, mode, semi_return)

    # Get stars in and out of cluster's radius.
    cl_region, stars_out, stars_in_rjct, stars_out_rjct = gio(kde_center,
        clust_rad, acpt_stars, rjct_stars)
    print "Stars separated in/out of cluster's boundaries."

    # Field regions around the cluster's center.
    flag_area_stronger, field_region = g_r(hist_lst, cent_bin,
        clust_rad, stars_out, fr_number)

    # Get the luminosity function and completeness level for each magnitude
    # bin. The completeness will be used by the isochrone/synthetic cluster
    # fitting algorithm.
    lum_func, completeness = lf(flag_area_stronger, mag_data, cl_region,
        field_region)
    print 'LF and Completeness magnitude levels obtained.'

    # Calculate integrated magnitude.
    integr_return = g_i_m(im_flag, cl_region, field_region, axes_params,
        flag_area_stronger)

    # Get physical cluster probability based on p_values distribution.
    if R_in_place:
        from functions.get_p_value import get_pval as g_pv
        pval_test_params, flag_pval_test = g_pv(cl_region, field_region,
            col1_data, mag_data, pv_params, flag_area_stronger)
    else:
        print 'Skipping KDE p-value function.'
        flag_pval_test, pval_test_params = False, [-1., [], [], [], [], [], []]

    # Apply decontamination algorithm if at least one equal-sized field region
    # was found around the cluster.
    print 'Applying decontamination algorithm.'
    decont_algor_return = dab(flag_area_stronger, cl_region, field_region,
                            memb_file, da_params)
    memb_prob_avrg_sort = decont_algor_return[0]

    # Create data file with membership probabilities.
    c_m_f(memb_file_out, memb_prob_avrg_sort)
    print 'Membership probabilities saved to file.'

    # Reduce number of stars in cluster according to a lower membership
    # probability or magnitude limit.
    red_return = rm(decont_algor_return, bf_params, rm_params)
    red_memb_prob = red_return[0]

    # Obtain exponential error function parameters to use by the
    # synthetic cluster creation function.
    err_lst = sce(phot_data, err_pck, bf_params, da_params)
    # Obtain best fitting parameters for cluster.
    bf_return = bfsc(err_lst, red_memb_prob, completeness, ip_list, bf_params,
        sc_params, ga_params, ps_params)

    # Create output data file in /output dir if it doesn't exist.
    out_file_name = c_o_d_f(output_dir)

    # Add cluster data and flags to output file
    a_d_o(out_file_name, write_name, center_params,
        radius_params, kp_params, cont_index, n_memb, pval_test_params[0],
        integr_return, axes_params, err_flags, flag_num_memb_low, bf_return)
    print 'Data added to output file.'

    # Make plots
    if pl_params[0]:
        mp(output_subdir, clust_name, x_data, y_data, gd_params,
            bin_width, center_params, rdp_params,
            field_dens, radius_params, cont_index, mag_data, col1_data,
            err_plot, err_flags, kp_params, cl_region, stars_out,
            stars_in_rjct, stars_out_rjct, integr_return, n_memb,
            flag_area_stronger, field_region, flag_pval_test,
            pval_test_params, memb_prob_avrg_sort, lum_func, completeness,
            ip_list, da_params, bf_params, red_return, err_lst, bf_return,
            ga_params, er_params, axes_params, pl_params)
        print 'Plots created.'

    # Move file to 'done' dir if flag is set.
    dm(flag_move_file, mypath, cl_file, done_dir, data_file, memb_file)

    elapsed = time.time() - start
    m, s = divmod(elapsed, 60)
    print 'End of analysis for {} in {:.0f}m {:.0f}s.\n'.format(clust_name,
        m, s)

    # Force the Garbage Collector to release unreferenced memory.
    gc.collect()

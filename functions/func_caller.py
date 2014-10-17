# -*- coding: utf-8 -*-

import time
import gc  # Garbage collector.
import get_in_params as g
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


def asteca_funcs(mypath, cl_file, phot_params, ip_list, R_in_place):
    '''
    Container that holds the calls to all the functions.
    '''

    # Start timing this loop.
    start = time.time()

    # Get file names and paths.
    clust_name, data_file, memb_file, output_dir, output_subdir, \
    memb_file_out, write_name = n_p(mypath, cl_file)
    print 'Analizing cluster {} ({}).'.format(clust_name, g.mode)

    # Get data from semi-data input file.
    semi_return = g_s(clust_name)

    # Get cluster's photometric data from file.
    id_coords, phot_data = gd(data_file, phot_params)
    # If Manual mode is set, display frame and ask if it should be trimmed.
    id_coords, phot_data = t_f(id_coords, phot_data)

    # Obtain 2D histograms for the observed frame using several bin widths.
    hist_lst = g2dh(id_coords)
    bin_width = hist_lst[-1]
    print "Frame's 2D histogram obtained"

    # Get cluster's center coordinates and errors.
    center_params = g_c(id_coords, phot_data, hist_lst, semi_return)
    # Unpack values from list.
    cent_bin, kde_center = center_params[0], center_params[1]

    # Get density profile
    rdp_params = gdp(hist_lst, cent_bin)
    radii, ring_density = rdp_params[:2]
    print 'Radial density profile (RDP) calculated.'

    # Get field density value in stars/px^2.
    field_dens = gfd(ring_density)
    print 'Field density calculated ({:.1E} stars/{}^2).'.format(field_dens,
        g.gd_params[0][-1])

    # Get cluster radius
    radius_params = gcr(phot_data, field_dens, center_params, semi_return,
        bin_width, rdp_params)
    clust_rad = radius_params[0]

    # Get King profiles based on the density profiles.
    kp_params = gkp(clust_rad, field_dens, radii, ring_density)

    # Get approximate number of cluster's members.
    n_memb, flag_num_memb_low, a_clust, n_clust = g_m_n(field_dens, clust_rad,
        rdp_params, bin_width)
    print ("Approximate number of members in cluster obtained "
        "({:.0f}).".format(n_memb))

    # Get contamination index.
    cont_index = g_c_i(field_dens, a_clust, n_clust)
    print 'Contamination index obtained ({:.2f}).'.format(cont_index)

    # Accept and reject stars based on their errors.
    acpt_stars, rjct_stars, err_plot, err_flags, err_pck = ear(id_coords,
        phot_data, semi_return, phot_params[-1])

    # Get stars in and out of cluster's radius.
    cl_region, stars_out, stars_in_rjct, stars_out_rjct = gio(kde_center,
        clust_rad, acpt_stars, rjct_stars)
    print "Stars separated in/out of cluster's boundaries."

    # Field regions around the cluster's center.
    flag_area_stronger, field_regions = g_r(hist_lst, cent_bin, clust_rad,
        stars_out)

    # Get the luminosity function and completeness level for each magnitude
    # bin. The completeness will be used by the isochrone/synthetic cluster
    # fitting algorithm.
    lum_func, completeness = lf(flag_area_stronger, phot_data, cl_region,
        field_regions)
    print 'LF and Completeness magnitude levels obtained.'

    # Calculate integrated magnitude.
    integr_return = g_i_m(cl_region, field_regions, flag_area_stronger)

    # Get physical cluster probability based on p_values distribution.
    if R_in_place:
        from functions.get_p_value import get_pval as g_pv
        pval_test_params, flag_pval_test = g_pv(cl_region, field_regions,
            phot_data, flag_area_stronger)
    else:
        print 'Skipping KDE p-value function.'
        flag_pval_test, pval_test_params = False, [-1., [], [], [], [], [], []]

    # Apply decontamination algorithm if at least one equal-sized field region
    # was found around the cluster.
    decont_algor_return = dab(flag_area_stronger, cl_region, field_regions,
        memb_file)
    memb_prob_avrg_sort = decont_algor_return[0]

    # Create data file with membership probabilities.
    c_m_f(memb_file_out, memb_prob_avrg_sort)
    print 'Membership probabilities saved to file.'

    # Reduce number of stars in cluster according to a lower membership
    # probability or magnitude limit.
    red_return = rm(decont_algor_return)
    red_memb_prob = red_return[0]

    # Obtain exponential error function parameters to use by the
    # synthetic cluster creation function.
    err_lst = sce(phot_data, err_pck)
    # Obtain best fitting parameters for cluster.
    bf_return = bfsc(err_lst, red_memb_prob, completeness, ip_list)

    # Create output data file in /output dir if it doesn't exist.
    out_file_name = c_o_d_f(output_dir)

    # Add cluster data and flags to output file
    a_d_o(out_file_name, write_name, center_params, radius_params, kp_params,
        cont_index, n_memb, pval_test_params[0], integr_return, err_flags,
        flag_num_memb_low, bf_return)
    print 'Data added to output file.'

    # Make plots
    if g.pl_params[0]:
        mp(output_subdir, clust_name, phot_params, id_coords, phot_data,
            bin_width, center_params, rdp_params, field_dens, radius_params,
            cont_index, err_plot, err_flags, kp_params, cl_region, stars_out,
            stars_in_rjct, stars_out_rjct, integr_return, n_memb,
            flag_area_stronger, field_regions, flag_pval_test,
            pval_test_params, memb_prob_avrg_sort, lum_func, completeness,
            ip_list, red_return, err_lst, bf_return)
        print 'Plots created.'

    # Move file to 'done' dir if flag is set.
    dm(mypath, cl_file, data_file, memb_file)

    elapsed = time.time() - start
    m, s = divmod(elapsed, 60)
    print 'End of analysis for {} in {:.0f}m {:.0f}s.\n'.format(clust_name,
        m, s)

    # Force the Garbage Collector to release unreferenced memory.
    gc.collect()

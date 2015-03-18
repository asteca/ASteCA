# -*- coding: utf-8 -*-

import time
import gc  # Garbage collector.
#
import _in.get_in_params as g
from functions._in.get_names_paths import names_paths as n_p
from functions._in.get_data_semi import get_semi as g_s
from functions._in.get_data import get_data as gd
#
from functions.structure.trim_frame import trim_frame as t_f
from functions.structure.get_2d_histo import get_2d_histo as g2dh
from functions.structure.get_center import get_center as g_c
from functions.structure.get_field_dens import field_dens as gfd
from functions.structure.get_dens_prof import get_dens_prof as gdp
from functions.structure.get_radius import get_clust_rad as gcr
from functions.structure.get_cluster_area import get_cl_area as g_a
from functions.structure.get_king_prof import get_king_profile as gkp
from functions.structure.get_in_out import get_in_out as gio
from functions.structure.get_regions import get_regions as g_r
#
from functions.errors.err_accpt_rejct import err_accpt_rejct as ear
#
from functions.phot_analysis.get_integ_mag import integ_mag as g_i_m
from functions.phot_analysis.get_members_number import get_memb_num as g_m_n
from functions.phot_analysis.get_cont_index import cont_indx as g_c_i
from functions.phot_analysis.decont_algor_bys import bys_da as dab
from functions.phot_analysis.get_members_param import mp_members as m_m
from functions.phot_analysis.get_lf import lf
#
from functions.best_fit.reduce_membership import red_memb as rm
from functions.best_fit.synth_cl_err import synth_clust_err as sce
from functions.best_fit.best_fit_synth_cl import best_fit as bfsc
#
from functions.out.make_plots import make_plots as mp
from functions.out.create_out_data_file import create_out_data_file as c_o_d_f
from functions.out.top_tiers import get_top_tiers as g_t_t
from functions.out.add_data_output import add_data_output as a_d_o
from functions.out.cl_members_file import cluster_members_file as c_m_f
from functions.out.done_move import done_move as dm


def asteca_funcs(cl_file, ip_list, R_in_place):
    '''
    Container that holds the calls to all the functions.
    '''

    # Start timing this loop.
    start = time.time()

    # Get file names and paths.
    clust_name, data_file, memb_file, output_dir, output_subdir, dst_dir,\
    memb_file_out, write_name = n_p(cl_file)
    print 'Analizing cluster {} ({}).'.format(clust_name, g.mode)

    # Get data from semi-data input file.
    semi_return = g_s(clust_name)

    # Get cluster's photometric data from file.
    phot_data = gd(data_file)
    # If Manual mode is set, display frame and ask if it should be trimmed.
    phot_data = t_f(phot_data)
    # Unpack coordinates, magnitude and color.
    x_data, y_data, mag_data, col1_data = phot_data[1], phot_data[2], \
    phot_data[3], phot_data[5]

    # Obtain 2D histograms for the observed frame using several bin widths.
    hist_lst = g2dh(x_data, y_data)
    bin_width = hist_lst[-1]
    print "Frame's 2D histogram obtained"

    # Get cluster's center coordinates and errors.
    center_params = g_c(x_data, y_data, mag_data, hist_lst, semi_return)
    # Unpack values from list.
    cent_bin, kde_center = center_params[0], center_params[1]

    # Get density profile
    rdp_params = gdp(hist_lst, cent_bin)
    radii, ring_density = rdp_params[:2]
    print 'Radial density profile (RDP) calculated.'

    # Get field density value in stars/px^2.
    field_dens = gfd(ring_density)

    # Get cluster radius
    radius_params = gcr(phot_data, field_dens, center_params, rdp_params,
        semi_return, bin_width)
    clust_rad = radius_params[0]

    # Get King profiles based on the density profiles.
    kp_params = gkp(clust_rad, field_dens, radii, ring_density)

    # Accept and reject stars based on their errors.
    acpt_stars, rjct_stars, err_plot, err_flags, err_pck = ear(phot_data,
        semi_return)

    # Get stars in and out of cluster's radius.
    cl_region, stars_out, stars_in_rjct, stars_out_rjct = gio(kde_center,
        clust_rad, acpt_stars, rjct_stars)
    print "Stars separated in/out of cluster's boundaries."
    # Number of stars inside the cluster region (including stars with rejected
    # photometric errors)
    n_clust = len(cl_region) + len(stars_in_rjct)

    # Get cluster's area.
    cl_area, frac_cl_area = g_a(kde_center, clust_rad, x_data, y_data,
        rdp_params, bin_width)
    print "Area of cluster obtained."

    # Get approximate number of cluster's members.
    n_memb, flag_num_memb_low = g_m_n(n_clust, cl_area, field_dens)
    print ("Approximate number of members in cluster obtained "
        "({:.0f}).".format(n_memb))

    # Get contamination index.
    cont_index = g_c_i(field_dens, cl_area, n_clust)
    print 'Contamination index obtained ({:.2f}).'.format(cont_index)

    # Field regions around the cluster's center.
    flag_area_stronger, field_region = g_r(hist_lst, cent_bin, clust_rad,
        cl_area, stars_out)

    # Get the luminosity function and completeness level for each magnitude
    # bin. The completeness will be used by the isochrone/synthetic cluster
    # fitting algorithm.
    lum_func, completeness = lf(flag_area_stronger, mag_data, cl_region,
        field_region)
    print 'LF and Completeness magnitude levels obtained.'

    # Calculate integrated magnitude.
    integr_return = g_i_m(cl_region, field_region, flag_area_stronger)

    # Get physical cluster probability based on p_values distribution.
    if R_in_place:
        from phot_analysis.get_p_value import get_pval as g_pv
        pval_test_params, flag_pval_test = g_pv(cl_region, field_region,
            col1_data, mag_data, flag_area_stronger)
    else:
        print 'Missing package. Skipping KDE p-value test for cluster.'
        flag_pval_test, pval_test_params = False, [-1.]

    # Apply decontamination algorithm if at least one equal-sized field region
    # was found around the cluster.
    decont_algor_return = dab(flag_area_stronger, cl_region, field_region,
        memb_file)
    memb_prob_avrg_sort = decont_algor_return[0]

    # Create data file with membership probabilities.
    c_m_f(memb_file_out, memb_prob_avrg_sort)
    print 'Membership probabilities saved to file.'

    # Obtain members parameter.
    memb_par, n_memb_da = m_m(n_memb, decont_algor_return)

    # Reduce number of stars in cluster according to a lower membership
    # probability or magnitude limit.
    red_return = rm(n_memb, decont_algor_return)

    # Obtain exponential error function parameters to use by the
    # synthetic cluster creation function.
    err_lst = sce(phot_data, err_pck)
    # Obtain best fitting parameters for cluster.
    bf_return = bfsc(err_lst, red_return[0], completeness, ip_list)

    # Create output data file in /output dir if it doesn't exist.
    out_file_name = c_o_d_f(output_dir)

    # Add cluster data and flags to output file
    a_d_o(out_file_name, write_name, center_params, radius_params, kp_params,
        cont_index, n_memb, memb_par, n_memb_da, frac_cl_area,
        pval_test_params[0], integr_return, err_flags, flag_num_memb_low,
        bf_return)
    print 'Data added to output file.'

    # Output top tiers models if best fit parameters were obtained.
    g_t_t(clust_name, output_subdir, mag_data, col1_data, ip_list, err_lst,
        completeness, bf_return)

    # Make plots
    if g.pl_params[0]:
        mp(output_subdir, clust_name, x_data, y_data,
            bin_width, center_params, rdp_params,
            field_dens, radius_params, cont_index, mag_data, col1_data,
            err_plot, err_flags, kp_params, cl_region, stars_out,
            stars_in_rjct, stars_out_rjct, integr_return, n_memb, n_memb_da,
            flag_area_stronger, field_region, flag_pval_test,
            pval_test_params, decont_algor_return, lum_func, completeness,
            ip_list, red_return, err_lst, bf_return)
        print 'Plots created.'

    # Move file to 'done' dir if flag is set.
    dm(dst_dir, data_file, memb_file)

    elapsed = time.time() - start
    m, s = divmod(elapsed, 60)
    print 'End of analysis for {} in {:.0f}m {:.0f}s.\n'.format(clust_name,
        m, s)

    # Force the Garbage Collector to release unreferenced memory.
    gc.collect()

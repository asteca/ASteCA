
import time
import gc  # Garbage collector.
#
import inp.input_params as g
from packages.inp import names_paths
from packages.inp import get_data_semi
from packages.inp import get_data
#
from packages.structure import trim_frame
from packages.structure import histo_2d
from packages.structure import center
from packages.structure import radial_dens_prof
from packages.structure import field_density
from packages.structure import radius
from packages.structure import king_profile
from packages.errors import err_accpt_rejct
from packages.structure import stars_in_out_cl_reg
from packages.structure import cluster_area
from packages.phot_analysis import members_number
from packages.phot_analysis import contamination_index
from packages.structure import field_regions
from packages.phot_analysis import luminosity_func
from packages.phot_analysis import integrated_mag
# from phot_analysis import kde_pvalue
from packages.decont_algors import bayesian_da
from packages.phot_analysis import members_N_compare
from packages.decont_algors import membership_removal
from packages.out import cluster_members_file
from packages.best_fit import synth_cl_err
from packages.best_fit import best_fit_synth_cl
from packages.out import synth_cl_file
from packages.out import create_out_data_file
from packages.out import add_data_output
from packages.out import top_tiers
from packages.out import make_plots
from packages.out import done_move


def main(cl_file, ip_list, R_in_place):
    '''
    Container that holds the calls to all the modules and functions.
    '''

    # Start timing this loop.
    start = time.time()

    # Get file names and paths.
    clust_name, data_file, memb_file, output_dir, output_subdir, dst_dir,\
        memb_file_out, synth_file_out, write_name = names_paths.main(cl_file)
    print('Analyzing cluster {} ({} mode).'.format(clust_name, g.mode))

    # Get data from semi-data input file.
    semi_return = get_data_semi.main(clust_name)

    # Get cluster's photometric data from file.
    phot_data = get_data.main(data_file)
    # If Manual mode is set, display frame and ask if it should be trimmed.
    phot_data = trim_frame.main(phot_data)
    # Unpack coordinates, magnitude and color.
    x_data, y_data, mag_data, col1_data = phot_data[1], phot_data[2], \
        phot_data[3], phot_data[5]

    # Obtain 2D histograms for the observed frame using several bin widths.
    hist_lst = histo_2d.main(x_data, y_data)
    bin_width = hist_lst[-1]

    # Get cluster's center coordinates and errors.
    center_params = center.main(x_data, y_data, mag_data, hist_lst,
                                semi_return)
    # Unpack values from list.
    cent_bin, kde_center = center_params[0], center_params[1]

    # Get density profile
    rdp_params = radial_dens_prof.main(hist_lst, cent_bin)
    radii, rdp_points, square_rings, rdp_length = rdp_params[:2] + \
        rdp_params[3:]

    # Get field density value in stars/px^2.
    field_dens = field_density.main(rdp_points)

    # Get cluster radius
    radius_params = radius.main(phot_data, field_dens, center_params,
                                rdp_params, semi_return, bin_width)
    clust_rad = radius_params[0]

    # Get King profiles based on the density profiles.
    kp_params = king_profile.main(clust_rad, field_dens, radii, rdp_points)

    # Accept and reject stars based on their errors.
    acpt_stars, rjct_stars, err_plot, err_flags, err_pck =\
        err_accpt_rejct.main(phot_data, semi_return)

    # Get stars in and out of cluster's radius.
    cl_region, stars_out, stars_in_rjct, stars_out_rjct =\
        stars_in_out_cl_reg.main(kde_center, clust_rad, acpt_stars, rjct_stars)
    # Number of stars inside the cluster region (including stars with rejected
    # photometric errors)
    n_clust = len(cl_region) + len(stars_in_rjct)

    # Get cluster's area.
    cl_area, frac_cl_area = cluster_area.main(
        kde_center, clust_rad, x_data, y_data, square_rings, bin_width)

    # Get approximate number of cluster's members.
    n_memb, flag_num_memb_low = members_number.main(
        n_clust, cl_area, field_dens, clust_rad, rdp_length)

    # Get contamination index.
    cont_index = contamination_index.main(n_clust, cl_area, field_dens,
                                          clust_rad, rdp_length)

    # Field regions around the cluster's center.
    flag_no_fl_regs, field_region = field_regions.main(
        semi_return, hist_lst, cent_bin, clust_rad, cl_area, stars_out)

    # Get the luminosity function and completeness level for each magnitude
    # bin. The completeness will be used by the isochrone/synthetic cluster
    # fitting algorithm.
    lum_func, completeness = luminosity_func.main(
        flag_no_fl_regs, mag_data, cl_region, field_region)

    # Calculate integrated magnitude.
    integr_return = integrated_mag.main(cl_region, field_region,
                                        flag_no_fl_regs)

    # Get physical cluster probability based on p_values distribution.
    if R_in_place:
        from phot_analysis import kde_pvalue
        pval_test_params, flag_pval_test = kde_pvalue.main(
            cl_region, field_region, col1_data, mag_data, flag_no_fl_regs)
    else:
        print('Missing package. Skipping KDE p-value test for cluster.')
        flag_pval_test, pval_test_params = False, [-1.]

    # Apply decontamination algorithm if at least one equal-sized field region
    # was found around the cluster.
    bayes_da_return = bayesian_da.main(flag_no_fl_regs, cl_region,
                                       field_region, memb_file)

    # Obtain members parameter.
    memb_par, n_memb_da, flag_memb_par = members_N_compare.main(
        n_memb, bayes_da_return)

    # Reduce number of stars in cluster according to a lower membership
    # probability or magnitude limit.
    memb_remove = membership_removal.main(n_memb, flag_no_fl_regs,
                                          bayes_da_return, field_region)

    # Create data file with membership probabilities.
    cluster_members_file.main(memb_file_out, memb_remove)

    # Obtain exponential error function parameters to use by the
    # synthetic cluster creation function.
    err_lst = synth_cl_err.main(phot_data, err_pck)
    # Obtain best fitting parameters for cluster.
    bf_return = best_fit_synth_cl.main(err_lst, memb_remove[0], completeness,
                                       ip_list)

    # Create output synthetic cluster file if one was found
    synth_cl_file.main(bf_return[3], synth_file_out)

    # Create output data file in /output dir if it doesn't exist.
    out_file_name = create_out_data_file.main(output_dir)

    # Add cluster data and flags to output file
    add_data_output.main(
        out_file_name, write_name, center_params, radius_params, kp_params,
        cont_index, n_memb, memb_par, n_memb_da, flag_memb_par, frac_cl_area,
        pval_test_params[0], integr_return, err_flags, flag_num_memb_low,
        bf_return)

    # Output top tiers models if best fit parameters were obtained.
    top_tiers.main(clust_name, output_subdir, mag_data, col1_data, ip_list,
                   err_lst, completeness, bf_return)

    # Make plots
    if g.pl_params[0]:
        make_plots.main(
            output_subdir, clust_name, x_data, y_data,
            bin_width, center_params, rdp_params,
            field_dens, radius_params, cont_index, mag_data, col1_data,
            err_plot, err_flags, kp_params, cl_region, stars_out,
            stars_in_rjct, stars_out_rjct, integr_return, n_memb, n_memb_da,
            flag_no_fl_regs, field_region, flag_pval_test,
            pval_test_params, bayes_da_return, lum_func, completeness,
            ip_list, memb_remove, err_lst, bf_return)

    # Move file to 'done' dir if flag is set.
    done_move.main(dst_dir, data_file, memb_file)

    elapsed = time.time() - start
    m, s = divmod(elapsed, 60)
    print('End of analysis for {} in {:.0f}m {:.0f}s.\n'.format(
        clust_name, m, s))

    # Force the Garbage Collector to release unreferenced memory.
    gc.collect()

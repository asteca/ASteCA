"""
@author: gabriel
"""

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from os.path import join
import warnings
# Custom functions.
from functions import __version__
from .._in import get_in_params as g
from ..errors import error_round as err_r
import mp_structure
import mp_phot_analysis
import mp_decont_algor
import mp_best_fit
import prep_plots


#############################################################
## Timer function: http://stackoverflow.com/a/21860100/1391441
#from contextlib import contextmanager
#import time
#@contextmanager
#def timeblock(label):
    #start = time.clock()
    #try:
        #yield
    #finally:
        #end = time.clock()
        #print ('{} elapsed: {}'.format(label, end - start))
#############################################################


def make_plots(output_subdir, clust_name, x_data, y_data,
    bin_width, center_params, rdp_params, field_dens, radius_params,
    cont_index, mag_data, col_data, err_plot, err_flags, kp_params,
    cl_region, stars_out, stars_in_rjct, stars_out_rjct, integr_return, n_memb,
    n_memb_da, flag_no_fl_regs, field_regions, flag_pval_test,
    pval_test_params, decont_algor_return, lum_func, completeness, ip_list,
    red_return, err_lst, bf_return):
    '''
    Make all plots.
    '''

    # Unpack params.
    # Parameters from get_center function.
    cent_bin, center_cl, e_cent, approx_cents, st_dev_lst, hist_2d_g, \
    kde_pl = center_params[:7]
    # RDP params.
    radii, ring_density, poisson_error = rdp_params[:3]
    # Parameters from get_radius function.
    clust_rad, e_rad = radius_params[:2]
    # Luminosity functions.
    x_cl, y_cl, x_fl, y_fl = lum_func
    # Best fitting process results.
    isoch_fit_params, isoch_fit_errors, shift_isoch, synth_clst, syn_b_edges = \
    bf_return

    # Plot all outputs
    # figsize(x1, y1), GridSpec(y2, x2) --> To have square plots: x1/x2 =
    # y1/y2 = 2.5
    fig = plt.figure(figsize=(30, 25))  # create the top-level container
    gs = gridspec.GridSpec(10, 12)      # create a GridSpec object
    # Add version number to top left.
    ver = '[ASteCA ' + __version__ + ']'
    x_coord = 0.957 - (len(__version__) - 6) * 0.001
    plt.figtext(x_coord, .988, ver, fontsize=9, color='#585858')

    # Obtain plotting parameters and data.
    x_min, x_max, y_min, y_max = prep_plots.frame_max_min(x_data, y_data)
    coord, x_name, y_name = prep_plots.coord_syst()
    x_zmin, x_zmax, y_zmin, y_zmax = prep_plots.frame_zoomed(x_min, x_max,
        y_min, y_max, center_cl, clust_rad)
    x_ax, y_ax, x_ax0, y_axis = prep_plots.ax_names()
    phot_x, phot_y = prep_plots.ax_data(mag_data, col_data)
    x_max_cmd, x_min_cmd, y_min_cmd, y_max_cmd = prep_plots.diag_limits(y_axis,
        phot_x, phot_y)
    x_data_z, y_data_z, mag_data_z, stars_f_rjct, stars_f_acpt = \
    prep_plots.separate_stars(x_data, y_data, mag_data, x_zmin, x_zmax,
        y_zmin, y_zmax, stars_out_rjct, field_regions)
    st_sizes_arr = prep_plots.star_size(mag_data)
    st_sizes_arr_z = prep_plots.star_size(mag_data_z)
    f_sz_pt = prep_plots.phot_diag_st_size(len(stars_f_acpt[0]))
    cl_sz_pt = prep_plots.phot_diag_st_size(len(cl_region))
    v_min_mp, v_max_mp, plot_colorbar, chart_fit_inv, chart_no_fit_inv, \
    out_clust_rad, diag_fit_inv, diag_no_fit_inv, err_bar = \
        prep_plots.da_plots(center_cl, clust_rad, stars_out, x_zmin, x_zmax,
        y_zmin, y_zmax, x_max_cmd, col_data, err_lst, red_return)

    #
    # Structure plots.
    arglist = [
        [gs, fig, x_name, y_name, coord, cent_bin, clust_rad, bin_width,
            err_r, hist_2d_g],
        [gs, x_name, y_name, coord, x_min, x_max, y_min, y_max, approx_cents,
            bin_width, st_dev_lst],
        [gs, fig, x_name, y_name, coord, x_min, x_max, y_min, y_max,
            center_cl, clust_rad, e_cent, kp_params, mag_data,
            x_data, y_data, st_sizes_arr],
        [gs, radii, ring_density, field_dens, coord, clust_name, kp_params,
            clust_rad, e_rad, poisson_error, bin_width],
        [gs, fig, x_name, y_name, coord, x_zmin, x_zmax, y_zmin, y_zmax,
            cont_index, kde_pl, x_data_z, y_data_z,
            st_sizes_arr_z, center_cl, clust_rad],
        [gs, fig, x_name, y_name, coord, x_min, x_max, y_min, y_max,
            center_cl, clust_rad, field_regions, cl_region,
            flag_no_fl_regs]
        ]
    for n, args in enumerate(arglist):
        #with timeblock("{}".format(n)):
        mp_structure.plot(n, *args)

    #
    # Photometric analysis plots.

    arglist = [
        [gs, fig, 'up', x_ax, y_ax, mag_data, err_plot, err_flags,
            cl_region, stars_in_rjct, stars_out, stars_out_rjct],
        [gs, fig, 'low', x_ax, y_ax, mag_data, err_plot, err_flags,
            cl_region, stars_in_rjct, stars_out, stars_out_rjct],
        [gs, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd, x_ax, y_ax,
            stars_f_rjct, stars_f_acpt, f_sz_pt],
        [gs, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd, x_ax, y_ax,
            stars_in_rjct, cl_region, n_memb, cl_sz_pt],
        [gs, mag_data, y_ax, x_cl, y_cl, flag_no_fl_regs, x_fl, y_fl,
            completeness],
        [gs, integr_return, y_ax, x_ax0, flag_no_fl_regs],
        [gs, flag_pval_test, pval_test_params]
        ]
    for n, args in enumerate(arglist):
        mp_phot_analysis.plot(n, *args)

    #
    # Decontamination algorithm plots.
    flag_decont_skip = decont_algor_return[1]
    mode_red_memb = g.rm_params[0]
    bf_flag = g.bf_params[0]

    # If the DA and the best fit functions were skipped and the reduced
    # membership mode is any mode but 'local', do not plot.
    if flag_decont_skip and bf_flag is False and mode_red_memb != 'local':
        pass
    else:
        arglist = [
            [gs, n_memb_da, red_return, decont_algor_return],
            [gs, fig, x_name, y_name, coord, x_zmin, x_zmax, y_zmin, y_zmax,
                center_cl, clust_rad, field_dens, flag_decont_skip,
                v_min_mp, v_max_mp, chart_fit_inv, chart_no_fit_inv,
                out_clust_rad]
            ]
        for n, args in enumerate(arglist):
            mp_decont_algor.plot(n, *args)

        # This function is called separately since we need to retrieve some
        # information from it to plot that #$%&! colorbar.
        try:
            sca, trans = mp_decont_algor.pl_mps_phot_diag(gs, fig,
                x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd, x_ax, y_ax,
                v_min_mp, v_max_mp, red_return, diag_fit_inv, diag_no_fit_inv,
                shift_isoch, err_bar)
        except:
            #import traceback
            #print traceback.format_exc()
            print("  WARNING: error when plotting MPs on cluster's "
            "photom diagram.")

    #
    # Best fit plots.
    if bf_flag:
        arglist = [
            [gs, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd, x_ax, y_ax,
                synth_clst, syn_b_edges, isoch_fit_params[0], isoch_fit_errors,
                shift_isoch]
            ]
        for n, args in enumerate(arglist):
            mp_best_fit.plot(n, *args)

    # Best fitting process plots for GA.
    best_fit_algor = g.bf_params[1]
    if bf_flag and best_fit_algor == 'genet':

        min_max_p = prep_plots.param_ranges(ip_list)
        # Get special axis ticks for metallicity.
        xp_min, xp_max = min_max_p[0]
        # The maximum number of characters in the axis, '30', is HARD-CODED.
        # Add values to the end of this list.
        min_max_p.append(prep_plots.BestTick(xp_min, xp_max, 30))

        # Unpack.
        lkl_old, new_bs_indx, model_done = isoch_fit_params[1:4]
        l_min_max = prep_plots.likl_y_range(lkl_old)

        arglist = [
            [gs, l_min_max, lkl_old, model_done, new_bs_indx],
            [gs, 'age-metal', min_max_p, isoch_fit_params, isoch_fit_errors,
                model_done],
            [gs, 'dist-ext', min_max_p, isoch_fit_params, isoch_fit_errors,
                model_done],
            [gs, 'metal-dist', min_max_p, isoch_fit_params, isoch_fit_errors,
                model_done],
            [gs, 'mass-binar', min_max_p, isoch_fit_params, isoch_fit_errors,
                model_done],
            [gs, '$z$', l_min_max[1], lkl_old, min_max_p, isoch_fit_params,
                isoch_fit_errors, model_done],
            [gs, '$log(age)$', l_min_max[1], lkl_old, min_max_p,
                isoch_fit_params, isoch_fit_errors, model_done],
            [gs, '$E_{{(B-V)}}$', l_min_max[1], lkl_old, min_max_p,
                isoch_fit_params, isoch_fit_errors, model_done],
            [gs, '$(m-M)_o$', l_min_max[1], lkl_old, min_max_p,
                isoch_fit_params, isoch_fit_errors, model_done],
            [gs, '$M_{{\odot}}$', l_min_max[1], lkl_old, min_max_p,
                isoch_fit_params, isoch_fit_errors, model_done],
            [gs, '$b_{{frac}}$', l_min_max[1], lkl_old, min_max_p,
                isoch_fit_params, isoch_fit_errors, model_done]
            ]
        for n, args in enumerate(arglist, 1):
            mp_best_fit.plot(n, *args)

    # Ignore warning issued by colorbar plotted in photom diagram with
    # membership probabilities.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fig.tight_layout()

    # Plot colorbar down here so tight_layout won't move it around.
    try:
        if plot_colorbar is True:
            import matplotlib.transforms as mts
            # Position and dimensions relative to the axes.
            x0, y0, width, height = [0.67, 0.92, 0.2, 0.04]
            # Transform them to get the ABSOLUTE POSITION AND DIMENSIONS
            Bbox = mts.Bbox.from_bounds(x0, y0, width, height)
            l, b, w, h = mts.TransformedBbox(Bbox, trans).bounds
            # Create the axes and the colorbar.
            cbaxes = fig.add_axes([l, b, w, h])
            cbar = plt.colorbar(sca, cax=cbaxes, ticks=[v_min_mp, v_max_mp],
                orientation='horizontal')
            cbar.ax.tick_params(labelsize=9)
    except:
        #import traceback
        #print traceback.format_exc()
        print("  WARNING: error when plotting colorbar on cluster's "
        "photom diagram.")

    # Generate output file for each data file.
    pl_fmt, pl_dpi = g.pl_params[1:3]
    plt.savefig(join(output_subdir, str(clust_name) + '.' + pl_fmt), dpi=pl_dpi)

    # Close to release memory.
    plt.clf()
    plt.close()

    print 'Plots created.'
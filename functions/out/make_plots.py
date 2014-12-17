"""
@author: gabriel
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from os.path import join
import warnings
# Custom functions.
from .._in import get_in_params as g
from ..errors import error_round as err_r
import mp_structure
import mp_phot_analysis
import mp_decont_algor
import mp_best_fit
from mp_star_size import star_size


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


def line(x, slope, intercept):
    '''
    Linar function.
    '''
    y = slope * x + intercept
    return y


def reject_outliers(data, m=6.5):
    '''
    Reject outliers from array.
    http://stackoverflow.com/a/16562028/1391441
    '''
    d = np.abs(data - np.median(data))
    mdev = np.median(d)
    s = d / mdev if mdev else 0.
    return data[s < m]


def make_plots(output_subdir, clust_name, x_data, y_data,
    bin_width, center_params, rdp_params, field_dens, radius_params,
    cont_index, mag_data, col1_data, err_plot, err_flags, kp_params,
    cl_region, stars_out, stars_in_rjct, stars_out_rjct, integr_return, n_memb,
    flag_area_stronger, field_regions, flag_pval_test,
    pval_test_params, memb_prob_avrg_sort, lum_func, completeness, ip_list,
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
    isoch_fit_params, isoch_fit_errors, shift_isoch, synth_clst = bf_return

    # Plot all outputs
    # figsize(x1, y1), GridSpec(y2, x2) --> To have square plots: x1/x2 =
    # y1/y2 = 2.5
    fig = plt.figure(figsize=(30, 25))  # create the top-level container
    gs = gridspec.GridSpec(10, 12)       # create a GridSpec object

    # Get max and min values in x,y coordinates.
    x_min, x_max = min(x_data), max(x_data)
    y_min, y_max = min(y_data), max(y_data)

    # Define system of coordinates used.
    px_deg = g.gd_params[-1]
    coord_lst = ['px', 'x', 'y'] if px_deg == 'px' else ['deg', 'ra', 'dec']
    coord, x_name, y_name = coord_lst

    # If possible, define zoomed frame.
    x_zmin, x_zmax = max(x_min, (center_cl[0] - 1.5 * clust_rad)), \
    min(x_max, (center_cl[0] + 1.5 * clust_rad))
    y_zmin, y_zmax = max(y_min, (center_cl[1] - 1.5 * clust_rad)), \
    min(y_max, (center_cl[1] + 1.5 * clust_rad))
    # Prevent axis stretching.
    if (x_zmax - x_zmin) != (y_zmax - y_zmin):
        lst = [(x_zmax - x_zmin), (y_zmax - y_zmin)]
        val, idx = min((val, idx) for (idx, val) in enumerate(lst))
        if idx == 0:
            x_zmax = x_zmin + lst[1]
        else:
            y_zmax = y_zmin + lst[0]

    # Define names for photometric diagram axes.
    y_axis = 0
    y_ax, x_ax0, m_ord = g.axes_params[0:3]
    if m_ord == 21:
        x_ax = '(' + x_ax0 + '-' + y_ax + ')'
    elif m_ord == 12:
        x_ax = '(' + y_ax + '-' + x_ax0 + ')'

    # Unpack coordinates and photometric data.
    #x_data, y_data = id_coords[1:]
    phot_x = col1_data
    phot_y = mag_data

    # Define plot limits for *all* photometric diagrams.
    phot_x_s, phot_y_s = reject_outliers(phot_x), reject_outliers(phot_y)
    x_max_cmd, x_min_cmd = max(phot_x_s) + 0.5, min(phot_x_s) - 0.5
    y_min_cmd, y_max_cmd = max(phot_y_s) + 0.5, min(phot_y_s) - 0.5
    # If photometric axis y is a magnitude, make sure the brighter stars
    # are plotted.
    y_max_cmd = (min(phot_y) - 1.) if y_axis == 0 else y_max_cmd

    # Obtain proper sizes for plotted stars.
    st_sizes_arr = star_size(mag_data)

    #
    # Structure plots.
    arglist = [
        [gs, fig, x_name, y_name, cent_bin, clust_rad, bin_width,
    err_r, coord, hist_2d_g],
        [gs, x_min, x_max, y_min, y_max, x_name, y_name, coord, approx_cents,
            bin_width, st_dev_lst],
        [gs, fig, x_name, y_name, coord, x_min, x_max, y_min, y_max,
            center_cl, clust_rad, e_cent, kp_params, mag_data, x_data,
            y_data, st_sizes_arr],
        [gs, radii, ring_density, field_dens, coord, clust_name, kp_params,
            clust_rad, e_rad, poisson_error, bin_width],
        [gs, fig, x_zmin, x_zmax, y_zmin, y_zmax, x_name, y_name, coord,
            cont_index, kde_pl, x_data, y_data, st_sizes_arr, center_cl,
            clust_rad],
        [gs, fig, x_min, x_max, y_min, y_max, x_name, y_name, coord,
            center_cl, clust_rad, field_regions, cl_region, flag_area_stronger]
        ]
    for n, args in enumerate(arglist, 1):
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
            stars_out_rjct, field_regions],
        [gs, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd, x_ax, y_ax,
            stars_in_rjct, cl_region, n_memb],
        [gs, mag_data, y_ax, x_cl, y_cl, flag_area_stronger, x_fl, y_fl,
            completeness],
        [gs, integr_return, y_ax, x_ax0, flag_area_stronger, m_ord],
        [gs, flag_pval_test, flag_area_stronger, pval_test_params]
        ]
    for n, args in enumerate(arglist, 1):
        mp_phot_analysis.plot(n, *args)

    #
    # Decontamination algorithm plots.
    da_flag, bf_flag = g.da_params[0], g.bf_params[0]

    if g.da_params[0] != 'skip':

        arglist = [
            [gs, red_return, memb_prob_avrg_sort],
            [gs, fig, x_zmin, x_zmax, y_zmin, y_zmax, center_cl, clust_rad,
                field_dens, memb_prob_avrg_sort, stars_out]
            ]
        for n, args in enumerate(arglist, 1):
            mp_decont_algor.plot(n, *args)

    plot_colorbar = False
    if da_flag != 'skip' or bf_flag:
            # This function is called separately since we need to retrieve some
            # information from it to plot that #$%&! colorbar.
            try:
                plot_colorbar, v_min_mp, v_max_mp, sca, trans = \
                mp_decont_algor.pl_mps_phot_diag(gs, fig, x_min_cmd, x_max_cmd,
                    y_min_cmd, y_max_cmd, x_ax, y_ax, memb_prob_avrg_sort,
                    shift_isoch, col1_data, err_lst)
            except:
                print("  WARNING: error when plotting MPs on cluster's "
                "photom diagram.")

    #
    # Best fit plots.
    best_fit_algor, N_b = g.bf_params[1], g.bf_params[-1]

    if bf_flag:
        arglist = [
            [gs, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd, x_ax, y_ax,
                synth_clst, isoch_fit_params[0], isoch_fit_errors, shift_isoch]
            ]
        for n, args in enumerate(arglist, 1):
            mp_best_fit.plot(n, *args)

    # Best fitting process plots for GA.
    if bf_flag and best_fit_algor == 'genet':

        # Set errors to zero if bootstrap was not used, for plotting purposes.
        if N_b >= 2:
            p_errs = isoch_fit_errors
        else:
            p_errs = [0.] * len(isoch_fit_errors)

        # Set parameter ranges used by GA plots.
        min_max_p = []
        for param in ip_list[1]:
            # Set the delta for the parameter range. If only one value was
            # used, set a very small delta value.
            delta_p = (max(param) - min(param)) * 0.05 \
            if max(param) != min(param) else 0.001
            # Store parameter range.
            min_max_p.append([min(param) - delta_p, max(param) + delta_p])

        # Unpack.
        lkl_old, new_bs_indx, model_done = isoch_fit_params[1], \
        isoch_fit_params[2], isoch_fit_params[3]
        # Obtain here the minimum likelihood axis value.
        y_min_edge = max(0, min(lkl_old[0]) - 0.3 * min(lkl_old[0]))

        arglist = [
            [gs, 'age-metal', min_max_p, isoch_fit_params, p_errs, model_done],
            [gs, 'dist-ext', min_max_p, isoch_fit_params, p_errs, model_done],
            [gs, 'metal-dist', min_max_p, isoch_fit_params, p_errs, model_done],
            [gs, 'mass-binar', min_max_p, isoch_fit_params, p_errs, model_done],
            [gs, lkl_old, model_done, new_bs_indx],
            [gs, '$z$', lkl_old, min_max_p, isoch_fit_params, p_errs,
                model_done, y_min_edge],
            [gs, '$log(age)$', lkl_old, min_max_p, isoch_fit_params, p_errs,
                model_done, y_min_edge],
            [gs, '$E_{{(B-V)}}$', lkl_old, min_max_p, isoch_fit_params, p_errs,
                model_done, y_min_edge],
            [gs, '$(m-M)_o$', lkl_old, min_max_p, isoch_fit_params, p_errs,
                model_done, y_min_edge],
            [gs, '$M_{{\odot}}$', lkl_old, min_max_p, isoch_fit_params, p_errs,
                model_done, y_min_edge],
            [gs, '$b_{{frac}}$', lkl_old, min_max_p, isoch_fit_params, p_errs,
                model_done, y_min_edge]
            ]
        for n, args in enumerate(arglist, 2):
            mp_best_fit.plot(n, *args)

    # Ignore warning issued by colorbar plotted in photom diagram with
    # membership probabilities.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fig.tight_layout()

    # Plot colorbar down here so tight_layout won't move it around.
    if plot_colorbar is True:
        import matplotlib.transforms as mts
        # Position and dimensions relative to the axes.
        x0, y0, width, height = [0.67, 0.85, 0.2, 0.04]
        # Transform them to get the ABSOLUTE POSITION AND DIMENSIONS
        Bbox = mts.Bbox.from_bounds(x0, y0, width, height)
        l, b, w, h = mts.TransformedBbox(Bbox, trans).bounds
        # Create the axes and the colorbar.
        cbaxes = fig.add_axes([l, b, w, h])
        cbar = plt.colorbar(sca, cax=cbaxes, ticks=[v_min_mp, v_max_mp],
            orientation='horizontal')
        cbar.ax.tick_params(labelsize=9)

    # Generate output file for each data file.
    pl_fmt, pl_dpi = g.pl_params[1:3]
    plt.savefig(join(output_subdir, str(clust_name) + '.' + pl_fmt), dpi=pl_dpi)

    # Close to release memory.
    plt.clf()
    plt.close()
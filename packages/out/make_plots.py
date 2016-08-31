
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from os.path import join
import warnings
from .._version import __version__
import mp_phot_analysis
import mp_decont_algor
import mp_best_fit
import prep_plots


#############################################################
# # Timer function: http://stackoverflow.com/a/21860100/1391441
# from contextlib import contextmanager
# import time
# @contextmanager
# def timeblock(label):
#     start = time.clock()
#     try:
#         yield
#     finally:
#         end = time.clock()
#         print ('{} elapsed: {}'.format(label, end - start))
#############################################################


def main(
        npd, cld, pd, bin_width, cent_bin, clust_cent, e_cent, approx_cents,
        st_dev_lst, hist_2d_g, kde_plot, radii, rdp_points, poisson_error,
        field_dens, clust_rad, e_rad, cont_index, err_plot, err_flags,
        core_rad, e_core, tidal_rad, e_tidal, K_conct_par, K_cent_dens,
        flag_2pk_conver, flag_3pk_conver, cl_region, stars_out, cl_region_rjct,
        stars_out_rjct, integr_return, n_memb, n_memb_da, flag_no_fl_regs,
        field_regions, flag_pval_test, pval_test_params, lum_func,
        completeness, memb_prob_avrg_sort, flag_decont_skip, cl_reg_fit,
        cl_reg_no_fit, cl_reg_clean_plot, err_lst, isoch_fit_params,
        isoch_fit_errors, shift_isoch, synth_clst, syn_b_edges, **kwargs):
    '''
    Make all plots.
    '''

    # flag_make_plot = pd['pl_params'][0]
    if pd['pl_params'][0]:

        # Unpack params.
        x, y, mags, cols = cld['x'], cld['y'], cld['mags'], cld['cols']
        # Luminosity functions.
        x_cl, y_cl, x_fl, y_fl = lum_func

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
        x_min, x_max, y_min, y_max = prep_plots.frame_max_min(x, y)
        coord, x_name, y_name = prep_plots.coord_syst(pd['coords'])
        x_zmin, x_zmax, y_zmin, y_zmax = prep_plots.frame_zoomed(
            x_min, x_max, y_min, y_max, clust_cent, clust_rad)
        x_ax, y_ax, x_ax0, y_axis = prep_plots.ax_names(pd['axes_params'])
        phot_x, phot_y = prep_plots.ax_data(mags, cols)
        x_max_cmd, x_min_cmd, y_min_cmd, y_max_cmd = prep_plots.diag_limits(
            y_axis, phot_x, phot_y)
        stars_f_rjct, stars_f_acpt = prep_plots.field_region_stars(
            stars_out_rjct, field_regions)
        f_sz_pt = prep_plots.phot_diag_st_size(len(stars_f_acpt[0]))
        cl_sz_pt = prep_plots.phot_diag_st_size(len(cl_region))
        v_min_mp, v_max_mp, plot_colorbar, chart_fit_inv, chart_no_fit_inv, \
            out_clust_rad, diag_fit_inv, diag_no_fit_inv, err_bar = \
            prep_plots.da_plots(
                clust_cent, clust_rad, stars_out, x_zmin, x_zmax, y_zmin,
                y_zmax, x_max_cmd, cols, err_lst, cl_reg_fit, cl_reg_no_fit)

        #
        # Photometric analysis plots.
        arglist = [
            # pl_phot_err: Photometric error rejection.
            [gs, fig, pd['er_params'], 'up', x_ax, y_ax, mags, err_plot,
             err_flags, cl_region, cl_region_rjct, stars_out, stars_out_rjct],
            [gs, fig, pd['er_params'], 'low', x_ax, y_ax, mags, err_plot,
             err_flags, cl_region, cl_region_rjct, stars_out, stars_out_rjct],
            # pl_fl_diag: Field stars CMD/CCD diagram.
            [gs, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd, x_ax, y_ax,
                stars_f_rjct, stars_f_acpt, f_sz_pt],
            # pl_cl_diag: Cluster's stars diagram (stars inside cluster's rad)
            [gs, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd, x_ax, y_ax,
                cl_region_rjct, cl_region, n_memb, cl_sz_pt],
            # pl_lum_func: LF of stars in cluster region and outside.
            [gs, mags, y_ax, x_cl, y_cl, flag_no_fl_regs, x_fl, y_fl,
                completeness],
            # pl_integ_mag: Integrated magnitudes.
            [gs, pd['axes_params'], integr_return, y_ax, x_ax0,
             flag_no_fl_regs],
            # pl_p_vals: Distribution of KDE p_values.
            [gs, flag_pval_test, pval_test_params]
        ]
        for n, args in enumerate(arglist):
            mp_phot_analysis.plot(n, *args)

        #
        # Decontamination algorithm plots.
        mode_red_memb, local_bin = pd['rm_params'][0], pd['rm_params'][1]
        lkl_method, bin_method, N_b = pd['bf_params'][2], pd['bf_params'][3],\
            pd['bf_params'][4]
        min_prob, bin_edges = cl_reg_clean_plot
        bf_flag = pd['bf_flag']

        # If the DA and the best fit functions were skipped and the reduced
        # membership mode is any mode but 'local', do not plot.
        if flag_decont_skip and bf_flag is False and mode_red_memb != 'local':
            pass
        else:
            arglist = [
                # pl_mp_histo
                [gs, n_memb_da, memb_prob_avrg_sort, flag_decont_skip,
                 cl_reg_fit, min_prob, mode_red_memb, local_bin],
                # pl_chart_mps
                [gs, fig, x_name, y_name, coord, x_zmin, x_zmax, y_zmin,
                 y_zmax, clust_cent, clust_rad, field_dens, flag_decont_skip,
                    v_min_mp, v_max_mp, chart_fit_inv, chart_no_fit_inv,
                    out_clust_rad, mode_red_memb, local_bin]
            ]
            for n, args in enumerate(arglist):
                mp_decont_algor.plot(n, *args)

            # This function is called separately since we need to retrieve some
            # information from it to plot that #$%&! colorbar.
            try:
                # pl_mps_phot_diag
                sca, trans = mp_decont_algor.pl_mps_phot_diag(
                    gs, fig, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd, x_ax,
                    y_ax, v_min_mp, v_max_mp, diag_fit_inv, diag_no_fit_inv,
                    shift_isoch, err_bar, mode_red_memb, bin_edges, bf_flag)
            except:
                # import traceback
                # print traceback.format_exc()
                print("  WARNING: error when plotting MPs on cluster's "
                      "photometric diagram.")

        #
        # Best fit plots.
        if bf_flag:
            arglist = [
                # pl_bf_synth_cl: Best fit synthetic cluster obtained.
                [gs, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd, x_ax, y_ax,
                 synth_clst, syn_b_edges, isoch_fit_params[0],
                 isoch_fit_errors, shift_isoch, lkl_method, bin_method,
                 pd['cmd_evol_tracks'], pd['evol_track']]
            ]
            for n, args in enumerate(arglist):
                mp_best_fit.plot(n, *args)

        # Best fitting process plots for GA.
        best_fit_algor = pd['bf_params'][1]
        if bf_flag and best_fit_algor == 'genet':

            min_max_p = prep_plots.param_ranges(pd['ip_list'])
            # Get special axis ticks for metallicity.
            xp_min, xp_max = min_max_p[0]
            # The max number of characters in the axis, '30', is HARD-CODED.
            # Add values to the end of this list.
            min_max_p.append(prep_plots.BestTick(xp_min, xp_max, 30))

            # Unpack.
            lkl_old, new_bs_indx, model_done = isoch_fit_params[1:4]
            l_min_max = prep_plots.likl_y_range(lkl_old)

            arglist = [
                # pl_ga_lkl: Likelihood evolution for the GA.
                [gs, l_min_max, lkl_old, model_done, new_bs_indx,
                 pd['ga_params'], N_b],
                # pl_2_param_dens: Param vs param solutions scatter map.
                [gs, 'age-metal', min_max_p, isoch_fit_params,
                 isoch_fit_errors, model_done],
                [gs, 'dist-ext', min_max_p, isoch_fit_params, isoch_fit_errors,
                    model_done],
                [gs, 'metal-dist', min_max_p, isoch_fit_params,
                 isoch_fit_errors, model_done],
                [gs, 'mass-binar', min_max_p, isoch_fit_params,
                 isoch_fit_errors, model_done],
                # pl_lkl_scatt: Parameter likelihood density plot.
                [gs, '$z$', min_max_p, isoch_fit_params, isoch_fit_errors,
                    model_done],
                [gs, '$log(age)$', min_max_p, isoch_fit_params,
                 isoch_fit_errors, model_done],
                [gs, '$E_{{(B-V)}}$', min_max_p, isoch_fit_params,
                    isoch_fit_errors, model_done],
                [gs, '$(m-M)_o$', min_max_p, isoch_fit_params,
                 isoch_fit_errors, model_done],
                [gs, '$M\,(M_{{\odot}})$', min_max_p, isoch_fit_params,
                    isoch_fit_errors, model_done],
                [gs, '$b_{{frac}}$', min_max_p, isoch_fit_params,
                 isoch_fit_errors, model_done]
            ]
            for n, args in enumerate(arglist, 1):
                mp_best_fit.plot(n, *args)

        # Ignore warning issued by colorbar plotted in photometric diagram with
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
                cbar = plt.colorbar(
                    sca, cax=cbaxes, ticks=[v_min_mp, v_max_mp],
                    orientation='horizontal')
                cbar.ax.tick_params(labelsize=9)
        except:
            # import traceback
            # print traceback.format_exc()
            print("  WARNING: error when plotting colorbar on cluster's "
                  "photometric diagram.")

        # Generate output file for each data file.
        pl_fmt, pl_dpi = pd['pl_params'][1:3]
        plt.savefig(
            join(npd['output_subdir'], str(npd['clust_name']) + '.' + pl_fmt),
            dpi=pl_dpi)

        # Close to release memory.
        plt.clf()
        plt.close()

        print('Plots created.')

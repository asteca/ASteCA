
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from os.path import join
from .._version import __version__
import mp_best_fit
import prep_plots
from ..synth_clust import synth_cl_plot


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
        isoch_fit_errors, syn_b_edges, **kwargs):
    '''
    Make all plots.
    '''

    # flag_make_plot = pd['pl_params'][0]
    if pd['pl_params'][0]:

        # Unpack params.
        x, y, mags, cols = cld['x'], cld['y'], cld['mags'], cld['cols']

        # Plot all outputs
        # figsize(x1, y1), GridSpec(y2, x2) --> To have square plots: x1/x2 =
        # y1/y2 = 2.5
        plt.figure(figsize=(30, 25))  # create the top-level container
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
        v_min_mp, v_max_mp, plot_colorbar, chart_fit_inv, chart_no_fit_inv, \
            out_clust_rad, diag_fit_inv, diag_no_fit_inv, err_bar = \
            prep_plots.da_plots(
                clust_cent, clust_rad, stars_out, x_zmin, x_zmax, y_zmin,
                y_zmax, x_max_cmd, cols, err_lst, cl_reg_fit, cl_reg_no_fit)

        # Obtain best fit synthetic cluster, and its isochrone.
        shift_isoch, synth_clst = synth_cl_plot.main(
            ip_list, isoch_fit_params, err_lst, completeness, st_dist_mass,
            e_max, bin_mr, cmd_sel)
        if not synth_clst.any():
            print("  WARNING: best fit synthetic cluster found is empty.")

        #
        # Best fit plots.
        if pd['bf_flag']:
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
        if pd['bf_flag'] and best_fit_algor == 'genet':

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

        # Generate output file for each data file.
        pl_fmt, pl_dpi = pd['pl_params'][1:3]
        plt.savefig(
            join(npd['output_subdir'], str(npd['clust_name']) + '.' + pl_fmt),
            dpi=pl_dpi)

        # Close to release memory.
        plt.clf()
        plt.close()

        print('Plots created.')

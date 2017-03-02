
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from os.path import join
import warnings
import add_version_plot
import mp_best_fit
import prep_plots


def plot_observed_cluster(
    cld, pd, fig, gs, cl_reg_fit, cl_reg_no_fit, err_lst, v_min_mp, v_max_mp,
        plot_colorbar, diag_fit_inv, lkl_method, hess_data, shift_isoch):
    """
    This function is called separately since we need to retrieve some
    information from it to plot that #$%&! colorbar.
    """
    x_ax, y_ax, y_axis = prep_plots.ax_names(pd['filters'], pd['colors'])
    x_max_cmd, x_min_cmd, y_min_cmd, y_max_cmd = prep_plots.diag_limits(
        y_axis, cld['cols'], cld['mags'])
    err_bar = prep_plots.error_bars(
        cl_reg_fit + cl_reg_no_fit, x_min_cmd, err_lst)

    try:
        # pl_mps_phot_diag
        sca, trans = mp_best_fit.pl_mps_phot_diag(
            gs, fig, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd,
            x_ax, y_ax, v_min_mp, v_max_mp, diag_fit_inv,
            err_bar, lkl_method, hess_data, shift_isoch)
    except:
        import traceback
        print traceback.format_exc()
        print("  WARNING: error when plotting MPs on cluster's "
              "photometric diagram.")

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
            x0, y0, width, height = [0.74, 0.93, 0.2, 0.04]
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


def main(
        npd, cld, pd, synth_clst, shift_isoch, isoch_fit_params,
        isoch_fit_errors, cl_reg_fit, cl_reg_no_fit, err_lst,
        **kwargs):
    '''
    Make D block plots.
    '''

    # flag_make_plot = pd['pl_params'][0]
    if pd['pl_params'][0]:
        fig = plt.figure(figsize=(30, 25))
        gs = gridspec.GridSpec(10, 12)
        add_version_plot.main()

        best_fit_algor, lkl_method, bin_method, N_b = pd['bf_params']
        x_ax, y_ax, y_axis = prep_plots.ax_names(pd['filters'], pd['colors'])
        x_max_cmd, x_min_cmd, y_min_cmd, y_max_cmd = prep_plots.diag_limits(
            y_axis, cld['cols'], cld['mags'])
        hess_data = prep_plots.get_hess(
            lkl_method, bin_method, cl_reg_fit, synth_clst)

        # Best fit plots.
        if pd['bf_flag']:
            arglist = [
                # hess_diag_pl: Hess diagram 'observed - synthetic'
                [gs, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd, x_ax, y_ax,
                 lkl_method, hess_data],
                # pl_bf_synth_cl: Best fit synthetic cluster obtained.
                [gs, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd, x_ax, y_ax,
                 synth_clst, hess_data, pd['IMF_name'], pd['R_V'],
                 isoch_fit_params[0], isoch_fit_errors, shift_isoch,
                 lkl_method, bin_method, pd['cmd_evol_tracks'],
                 pd['evol_track']]
            ]
            for n, args in enumerate(arglist):
                mp_best_fit.plot(n, *args)

            # Best fitting process plots for GA.
            if best_fit_algor == 'genet':

                min_max_p = prep_plots.param_ranges(pd['fundam_params'])
                # Get special axis ticks for metallicity.
                xp_min, xp_max = min_max_p[0]
                # The max number of characters in the axis '30', is HARD-CODED.
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
                    [gs, 'dist-ext', min_max_p, isoch_fit_params,
                        isoch_fit_errors, model_done],
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
                for n, args in enumerate(arglist, 2):
                    mp_best_fit.plot(n, *args)

            # tight_layout is called here
            v_min_mp, v_max_mp = prep_plots.da_colorbar_range(
                cl_reg_fit, cl_reg_no_fit)
            plot_colorbar, diag_fit_inv, diag_no_fit_inv =\
                prep_plots.da_phot_diag(
                    cl_reg_fit, cl_reg_no_fit, v_min_mp, v_max_mp)
            plot_observed_cluster(
                cld, pd, fig, gs, cl_reg_fit, cl_reg_no_fit, err_lst, v_min_mp,
                v_max_mp, plot_colorbar, diag_fit_inv, lkl_method, hess_data,
                shift_isoch)

        # Generate output file.
        pl_fmt, pl_dpi = pd['pl_params'][1:3]
        plt.savefig(
            join(npd['output_subdir'], str(npd['clust_name']) +
                 '_D.' + pl_fmt), dpi=pl_dpi, bbox_inches='tight')
        # Close to release memory.
        plt.clf()
        plt.close()

        print("<<Plots from 'D' block created>>")


import itertools
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from os.path import join
import warnings
import add_version_plot
import mp_best_fit2
import prep_plots


def plot_observed_cluster(
    cld, pd, fig, gs, cl_max_mag, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd,
    err_lst, v_min_mp, v_max_mp, plot_colorbar, diag_fit_inv, lkl_method,
        hess_data, shift_isoch):
    """
    This function is called separately since we need to retrieve some
    information from it to plot that #$%&! colorbar.
    """
    x_ax, y_ax, y_axis = prep_plots.ax_names(pd['filters'], pd['colors'])
    err_bar = prep_plots.error_bars(cl_max_mag, x_min_cmd, err_lst)

    try:
        # pl_mps_phot_diag
        sca, trans = mp_best_fit2.pl_mps_phot_diag(
            gs, fig, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd,
            x_ax, y_ax, v_min_mp, v_max_mp, diag_fit_inv,
            err_bar, lkl_method, hess_data, shift_isoch)
    except Exception:
        # import traceback
        # print traceback.format_exc()
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
    except Exception:
        # import traceback
        # print traceback.format_exc()
        print("  WARNING: error when plotting colorbar on cluster's "
              "photometric diagram.")


def main(npd, cld, pd, synth_clst, shift_isoch, fit_params_r, fit_errors_r,
         cl_max_mag, err_lst, **kwargs):
    '''
    Make D2 block plots.
    '''
    if pd['flag_make_plot'] and pd['bf_flag']:
        fig = plt.figure(figsize=(30, 25))
        gs = gridspec.GridSpec(10, 12)
        add_version_plot.main()

        x_ax, y_ax, y_axis = prep_plots.ax_names(pd['filters'], pd['colors'])
        # TODO using first magnitude and color defined
        firts_col = list(itertools.chain.from_iterable(zip(*cl_max_mag)[5]))
        first_mag = list(itertools.chain.from_iterable(zip(*cl_max_mag)[3]))
        x_max_cmd, x_min_cmd, y_min_cmd, y_max_cmd =\
            prep_plots.diag_limits(y_axis, firts_col, first_mag)
        hess_data = prep_plots.get_hess(
            pd['lkl_method'], pd['lkl_binning'], cl_max_mag, synth_clst)

        arglist = [
            # hess_diag_pl: Hess diagram 'observed - synthetic'
            [gs, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd, x_ax, y_ax,
             pd['lkl_method'], hess_data],
            # pl_bf_synth_cl: Best fit synthetic cluster obtained.
            [gs, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd, x_ax, y_ax,
             synth_clst, hess_data, pd['IMF_name'], pd['R_V'],
             fit_params_r, fit_errors_r, shift_isoch,
             pd['lkl_method'], pd['lkl_binning'], pd['cmd_evol_tracks'],
             pd['evol_track']]
        ]
        for n, args in enumerate(arglist):
            mp_best_fit2.plot(n, *args)

        # tight_layout is called here
        v_min_mp, v_max_mp = prep_plots.da_colorbar_range(cl_max_mag, [])
        plot_colorbar, diag_fit_inv, dummy = prep_plots.da_phot_diag(
            cl_max_mag, [], v_min_mp, v_max_mp)
        # Main photometric diagram of observed cluster.
        plot_observed_cluster(
            cld, pd, fig, gs, cl_max_mag, x_min_cmd, x_max_cmd, y_min_cmd,
            y_max_cmd, err_lst, v_min_mp, v_max_mp, plot_colorbar,
            diag_fit_inv, pd['lkl_method'], hess_data, shift_isoch)

        # Generate output file.
        plt.savefig(
            join(npd['output_subdir'], str(npd['clust_name']) +
                 '_D2.' + pd['plot_frmt']), dpi=pd['plot_dpi'],
            bbox_inches='tight')
        # Close to release memory.
        plt.clf()
        plt.close()

        print("<<Plots for D2 block created>>")


import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from os.path import join
import warnings
import add_version_plot
import mp_decont_algor
import prep_plots


def main(
        npd, cld, pd, cl_reg_clean_plot, flag_decont_skip, n_memb_da,
        memb_prob_avrg_sort, cl_reg_fit, cl_reg_no_fit, clust_cent,
        clust_rad, field_dens, stars_out, err_lst, **kwargs):
    '''
    Make C block plots.
    '''
    # flag_make_plot = pd['pl_params'][0]
    if pd['pl_params'][0]:
        # Plot all outputs
        # figsize(x1, y1), GridSpec(y2, x2) --> To have square plots: x1/x2 =
        # y1/y2 = 2.5
        fig = plt.figure(figsize=(30, 25))  # create the top-level container
        gs = gridspec.GridSpec(10, 12)      # create a GridSpec object
        add_version_plot.main()

        # Obtain plotting parameters and data.
        x_min, x_max, y_min, y_max = prep_plots.frame_max_min(
            cld['x'], cld['y'])
        x_zmin, x_zmax, y_zmin, y_zmax = prep_plots.frame_zoomed(
            x_min, x_max, y_min, y_max, clust_cent, clust_rad)
        coord, x_name, y_name = prep_plots.coord_syst(pd['coords'])
        x_ax, y_ax, y_axis = prep_plots.ax_names(pd['filters'], pd['colors'])
        x_max_cmd, x_min_cmd, y_min_cmd, y_max_cmd = prep_plots.diag_limits(
            y_axis, cld['cols'], cld['mags'])
        v_min_mp, v_max_mp, plot_colorbar, chart_fit_inv, chart_no_fit_inv, \
            out_clust_rad, diag_fit_inv, diag_no_fit_inv = prep_plots.da_plots(
                clust_cent, clust_rad, stars_out, x_zmin, x_zmax, y_zmin,
                y_zmax, cl_reg_fit, cl_reg_no_fit)
        err_bar = prep_plots.error_bars(
            cl_reg_fit + cl_reg_no_fit, x_min_cmd, err_lst)

        # Decontamination algorithm plots.
        mode_red_memb, local_bin = pd['rm_params'][0], pd['rm_params'][1]
        min_prob, bin_edges = cl_reg_clean_plot

        arglist = [
            # pl_mp_histo
            [gs, n_memb_da, memb_prob_avrg_sort, flag_decont_skip,
             cl_reg_fit, min_prob, mode_red_memb, local_bin],
            # pl_chart_mps
            [gs, fig, x_name, y_name, coord, x_zmin, x_zmax, y_zmin,
             y_zmax, clust_cent, clust_rad, field_dens, flag_decont_skip,
             v_min_mp, v_max_mp, chart_fit_inv, chart_no_fit_inv,
             out_clust_rad, mode_red_memb, local_bin]]
        for n, args in enumerate(arglist):
            mp_decont_algor.plot(n, *args)

        # This function is called separately since we need to retrieve some
        # information from it to plot that #$%&! colorbar.
        try:
            # pl_mps_phot_diag
            sca, trans = mp_decont_algor.pl_mps_phot_diag(
                gs, fig, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd, x_ax,
                y_ax, v_min_mp, v_max_mp, diag_fit_inv, diag_no_fit_inv,
                err_bar, mode_red_memb, bin_edges)
        except:
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
        except:
            # import traceback
            # print traceback.format_exc()
            print("  WARNING: error when plotting colorbar on cluster's "
                  "photometric diagram.")

        # Generate output file for each data file.
        pl_fmt, pl_dpi = pd['pl_params'][1:3]
        plt.savefig(
            join(npd['output_subdir'], str(npd['clust_name']) +
                 '_C.' + pl_fmt), dpi=pl_dpi, bbox_inches='tight')

        # Close to release memory.
        plt.clf()
        plt.close()

        print("<<Plots from 'C' block created>>")

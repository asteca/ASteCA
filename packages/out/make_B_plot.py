
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from os.path import join
import warnings
import mp_phot_analysis
import add_version_plot
import prep_plots


def main(
        npd, cld, pd, err_plot, err_flags, cl_region, cl_region_rjct,
        stars_out, stars_out_rjct, err_lst, field_regions, n_memb,
        flag_no_fl_regs, lum_func, completeness, cl_reg_imag, fl_reg_imag,
        integ_mag, flag_pval_test, pval_test_params, **kwargs):
    '''
    Make B block plots.
    '''
    # flag_make_plot = pd['pl_params'][0]
    if pd['pl_params'][0]:
        # figsize(x1, y1), GridSpec(y2, x2) --> To have square plots: x1/x2 =
        # y1/y2 = 2.5
        fig = plt.figure(figsize=(30, 25))  # create the top-level container
        gs = gridspec.GridSpec(10, 12)      # create a GridSpec object
        add_version_plot.main()

        # # Obtain plotting parameters and data.
        x_ax, y_ax, y_axis = prep_plots.ax_names(pd['filters'], pd['colors'])
        x_max_cmd, x_min_cmd, y_min_cmd, y_max_cmd = prep_plots.diag_limits(
            y_axis, cld['cols'], cld['mags'])
        stars_f_rjct, stars_f_acpt = prep_plots.field_region_stars(
            stars_out_rjct, field_regions)
        f_sz_pt = prep_plots.phot_diag_st_size(len(stars_f_acpt[0]))
        cl_sz_pt = prep_plots.phot_diag_st_size(len(cl_region))
        err_bar_fl = prep_plots.error_bars(
            stars_out, cld['cols'], x_min_cmd, err_lst)
        err_bar_cl = prep_plots.error_bars(
            cl_region, cld['cols'], x_min_cmd, err_lst)

        # Photometric analysis plots.
        arglist = [
            # pl_phot_err: Photometric error rejection.
            [gs, fig, pd['er_params'], 'up', x_ax, y_ax, cld['mags'],
             err_plot, err_flags, cl_region, cl_region_rjct, stars_out,
             stars_out_rjct],
            [gs, fig, pd['er_params'], 'low', x_ax, y_ax, cld['mags'],
             err_plot, err_flags, cl_region, cl_region_rjct, stars_out,
             stars_out_rjct],
            # pl_fl_diag: Field stars CMD/CCD diagram.
            [gs, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd, x_ax, y_ax,
                stars_f_rjct, stars_f_acpt, f_sz_pt, err_bar_fl],
            # pl_cl_diag: Cluster's stars diagram (stars inside cluster's rad)
            [gs, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd, x_ax, y_ax,
                cl_region_rjct, cl_region, n_memb, cl_sz_pt, err_bar_cl],
            # pl_lum_func: LF of stars in cluster region and outside.
            [gs, cld['mags'], y_ax, flag_no_fl_regs, lum_func, completeness],
            # pl_integ_mag: Integrated magnitudes.
            [gs, cl_reg_imag, fl_reg_imag, integ_mag, y_ax, flag_no_fl_regs],
            # pl_p_vals: Distribution of KDE p_values.
            [gs, flag_pval_test, pval_test_params]
        ]
        for n, args in enumerate(arglist):
            mp_phot_analysis.plot(n, *args)

        # Ignore warning issued by colorbar plotted in photometric diagram with
        # membership probabilities.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fig.tight_layout()

        # Generate output file for each data file.
        pl_fmt, pl_dpi = pd['pl_params'][1:3]
        plt.savefig(
            join(npd['output_subdir'], str(npd['clust_name']) +
                 '_B.' + pl_fmt), dpi=pl_dpi, bbox_inches='tight')

        # Close to release memory.
        plt.clf()
        plt.close()

        print("<<Plots from 'B' block created>>")

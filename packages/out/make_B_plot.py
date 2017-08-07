
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from os.path import join
import mp_phot_analysis
import add_version_plot
import prep_plots


def main(
        npd, cld, pd, err_max, cl_region, cl_region_rjct,
        stars_out, stars_out_rjct, err_lst, field_regions, n_memb,
        flag_no_fl_regs, lum_func, completeness, cl_reg_imag, fl_reg_imag,
        integ_mag, flag_pval_test, pval_test_params, **kwargs):
    '''
    Make B block plots.
    '''
    if 'B' in pd['flag_make_plot']:
        fig = plt.figure(figsize=(30, 25))
        gs = gridspec.GridSpec(10, 12)
        add_version_plot.main()

        # # Obtain plotting parameters and data.
        x_ax, y_ax, y_axis = prep_plots.ax_names(pd['filters'], pd['colors'])
        # TODO using first magnitude and color defined
        x_max_cmd, x_min_cmd, y_min_cmd, y_max_cmd = prep_plots.diag_limits(
            y_axis, cld['cols'][0], cld['mags'][0])
        stars_f_rjct, stars_f_acpt = prep_plots.field_region_stars(
            stars_out_rjct, field_regions)
        f_sz_pt = prep_plots.phot_diag_st_size(len(stars_f_acpt[0]))
        cl_sz_pt = prep_plots.phot_diag_st_size(len(cl_region))
        err_bar_fl = prep_plots.error_bars(stars_out, x_min_cmd, err_lst)
        err_bar_cl = prep_plots.error_bars(cl_region, x_min_cmd, err_lst)

        # Photometric analysis plots.
        arglist = [
            # pl_phot_err: Photometric error rejection.
            [gs, fig, 'up', x_ax, y_ax, cld['mags'], err_max, cl_region,
             cl_region_rjct, stars_out, stars_out_rjct],
            [gs, fig, 'low', x_ax, y_ax, cld['mags'], err_max, cl_region,
             cl_region_rjct, stars_out, stars_out_rjct],
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

        # Generate output file.
        fig.tight_layout()
        plt.savefig(
            join(npd['output_subdir'], str(npd['clust_name']) +
                 '_B.' + pd['plot_frmt']), dpi=pd['plot_dpi'],
            bbox_inches='tight')
        # Close to release memory.
        plt.clf()
        plt.close()

        print("<<Plots for B block created>>")
    else:
        print("<<Skip B block plot>>")

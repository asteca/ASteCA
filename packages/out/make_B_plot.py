
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from os.path import join
from . import mp_data_analysis
from . import add_version_plot
from . import prep_plots


def main(
        npd, cld_c, pd, em_float, err_lst, cl_region_c, cl_region_rjct_c,
        stars_out_c, stars_out_rjct_c, field_regions_c, flag_no_fl_regs_c,
        field_regions_rjct_c, n_memb, lum_func, phot_analy_compl,
        phot_data_compl, err_rm_data, completeness, flag_ad_test, ad_cl,
        ad_fr, ad_cl_fr_p, ad_cl_fr_pk, stars_f_acpt, stars_f_rjct,
        col_0_comb, col_1_comb, mag_0_comb, **kwargs):
    '''
    Make B block plots.
    '''
    if 'B' in pd['flag_make_plot']:
        fig = plt.figure(figsize=(30, 25))
        gs = gridspec.GridSpec(10, 12)
        y_ver = .999
        # If no kinematic data was used.
        if np.array([_ == 'n' for _ in (
                pd['id_kinem'][0], pd['id_kinem'][2],
                pd['id_kinem'][6])]).all():
            y_ver = .905
        add_version_plot.main(y_fix=y_ver)

        # Obtain plotting parameters and data.
        x_ax0, y_ax = prep_plots.ax_names(
            pd['colors'][0], pd['filters'][0], 'mag')
        x_max_cmd0, x_min_cmd0, y_min_cmd0, y_max_cmd0 =\
            prep_plots.diag_limits('mag', col_0_comb, mag_0_comb)
        err_bar_fl0 = prep_plots.error_bars(stars_out_c, x_min_cmd0, err_lst)
        err_bar_cl0 = prep_plots.error_bars(cl_region_c, x_min_cmd0, err_lst)
        # For the error curve 'x_min_cmd' is not used
        err_bar_all = prep_plots.error_bars(
            cld_c['mags'][0], np.nan, err_lst, 'all')

        x_ax1, x_min_cmd1, x_max_cmd1, y_min_cmd1, y_max_cmd1,\
            err_bar_fl1, err_bar_cl1 = '', .0, .0, .0, .0, [], []
        if len(pd['colors']) > 1:
            x_ax1, dummy = prep_plots.ax_names(
                pd['colors'][1], pd['filters'][0], 'mag')
            x_max_cmd1, x_min_cmd1, y_min_cmd1, y_max_cmd1 =\
                prep_plots.diag_limits('mag', col_1_comb, mag_0_comb)
            err_bar_fl1 = prep_plots.error_bars(
                stars_out_c, x_min_cmd1, err_lst)
            err_bar_cl1 = prep_plots.error_bars(
                cl_region_c, x_min_cmd1, err_lst)

        f_sz_pt = prep_plots.phot_diag_st_size(len(stars_f_acpt[0]))
        cl_sz_pt = prep_plots.phot_diag_st_size(len(cl_region_c))
        x_min, x_max, y_min, y_max = prep_plots.frame_max_min(
            cld_c['x'], cld_c['y'])
        asp_ratio = prep_plots.aspect_ratio(x_min, x_max, y_min, y_max)
        coord, x_name, y_name = prep_plots.coord_syst(pd['coords'])

        b = 2 if len(pd['colors']) > 1 else 0

        # Photometric analysis plots.
        arglist = [
            # pl_phot_err: Photometric error rejection.
            [gs, pd['colors'], pd['filters'], pd['id_kinem'], cld_c['mags'],
             em_float, cl_region_c, cl_region_rjct_c, stars_out_c,
             stars_out_rjct_c, err_bar_all],
            # pl_cl_fl_regions
            [gs, x_name, y_name, coord, x_min, x_max, y_min, y_max, asp_ratio,
             field_regions_rjct_c, cl_region_rjct_c, flag_no_fl_regs_c],
            # pl_fl_diag: Field stars CMD/CCD diagram.
            [gs, x_ax0, y_ax, x_min_cmd0, x_max_cmd0, y_min_cmd0, y_max_cmd0,
             x_ax1, x_min_cmd1, x_max_cmd1, y_min_cmd1, y_max_cmd1,
             field_regions_c, stars_f_rjct, stars_f_acpt, f_sz_pt,
             err_bar_fl0, err_bar_fl1],
            # pl_cl_diag: Cluster's stars diagram (stars inside cluster's rad)
            [gs, x_ax0, y_ax, x_min_cmd0, x_max_cmd0, y_min_cmd0, y_max_cmd0,
             x_ax1, x_min_cmd1, x_max_cmd1, y_min_cmd1, y_max_cmd1,
             err_bar_cl0, err_bar_cl1, cl_region_rjct_c, cl_region_c, n_memb,
             cl_sz_pt],
            # pl_hess_cmd
            [gs, x_ax0, x_ax1, y_ax, x_max_cmd0, x_min_cmd0, y_min_cmd0,
             y_max_cmd0, x_max_cmd1, x_min_cmd1, y_min_cmd1, y_max_cmd1,
             stars_f_acpt, cl_region_c],
            # pl_lum_func: LF of stars in cluster region and outside.
            [gs, y_ax, flag_no_fl_regs_c, lum_func],
            # pl_data_rm_perc
            [gs, y_ax, phot_analy_compl, phot_data_compl, err_rm_data,
             completeness],
            # pl_ad_test
            [gs, b, flag_ad_test, ad_cl, ad_fr, pd['id_kinem'],
             pd['ad_k_comb']],
            # pl_ad_pvals_phot
            [gs, b, flag_ad_test, ad_cl_fr_p],
            # pl_ad_pvals_pk
            [gs, b, flag_ad_test, ad_cl_fr_pk, pd['id_kinem'], pd['ad_k_comb']]
        ]
        for n, args in enumerate(arglist):
            mp_data_analysis.plot(n, *args)

        # Generate output file.
        fig.tight_layout()
        plt.savefig(
            join(npd['output_subdir'], str(npd['clust_name']) +
                 '_B.' + pd['plot_frmt']), dpi=pd['plot_dpi'],
            bbox_inches='tight')
        # Close to release memory.
        plt.clf()
        plt.close("all")

        print("<<Plots for B block created>>")
    else:
        print("<<Skip B block plot>>")

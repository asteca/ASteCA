
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from os.path import join
from . import mp_data_analysis
from . import add_version_plot
from . import prep_plots
from . prep_plots import figsize_x, figsize_y, grid_x, grid_y


def main(
    npd, cld_i, pd, err_lst, cl_region_c, cl_region_rjct_c, stars_out_c,
    stars_out_rjct_c, field_regions_c, flag_no_fl_regs_c, field_regions_rjct_c,
    n_memb, lum_func, phot_analy_compl, phot_data_compl, err_rm_data,
    completeness, stars_f_acpt, stars_f_rjct, col_0_comb, col_1_comb,
        mag_0_comb, **kwargs):
    """
    Make B2 block plots.
    """
    fig = plt.figure(figsize=(figsize_x, figsize_y))
    gs = gridspec.GridSpec(grid_y, grid_x)
    add_version_plot.main(y_fix=.999)

    # Obtain plotting parameters and data.
    x_ax0, y_ax = prep_plots.ax_names(
        pd['colors'][0], pd['filters'][0], 'mag')
    x_max_cmd0, x_min_cmd0, y_min_cmd0, y_max_cmd0 =\
        prep_plots.diag_limits('mag', col_0_comb, mag_0_comb)
    err_bar_fl0 = prep_plots.error_bars(stars_out_c, x_min_cmd0, err_lst)
    err_bar_cl0 = prep_plots.error_bars(cl_region_c, x_min_cmd0, err_lst)

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

    f_sz_pt = prep_plots.phot_diag_st_size(stars_f_acpt[0])
    cl_sz_pt = prep_plots.phot_diag_st_size(cl_region_c)
    x_min, x_max, y_min, y_max = prep_plots.frame_max_min(
        cld_i['x'], cld_i['y'])
    asp_ratio = prep_plots.aspect_ratio(x_min, x_max, y_min, y_max)
    coord, x_name, y_name = "deg", "ra", "dec"

    # Photometric analysis plots.
    arglist = [
        # pl_cl_fl_regions
        [gs, pd['plot_style'], x_name, y_name, coord, x_min, x_max, y_min,
         y_max, asp_ratio, field_regions_rjct_c, cl_region_rjct_c,
         flag_no_fl_regs_c],
        # pl_cl_diag: Cluster's stars diagram (stars inside cluster's rad)
        [gs, pd['plot_style'], x_ax0, y_ax, x_min_cmd0, x_max_cmd0,
         y_min_cmd0, y_max_cmd0, x_ax1, x_min_cmd1, x_max_cmd1, y_min_cmd1,
         y_max_cmd1, err_bar_cl0, err_bar_cl1, cl_region_rjct_c,
         cl_region_c, n_memb, cl_sz_pt],
        # pl_hess_cmd
        [gs, pd['plot_style'], x_ax0, x_ax1, y_ax, x_max_cmd0, x_min_cmd0,
         y_min_cmd0, y_max_cmd0, x_max_cmd1, x_min_cmd1, y_min_cmd1,
         y_max_cmd1, stars_f_acpt, cl_region_c],
        # pl_fl_diag: Field stars CMD/CCD diagram.
        [gs, pd['plot_style'], x_ax0, y_ax, x_min_cmd0, x_max_cmd0,
         y_min_cmd0, y_max_cmd0, x_ax1, x_min_cmd1, x_max_cmd1, y_min_cmd1,
         y_max_cmd1, field_regions_c, stars_f_rjct, stars_f_acpt, f_sz_pt,
         err_bar_fl0, err_bar_fl1],
        # pl_lum_func: LF of stars in cluster region and outside.
        [gs, pd['plot_style'], y_ax, flag_no_fl_regs_c, lum_func],
        # pl_data_rm_perc
        [gs, pd['plot_style'], y_ax, phot_analy_compl, phot_data_compl,
         err_rm_data, completeness]
    ]
    for n, args in enumerate(arglist):
        mp_data_analysis.plot(n, *args)

    # Generate output file.
    fig.tight_layout()
    plt.savefig(
        join(npd['output_subdir'], str(npd['clust_name']) + '_B2'
             + npd['ext']))
    # Close to release memory.
    plt.clf()
    plt.close("all")

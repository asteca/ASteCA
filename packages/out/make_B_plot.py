
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from os.path import join
from . import mp_phot_analysis
from . import add_version_plot
from . import prep_plots


def main(
        npd, cld_c, pd, em_float, err_lst, cl_region_c, cl_region_rjct_c,
        stars_out_c, stars_out_rjct_c, field_regions_c, flag_no_fl_regs_c,
        field_regions_rjct_c, n_memb, lum_func, phot_analy_compl,
        phot_data_compl, err_rm_data, completeness, **kwargs):
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
        x_ax, y_ax = prep_plots.ax_names(
            pd['colors'][0], pd['filters'][0], 'mag')
        # TODO using first magnitude and color defined
        x_max_cmd, x_min_cmd, y_min_cmd, y_max_cmd = prep_plots.diag_limits(
            'mag', cld_c['cols'][0], cld_c['mags'][0])
        stars_f_rjct, stars_f_acpt = prep_plots.field_region_stars(
            field_regions_c, field_regions_rjct_c)
        f_sz_pt = prep_plots.phot_diag_st_size(len(stars_f_acpt[0]))
        cl_sz_pt = prep_plots.phot_diag_st_size(len(cl_region_c))
        err_bar_fl = prep_plots.error_bars(stars_out_c, x_min_cmd, err_lst)
        err_bar_cl = prep_plots.error_bars(cl_region_c, x_min_cmd, err_lst)
        err_bar_all = prep_plots.error_bars(
            cld_c['mags'][0], x_min_cmd, err_lst, 'all')
        x_min, x_max, y_min, y_max = prep_plots.frame_max_min(
            cld_c['x'], cld_c['y'])
        asp_ratio = prep_plots.aspect_ratio(x_min, x_max, y_min, y_max)
        coord, x_name, y_name = prep_plots.coord_syst(pd['coords'])

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
            [gs, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd, x_ax, y_ax,
                field_regions_c, stars_f_rjct, stars_f_acpt, f_sz_pt,
                err_bar_fl],
            # pl_cl_diag: Cluster's stars diagram (stars inside cluster's rad)
            [gs, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd, x_ax, y_ax,
                cl_region_rjct_c, cl_region_c, n_memb, cl_sz_pt, err_bar_cl],
            # pl_lum_func: LF of stars in cluster region and outside.
            [gs, y_ax, flag_no_fl_regs_c, lum_func],
            # pl_data_rm_perc
            [gs, y_ax, phot_analy_compl, phot_data_compl, err_rm_data,
             completeness]
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
        plt.close("all")

        print("<<Plots for B block created>>")
    else:
        print("<<Skip B block plot>>")


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from os.path import join
from . import add_version_plot
from . import mp_kinem_pms
from . import prep_plots


def main(
    npd, pd, cld_i, PM_flag, pmMP, pmRA_DE, e_pmRA_DE, pmDE, e_pmDE, mmag_pm,
    pmRA_fl_DE, e_pmRA_fl_DE, pmDE_fl, e_pmDE_fl, pm_mag_fl, PM_cl_x, PM_cl_y,
    PM_cl_z, PM_fl_x, PM_fl_y, PM_fl_z, PMs_cl_cx, PMs_cl_cy, PMs_fl_cx,
    PMs_fl_cy, pm_dist_max, PM_flag_all, PM_kde_all, pmRA_all, pmDE_all,
    PMs_d_median, pmMag_all, xRA_all, yDE_all, kde_cent, clust_rad,
        flag_no_fl_regs_i, **kwargs):
    '''
    Make C3 block plots.
    '''

    if 'C3' in pd['flag_make_plot']:

        if PM_flag is False:
            print("  WARNING: nothing to plot in 'C3' block")
            print("<<Skip C3 block plot>>")
            return

        fig = plt.figure(figsize=(30, 25))
        gs = gridspec.GridSpec(10, 12)
        add_version_plot.main(y_fix=.999)

        coord, x_name, y_name = prep_plots.coord_syst(pd['coords'])
        # Uses first magnitude and color defined
        x_ax, y_ax = prep_plots.ax_names(
            pd['colors'][0], pd['filters'][0], 'mag')

        if PM_flag_all:
            arglist = [
                # pms_KDE_all
                [gs, coord, y_ax, pmRA_all, pmDE_all, pmMag_all, PM_kde_all,
                 pd['PM_KDE_std']],
                # pms_NN_all
                [gs, coord, y_ax, pmRA_all, pmDE_all, pmMag_all,
                 PMs_d_median, pd['PM_nnmax'], pd['PM_nnperc']],
                # pms_coords_all
                [fig, gs, cld_i, y_ax, x_name, y_name, coord, PMs_d_median,
                 pd['PM_nnperc'], xRA_all, yDE_all, pmMag_all, kde_cent,
                 clust_rad]
            ]
        for n, args in enumerate(arglist):
            mp_kinem_pms.plot(n, *args)

        # PMs data.
        pmMP, pmRA_DE, e_pmRA_DE, pmDE, e_pmDE, mmag_pm, pm_dist_max =\
            prep_plots.PMsPlot(
                pmMP, pmRA_DE, e_pmRA_DE, pmDE, e_pmDE, mmag_pm, pm_dist_max)
        raPMrng, dePMrng = prep_plots.PMsrange(pmRA_DE, pmDE)
        Nsigma = 2.
        PMs_cent, PMs_width, PMs_height, PMs_theta = prep_plots.SigmaEllipse(
            np.array([pmRA_DE, pmDE]).T, Nsigma)

        arglist = [
            # pms_vpd_mag
            [gs, coord, y_ax, pmRA_DE, pmDE, mmag_pm, pmRA_fl_DE,
             pmDE_fl, pm_mag_fl, raPMrng, dePMrng, flag_no_fl_regs_i],
            # pms_KDE_diag
            [gs, coord, PM_cl_x, PM_cl_y, PM_cl_z, PM_fl_x,
             PM_fl_y, PM_fl_z, PMs_cl_cx, PMs_cl_cy, PMs_fl_cx, PMs_fl_cy,
             raPMrng, dePMrng, PMs_cent, PMs_width, PMs_height, PMs_theta,
             Nsigma],
            # pms_vpd_mp
            [gs, coord, pmMP, pmRA_DE, e_pmRA_DE, pmDE, e_pmDE,
             pmRA_fl_DE, e_pmRA_fl_DE, pmDE_fl, e_pmDE_fl, raPMrng,
             dePMrng],
            # pms_vs_mag
            [gs, coord, y_ax, pmMP, pmRA_DE, pmDE, mmag_pm,
             pmRA_fl_DE, pmDE_fl, pm_mag_fl, raPMrng, dePMrng],
            # pms_dist
            [gs, y_ax, pmMP, pm_dist_max, mmag_pm]
        ]
        for n, args in enumerate(arglist):
            mp_kinem_pms.plot(n + 3, *args)

        # Generate output file.
        fig.tight_layout()
        plt.savefig(
            join(npd['output_subdir'], str(npd['clust_name']) +
                 '_C3.' + pd['plot_frmt']), dpi=pd['plot_dpi'],
            bbox_inches='tight')
        print("<<Plots for C3 block created>>")
        # Close to release memory.
        plt.clf()
        plt.close("all")
    else:
        print("<<Skip C3 block plot>>")

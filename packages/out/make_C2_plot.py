
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from os.path import join
from . import add_version_plot
from . import mp_decont_kinem
from . import prep_plots


def main(
    npd, pd, col_0_comb, mag_0_comb, cl_reg_clean_plot, plx_flag, plx_clrg,
    mmag_clp, mp_clp, plx_clp, e_plx_clp, flag_no_fl_regs_i, field_regions_i,
    PM_flag, pmMP, pmRA_DE, e_pmRA_DE, pmDE, e_pmDE, mmag_pm, pmRA_Bys,
        pmDE_Bys, plx_Bys, plx_wa, cl_reg_fit, **kwargs):
    '''
    Make C2 block plots.
    '''

    if 'C2' in pd['flag_make_plot']:

        if plx_flag is False and PM_flag is False:
            print("<<WARNING: nothing to plot in 'C2' block>>")
            return

        fig = plt.figure(figsize=(30, 25))
        gs = gridspec.GridSpec(10, 12)
        add_version_plot.main()

        coord, x_name, y_name = prep_plots.coord_syst(pd['coords'])
        # # Uses first magnitude and color defined
        x_max_cmd, x_min_cmd, y_min_cmd, y_max_cmd = prep_plots.diag_limits(
            'mag', col_0_comb, mag_0_comb)
        x_ax, y_ax = prep_plots.ax_names(
            pd['colors'][0], pd['filters'][0], 'mag')

        # Decontamination algorithm plots.
        min_prob, bin_edges = cl_reg_clean_plot
        # Parallax data.
        plx_x_kde, kde_pl, plx_flrg, mmag_clp, mp_clp, plx_clp, e_plx_clp =\
            prep_plots.plxPlot(
                plx_flag, plx_clrg, mmag_clp, mp_clp, plx_clp, e_plx_clp,
                flag_no_fl_regs_i, field_regions_i)
        # PMs data.
        pmMP, pmRA_DE, e_pmRA_DE, pmDE, e_pmDE, mmag_pm, pmRA_fl_DE,\
            e_pmRA_fl_DE, pmDE_fl, e_pmDE_fl, x_clpm, y_clpm, z_clpm,\
            pm_dist_max, x_flpm, y_flpm, z_flpm = prep_plots.PMsPlot(
                PM_flag, pmMP, pmRA_DE, e_pmRA_DE, pmDE, e_pmDE,
                mmag_pm, pmRA_Bys, pmDE_Bys, coord, flag_no_fl_regs_i,
                field_regions_i)

        arglist = [
            # plx_histo
            [gs, plx_flag, plx_clrg, plx_x_kde, kde_pl, plx_flrg,
             flag_no_fl_regs_i],
            # plx_chart
            [gs, plx_flag, x_name, y_name, coord, cl_reg_fit, plx_Bys],
            # plx_vs_mag
            [gs, y_min_cmd, y_max_cmd, y_ax, plx_flag, mmag_clp,
             mp_clp, plx_clp, e_plx_clp, plx_Bys, plx_wa],
            # pms_vpd
            [gs, coord, plx_flag, PM_flag, pmMP, pmRA_DE, e_pmRA_DE, pmDE,
             e_pmDE, pmRA_fl_DE, e_pmRA_fl_DE, pmDE_fl, e_pmDE_fl],
            # pms_KDE_diag
            [gs, coord, plx_flag, PM_flag, pmRA_DE, pmDE, x_clpm,
             y_clpm, z_clpm, x_flpm, y_flpm, z_flpm, pmRA_Bys, pmDE_Bys],
            # pms_vs_MP
            [gs, y_ax, plx_flag, PM_flag, pmMP, pm_dist_max, mmag_pm]
        ]
        for n, args in enumerate(arglist):
            mp_decont_kinem.plot(n, *args)

        # Generate output file.
        fig.tight_layout()
        plt.savefig(
            join(npd['output_subdir'], str(npd['clust_name']) +
                 '_C2.' + pd['plot_frmt']), dpi=pd['plot_dpi'],
            bbox_inches='tight')
        print("<<Plots for C2 block created>>")
        # Close to release memory.
        plt.clf()
        plt.close("all")
    else:
        print("<<Skip C2 block plot>>")

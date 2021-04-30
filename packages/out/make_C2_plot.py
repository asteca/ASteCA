
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from os.path import join
from . import add_version_plot
from . import mp_kinem_plx
from . import prep_plots
from .. import aux_funcs
from . prep_plots import figsize_x, figsize_y, grid_x, grid_y


def main(npd, pd, clp):
    """
    Make C2 block plots.
    """
    if clp['plx_flag_clp'] is False:
        print("  WARNING: nothing to plot in 'C2' block")
        print("<<Skip C2 plot>>")
        return

    fig = plt.figure(figsize=(figsize_x, figsize_y))
    gs = gridspec.GridSpec(grid_y, grid_x)
    add_version_plot.main()

    coord, x_name, y_name = prep_plots.coord_syst(pd['coords'])
    x_max_cmd, x_min_cmd, y_min_cmd, y_max_cmd = prep_plots.diag_limits(
        'mag', clp['col_0_comb'], clp['mag_0_comb'])
    x_ax, y_ax = prep_plots.ax_names(
        pd['colors'][0], pd['filters'][0], 'mag')

    # Parallax data.
    plx_flrg, mag_flrg = prep_plots.plxPlot(
        clp['flag_no_fl_regs_i'], clp['field_regions_i'])
    plx_cl_kde_x, plx_cl_kde = aux_funcs.kde1D(clp['plx_clrg'])

    arglist = [
        # plx_histo
        [gs, pd['plot_style'], pd['plx_offset'], clp['plx_clrg'], plx_cl_kde_x,
         plx_cl_kde, plx_flrg, clp['flag_no_fl_regs_i']],
        # plx_chart
        [gs, pd['plot_style'], x_name, y_name, coord, clp['cl_reg_fit'],
         clp['plx_Bys']],
        # plx_vs_mag
        [gs, pd['plot_style'], y_min_cmd, y_max_cmd, y_ax, clp['mmag_clp'],
         clp['plx_clp'], clp['e_plx_clp'], plx_flrg, mag_flrg,
         clp['plx_Bys'], clp['plx_wa']],
        # plx_bys_params
        [gs, clp['plx_bayes_flag_clp'], clp['plx_samples'],
         clp['plx_Bayes_kde'], clp['plx_Bys'], clp['plx_tau_autocorr'],
         clp['mean_afs'], clp['plx_ess']]
    ]
    for n, args in enumerate(arglist):
        mp_kinem_plx.plot(n, *args)

    # DEPRECATED 04/2021
    # if clp['plx_pm_flag']:
    #     # PMs data.
    #     raPMrng, dePMrng = prep_plots.PMsrange(
    #         clp['clreg_PMs']['pmRA'], clp['clreg_PMs']['pmDE'],)

    #     arglist = [
    #         # pms_vs_plx_mp_mag
    #         gs, pd['plot_style'], coord, pd['cosDE_flag'], y_ax,
    #         clp['plx_bayes_flag_clp'], clp['plx_clp'], clp['plx_Bys'],
    #         clp['clreg_PMs'], clp['fregs_PMs'], clp['pm_Plx_cl'],
    #         clp['pm_Plx_fr'], raPMrng, dePMrng]
    #     mp_kinem_plx.plot(4, *arglist)

    # Generate output file.
    fig.tight_layout()
    plt.savefig(
        join(npd['output_subdir'], str(npd['clust_name']) + '_C2'
             + npd['ext']))
    # Close to release memory.
    plt.clf()
    plt.close("all")


import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from os.path import join
from . import add_version_plot
from . import mp_kinem_plx
from . import prep_plots
from .. import aux_funcs
from . prep_plots import figsize_x, figsize_y, grid_x, grid_y


def main(
    npd, pd, col_0_comb, mag_0_comb, plx_flag_clp, plx_clrg,
    mmag_clp, mp_clp, plx_clp, e_plx_clp, flag_no_fl_regs_i, field_regions_i,
    cl_reg_fit, plx_bayes_flag_clp, plx_samples, plx_Bayes_kde, plx_Bys,
    plx_tau_autocorr, mean_afs, plx_ess, plx_wa, plx_pm_flag, clreg_PMs,
        fregs_PMs, pm_Plx_cl, pm_Plx_fr, **kwargs):
    """
    Make C2 block plots.
    """
    if plx_flag_clp is False:
        print("  WARNING: nothing to plot in 'C2' block")
        print("<<Skip C2 plot>>")
        return

    fig = plt.figure(figsize=(figsize_x, figsize_y))
    gs = gridspec.GridSpec(grid_y, grid_x)
    add_version_plot.main()

    coord, x_name, y_name = prep_plots.coord_syst(pd['coords'])
    x_max_cmd, x_min_cmd, y_min_cmd, y_max_cmd = prep_plots.diag_limits(
        'mag', col_0_comb, mag_0_comb)
    x_ax, y_ax = prep_plots.ax_names(
        pd['colors'][0], pd['filters'][0], 'mag')

    # Parallax data.
    plx_flrg, mag_flrg, mmag_clp, mp_clp, plx_clp, e_plx_clp =\
        prep_plots.plxPlot(
            mmag_clp, mp_clp, plx_clp, e_plx_clp, flag_no_fl_regs_i,
            field_regions_i)
    plx_cl_kde_x, plx_cl_kde = aux_funcs.kde1D(plx_clrg)

    arglist = [
        # plx_histo
        [gs, pd['plot_style'], pd['plx_offset'], plx_clrg, plx_cl_kde_x,
         plx_cl_kde, plx_flrg, flag_no_fl_regs_i],
        # plx_chart
        [gs, pd['plot_style'], x_name, y_name, coord, cl_reg_fit, plx_Bys],
        # plx_vs_mag
        [gs, pd['plot_style'], y_min_cmd, y_max_cmd, y_ax, mmag_clp,
         mp_clp, plx_clp, e_plx_clp, plx_flrg, mag_flrg, plx_Bys, plx_wa],
        # plx_bys_params
        [gs, plx_bayes_flag_clp, plx_samples, plx_Bayes_kde, plx_Bys,
         plx_tau_autocorr, mean_afs, plx_ess]
    ]
    for n, args in enumerate(arglist):
        mp_kinem_plx.plot(n, *args)

    if plx_pm_flag:
        # PMs data.
        raPMrng, dePMrng = prep_plots.PMsrange(
            clreg_PMs['pmRA'], clreg_PMs['pmDE'],)

        arglist = [
            # pms_vs_plx_mp_mag
            gs, pd['plot_style'], coord, pd['cosDE_flag'], y_ax,
            plx_bayes_flag_clp, plx_clp, plx_Bys, clreg_PMs, fregs_PMs,
            pm_Plx_cl, pm_Plx_fr, raPMrng, dePMrng]
        mp_kinem_plx.plot(4, *arglist)

    # Generate output file.
    fig.tight_layout()
    plt.savefig(
        join(npd['output_subdir'], str(npd['clust_name']) + '_C2'
             + npd['ext']))
    # Close to release memory.
    plt.clf()
    plt.close("all")

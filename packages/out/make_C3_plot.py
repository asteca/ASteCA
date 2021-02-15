
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from os.path import join
from . import add_version_plot
from . import mp_kinem_pms
from . import prep_plots
from . prep_plots import figsize_x, figsize_y, grid_x, grid_y


def main(
    npd, pd, cld_i, PM_flag, clreg_PMs, fregs_PMs, allfr_PMs, cr_KDE_PMs,
    fr_KDE_PMs, allr_KDE_PMs, pm_dist_max, kde_cent, clust_rad,
        flag_no_fl_regs_i, **kwargs):
    """
    Make C3 block plots.
    """
    if PM_flag is False:
        print("  WARNING: nothing to plot in 'C3' block")
        print("<<Skip C3 plot>>")
        return

    fig = plt.figure(figsize=(figsize_x, figsize_y))
    gs = gridspec.GridSpec(grid_y, grid_x)
    add_version_plot.main(y_fix=.999)

    coord, x_name, y_name = prep_plots.coord_syst(pd['coords'])
    # Uses first magnitude and color defined
    x_ax, y_ax = prep_plots.ax_names(
        pd['colors'][0], pd['filters'][0], 'mag')
    x_min, x_max, y_min, y_max = prep_plots.frame_max_min(
        cld_i['x'], cld_i['y'])
    asp_ratio = prep_plots.aspect_ratio(x_min, x_max, y_min, y_max)

    # PMs data.
    raPMrng, dePMrng = prep_plots.PMsrange(
        clreg_PMs['pmRA'], clreg_PMs['pmDE'])
    Nsigma = 2.
    PMs_cent, PMs_width, PMs_height, PMs_theta = prep_plots.SigmaEllipse(
        np.array([clreg_PMs['pmRA'], clreg_PMs['pmDE']]).T, Nsigma)
    xydelta, xyrang = prep_plots.pmRectangle(allfr_PMs)

    if pd['cosDE_flag'] is False and coord == 'px':
        xlabel = r"$\mu_{{\alpha}} \, \mathrm{[mas/yr]}$"
    else:
        xlabel = r"$\mu_{{\alpha}} \, cos \delta \, \mathrm{[mas/yr]}$"

    arglist = [
        # pms_VPD_all
        [gs, pd['plot_style'], xlabel, pd['PM_KDE_std'], coord, y_ax,
         allfr_PMs],
        # pms_VPD_KDE_all
        [gs, xlabel, coord, y_ax, allr_KDE_PMs, xydelta, xyrang],
        # pms_coords_all
        [fig, gs, pd['plot_style'], coord, x_min, x_max, y_min, y_max,
         asp_ratio, x_name, y_name, kde_cent, clust_rad, allfr_PMs,
         allr_KDE_PMs, xydelta],
        # pms_VPD_zoom
        [gs, pd['plot_style'], xlabel, coord, y_ax, clreg_PMs, fregs_PMs,
         raPMrng, dePMrng, flag_no_fl_regs_i],
        # pms_VPD_zoom_KDE
        [gs, pd['plot_style'], xlabel, coord, cr_KDE_PMs, fr_KDE_PMs,
         raPMrng, dePMrng, PMs_cent, PMs_width, PMs_height, PMs_theta,
         Nsigma],
        # # pms_VPD_zoom_MP
        [gs, pd['plot_style'], xlabel, coord, clreg_PMs, fregs_PMs,
         raPMrng, dePMrng],
        # pms_vs_mag
        [gs, pd['plot_style'], xlabel, coord, y_ax, clreg_PMs, fregs_PMs,
         raPMrng, dePMrng],
        # pms_dist
        [gs, pd['plot_style'], y_ax, clreg_PMs, pm_dist_max]
    ]
    for n, args in enumerate(arglist):
        mp_kinem_pms.plot(n, *args)

    # Generate output file.
    fig.tight_layout()
    plt.savefig(
        join(npd['output_subdir'], str(npd['clust_name']) + '_C3'),
        bbox_inches='tight')
    # Close to release memory.
    plt.clf()
    plt.close("all")

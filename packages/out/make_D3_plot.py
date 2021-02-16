
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from os.path import join
import warnings
from . import add_version_plot
from . import mp_bestfit_CMD
from . import prep_plots
from . prep_plots import figsize_x, figsize_y, grid_x, grid_y


def main(npd, pd, clp):
    """
    Make D3 block plots.
    """
    fig = plt.figure(figsize=(figsize_x, figsize_y))
    gs = gridspec.GridSpec(grid_y, grid_x)
    add_version_plot.main(y_fix=1.005)

    fit_pars = clp['isoch_fit_params']

    best_sol = fit_pars[pd['D3_sol'] + '_sol']

    shift_isoch, synthcl_Nsigma = prep_plots.isoch_sigmaNreg(
        pd['fundam_params'], pd['R_V'], pd['D3_sol'], pd['theor_tracks'],
        pd['m_ini_idx'], pd['ext_coefs'], pd['N_fc'], pd['ext_unif_rand'],
        fit_pars, clp['isoch_fit_errors'], clp['syntClustArgs'])

    # Plot one ore more rows of CMDs/CCDs.
    bf_bin_edges = clp['obs_clust'][0]
    hr_diags = prep_plots.packData(
        pd['lkl_method'], pd['colors'], pd['filters'], clp['cl_max_mag'],
        clp['synth_cl_phot'], clp['binar_idx'], clp['col_0_comb'],
        clp['mag_0_comb'], clp['col_1_comb'], bf_bin_edges,
        shift_isoch, synthcl_Nsigma)
    for (x_phot_all, y_phot_all, x_phot_obs, y_phot_obs, x_synth_phot,
         y_synth_phot, binar_idx, hess_xedges, hess_yedges, x_isoch,
         y_isoch, phot_Nsigma, x_name, y_name, yaxis, i_obs_x,
         i_obs_y, gs_y1, gs_y2) in hr_diags:

        hess_x, hess_y, HD = prep_plots.get_hess(
            [x_phot_obs, y_phot_obs], [x_synth_phot, y_synth_phot],
            hess_xedges, hess_yedges)
        x_ax, y_ax = prep_plots.ax_names(x_name, y_name, yaxis)
        x_max_cmd, x_min_cmd, y_min_cmd, y_max_cmd =\
            prep_plots.diag_limits(yaxis, x_phot_all, y_phot_all)
        sy_sz_pt = prep_plots.phot_diag_st_size(x_synth_phot)

        arglist = [
            # pl_hess_diag: Hess diagram 'observed - synthetic'
            [gs, gs_y1, gs_y2, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd,
             x_ax, y_ax, pd['lkl_method'], hess_xedges, hess_yedges,
             hess_x, hess_y, HD],
            # pl_bf_synth_cl: Best fit synthetic cluster obtained.
            [gs, gs_y1, gs_y2, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd,
             x_ax, y_ax, hess_xedges, hess_yedges, x_synth_phot,
             y_synth_phot, sy_sz_pt, binar_idx, pd['IMF_name'], pd['R_V'],
             best_sol, clp['isoch_fit_errors'], x_isoch, y_isoch,
             pd['lkl_method'], pd['lkl_binning'], pd['evol_track'],
             pd['D3_sol']]
        ]
        for n, args in enumerate(arglist):
            mp_bestfit_CMD.plot(n, *args)

        v_min_mp, v_max_mp = prep_plots.da_colorbar_range(
            clp['cl_max_mag'], [])
        diag_fit_inv, dummy = prep_plots.da_phot_diag(clp['cl_max_mag'], [])
        cl_sz_pt = prep_plots.phot_diag_st_size(diag_fit_inv)
        # Main photometric diagram of observed cluster.
        i_y = 0 if yaxis == 'mag' else 1
        # x axis is always a color so this the index is fixed to '1'.
        # y axis is not, so the 'i_y' index determines what goes there.
        obs_x, obs_y, obs_MPs = diag_fit_inv[1][i_obs_x],\
            diag_fit_inv[i_y][i_obs_y], diag_fit_inv[2]
        # tight_layout is called here
        plot_observed_cluster(
            fig, gs, gs_y1, gs_y2, x_ax, y_ax, clp['cl_max_mag'], x_min_cmd,
            x_max_cmd, y_min_cmd, y_max_cmd, clp['err_lst'], v_min_mp,
            v_max_mp, obs_x, obs_y, obs_MPs, cl_sz_pt, hess_xedges,
            hess_yedges, x_isoch, y_isoch, phot_Nsigma, pd['lkl_method'])

    # Generate output file.
    plt.savefig(
        join(npd['output_subdir'], str(npd['clust_name']) + '_D3'
             + npd['ext']))
    # Close to release memory.
    plt.clf()
    plt.close("all")


def plot_observed_cluster(
    fig, gs, gs_y1, gs_y2, x_ax, y_ax, cl_max_mag, x_min_cmd, x_max_cmd,
    y_min_cmd, y_max_cmd, err_lst, v_min_mp, v_max_mp, obs_x, obs_y, obs_MPs,
    cl_sz_pt, hess_xedges, hess_yedges, x_isoch, y_isoch, phot_Nsigma,
        lkl_method):
    """
    This function is called separately since we need to retrieve some
    information from it to plot that #$%&! colorbar.
    """
    err_bar = prep_plots.error_bars(cl_max_mag, x_min_cmd, err_lst)

    # pl_mps_phot_diag
    plot_colorbar, sca, trans = mp_bestfit_CMD.pl_mps_phot_diag(
        gs, gs_y1, gs_y2, fig, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd,
        x_ax, y_ax, v_min_mp, v_max_mp, obs_x, obs_y, obs_MPs, err_bar,
        cl_sz_pt, hess_xedges, hess_yedges, x_isoch, y_isoch, phot_Nsigma,
        lkl_method)

    # Ignore warning issued by colorbar plotted in photometric diagram with
    # membership probabilities.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fig.tight_layout()

    # Force to not plot colorbar after the first row.
    plot_colorbar = False if gs_y1 != 0 else plot_colorbar

    # Plot colorbar down here so tight_layout won't move it around.
    if plot_colorbar is True:
        import matplotlib.transforms as mts
        # Position and dimensions relative to the axes.
        x0, y0, width, height = [0.74, 0.93, 0.2, 0.04]
        # Transform them to get the ABSOLUTE POSITION AND DIMENSIONS
        Bbox = mts.Bbox.from_bounds(x0, y0, width, height)
        l, b, w, h = mts.TransformedBbox(Bbox, trans).bounds
        # Create the axes and the colorbar.
        cbaxes = fig.add_axes([l, b, w, h])
        cbar = plt.colorbar(
            sca, cax=cbaxes, ticks=[v_min_mp, v_max_mp],
            orientation='horizontal')
        cbar.ax.minorticks_off()
        # cbar.ax.set_title('MP', x=1.2, y=-.5)


import itertools
import warnings
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from os.path import join
from . import add_version_plot
from . import mp_decont_photom
from . import prep_plots
from . prep_plots import figsize_x, figsize_y, grid_x, grid_y


def main(
    npd, cld_c, pd, kde_cent, clust_rad, stars_out_c, cl_region_i,
    memb_probs_cl_region_i, memb_prob_avrg_sort, flag_decont_skip, n_memb_da,
    cl_reg_fit, cl_reg_no_fit, local_rm_edges, col_0_comb, mag_0_comb,
        err_lst, **kwargs):
    """
    Make C1 block plots.
    """
    fig = plt.figure(figsize=(figsize_x, figsize_y))
    gs = gridspec.GridSpec(grid_y, grid_x)
    add_version_plot.main(y_fix=.999)

    # Obtain plotting parameters and data.
    x_min, x_max, y_min, y_max = prep_plots.frame_max_min(
        cld_c['x'], cld_c['y'])
    x_zmin, x_zmax, y_zmin, y_zmax = prep_plots.frame_zoomed(
        x_min, x_max, y_min, y_max, kde_cent, clust_rad)
    coord, x_name, y_name = prep_plots.coord_syst(pd['coords'])
    v_min_mp_comp, v_max_mp_comp = prep_plots.da_colorbar_range(
        cl_reg_fit, cl_reg_no_fit)
    chart_fit_inv, chart_no_fit_inv, out_clust_rad =\
        prep_plots.da_find_chart(
            kde_cent, clust_rad, stars_out_c, x_zmin, x_zmax, y_zmin,
            y_zmax, cl_reg_fit, cl_reg_no_fit)

    # Decontamination algorithm plots.
    arglist = [
        # pl_mp_histo
        [gs, pd['plot_style'], n_memb_da, memb_prob_avrg_sort,
         flag_decont_skip, cl_reg_fit, pd['fld_clean_mode'],
         pd['fld_clean_bin']],
        # pl_chart_mps
        [gs, fig, x_name, y_name, coord, x_zmin, x_zmax, y_zmin,
         y_zmax, kde_cent, clust_rad, flag_decont_skip,
         v_min_mp_comp, v_max_mp_comp, chart_fit_inv, chart_no_fit_inv,
         out_clust_rad, pd['fld_clean_mode'], pd['fld_clean_bin']],
    ]
    for n, args in enumerate(arglist):
        mp_decont_photom.plot(n, *args)

    # Plot all 2D combinations of magnitudes and colors.
    mags_c, mags_i, cols_c, cols_i, colors_c, colors_i =\
        prep_plots.complSeparate(
            cl_region_i, memb_probs_cl_region_i, memb_prob_avrg_sort)
    mags_all = mags_c + mags_i
    j_gs = 0
    all_colorbars = []
    for i, col_c in enumerate(cols_c):
        col_i = cols_i[i]
        cols_all = list(col_c) + list(col_i)

        x_ax, y_ax = prep_plots.ax_names(
            pd['colors'][i], pd['filters'][0], 'mag')
        x_max_cmd, x_min_cmd, y_min_cmd, y_max_cmd =\
            prep_plots.diag_limits('mag', cols_all, mags_all)
        cl_c_sz_pt = prep_plots.phot_diag_st_size(mags_c)
        cl_i_sz_pt = prep_plots.phot_diag_st_size(mags_i)

        # Color indexes for the 'local_cell_clean' edges
        c1, c2 = i + 1, 0
        plot_colorbar, sca, trans, v_min_mp, v_max_mp =\
            mp_decont_photom.pl_mps_incomp_diags(
                gs, fig, pd['plot_style'], x_min_cmd, x_max_cmd, y_min_cmd,
                y_max_cmd, x_ax, y_ax, col_c, mags_c, colors_c, col_i,
                mags_i, colors_i, cl_c_sz_pt, cl_i_sz_pt,
                pd['fld_clean_mode'], local_rm_edges, c1, c2, j_gs)
        all_colorbars.append((
            plot_colorbar, sca, trans, v_min_mp, v_max_mp))
        j_gs += 1

    c_combs = list(itertools.combinations(range(len(cols_c)), 2))
    for i, j in c_combs:
        cols0_all = list(cols_c[i]) + list(cols_i[i])
        cols1_all = list(cols_c[j]) + list(cols_i[j])
        x_max_cmd, x_min_cmd, y_min_cmd, y_max_cmd =\
            prep_plots.diag_limits('col', cols0_all, cols1_all)
        x_ax, y_ax = prep_plots.ax_names(
            pd['colors'][i], pd['colors'][j], 'col')
        cl_c_sz_pt = prep_plots.phot_diag_st_size(cols_c[i])
        cl_i_sz_pt = prep_plots.phot_diag_st_size(cols_i[i])

        # Color indexes for the 'local_cell_clean' edges
        c1, c2 = i + 1, j + 1
        plot_colorbar, sca, trans, v_min_mp, v_max_mp =\
            mp_decont_photom.pl_mps_incomp_diags(
                gs, fig, pd['plot_style'], x_min_cmd, x_max_cmd, y_min_cmd,
                y_max_cmd, x_ax, y_ax, cols_c[i], cols_c[j], colors_c,
                cols_i[i], cols_i[j], colors_i, cl_c_sz_pt, cl_i_sz_pt,
                pd['fld_clean_mode'], local_rm_edges, c1, c2, j_gs)
        all_colorbars.append((
            plot_colorbar, sca, trans, v_min_mp, v_max_mp))
        j_gs += 1

    diag_fit_inv, diag_no_fit_inv = prep_plots.da_phot_diag(
        cl_reg_fit, cl_reg_no_fit)
    cl_f_sz_pt = prep_plots.phot_diag_st_size(diag_fit_inv)
    cl_nf_sz_pt = prep_plots.phot_diag_st_size(diag_no_fit_inv)
    # Uses first magnitude and color defined
    x_max_cmd, x_min_cmd, y_min_cmd, y_max_cmd = prep_plots.diag_limits(
        'mag', col_0_comb, mag_0_comb)
    x_ax, y_ax = prep_plots.ax_names(
        pd['colors'][0], pd['filters'][0], 'mag')
    err_bar = prep_plots.error_bars(
        cl_reg_fit + cl_reg_no_fit, x_min_cmd, err_lst)
    # pl_mps_phot_diag
    plot_colorbar, sca, trans = mp_decont_photom.pl_mps_phot_diag(
        gs, pd['plot_style'], fig, x_min_cmd, x_max_cmd, y_min_cmd,
        y_max_cmd, x_ax, y_ax, v_min_mp_comp, v_max_mp_comp, diag_fit_inv,
        diag_no_fit_inv, cl_f_sz_pt, cl_nf_sz_pt, err_bar,
        pd['fld_clean_mode'], local_rm_edges)
    all_colorbars.append((
        plot_colorbar, sca, trans, v_min_mp_comp, v_max_mp_comp))

    plot_colorbars(fig, all_colorbars)

    # Generate output file.
    plt.savefig(
        join(npd['output_subdir'], str(npd['clust_name']) + '_C1'
             + npd['ext']))
    # Close to release memory.
    plt.clf()
    plt.close("all")


def plot_colorbars(fig, all_colorbars):
    """
    Plot colorbars *after* tight_layout(), otherwise it will move them around.
    """

    # Ignore warning issued by colorbar.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fig.tight_layout()

    # Plot colorbar down here so tight_layout won't move it around.
    for (plot_colorbar, sca, trans, v_min_mp, v_max_mp) in all_colorbars:
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

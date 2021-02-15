
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from os.path import join
from . import mp_cent_dens
from . import add_version_plot
from . import prep_plots
from . prep_plots import grid_x, grid_y, figsize_x, figsize_y


def main(
    npd, cld_i, pd, x_offset, y_offset, bw_list, kde_cent, frame_kde_cent,
    integ_dists, integ_mags, xy_filtered, xy_cent_dist, NN_dd, NN_dist,
    fr_dens, fdens_min_d, fdens_lst, fdens_std_lst, field_dens_d, field_dens,
        field_dens_std, clust_rad, **kwargs):
    """
    Make A2 block plots.
    """
    fig = plt.figure(figsize=(figsize_x, figsize_y))
    gs = gridspec.GridSpec(grid_y, grid_x)
    add_version_plot.main()

    # Obtain plotting parameters and data.
    x_min, x_max, y_min, y_max = prep_plots.frame_max_min(
        cld_i['x'], cld_i['y'])
    asp_ratio = prep_plots.aspect_ratio(x_min, x_max, y_min, y_max)
    coord, x_name, y_name = prep_plots.coord_syst(pd['coords'])
    st_sizes_arr = prep_plots.star_size(cld_i['mags'][0])
    _, y_ax = prep_plots.ax_names(pd['colors'][0], pd['filters'][0], 'mag')

    # Structure plots.
    arglist = [
        # pl_full_frame: x,y finding chart of full frame.
        [gs, fig, pd['project'], x_offset, y_offset, x_name, y_name, coord,
         x_min, x_max, y_min, y_max, asp_ratio, kde_cent, cld_i['x'],
         cld_i['y'], st_sizes_arr, clust_rad],
        # pl_densmap: 2D Gaussian convolved histogram.
        [gs, fig, asp_ratio, x_name, y_name, coord, bw_list, kde_cent,
         frame_kde_cent, fr_dens, clust_rad],
        # pl_knn_dens
        [gs, fig, pd['plot_style'], asp_ratio, x_min, x_max, y_min, y_max,
         x_name, y_name, coord, NN_dd, xy_filtered, fr_dens, NN_dist,
         pd['project'], x_offset, y_offset, kde_cent, clust_rad],
        # pl_field_dens
        [gs, pd['plot_style'], coord, pd['fdens_method'], xy_cent_dist,
         fr_dens, fdens_min_d, fdens_lst, fdens_std_lst, field_dens_d,
         field_dens, field_dens_std],
        # pl_centdist_vs_mag
        [gs, fig, pd['plot_style'], y_ax, coord, cld_i['x'], cld_i['y'],
         cld_i['mags'][0], kde_cent, clust_rad, integ_dists, integ_mags]
    ]
    for n, args in enumerate(arglist):
        mp_cent_dens.plot(n, *args)

    fig.tight_layout()
    plt.savefig(
        join(npd['output_subdir'], str(npd['clust_name']) + '_A2'),
        bbox_inches='tight')
    # Close to release memory.
    plt.clf()
    plt.close("all")

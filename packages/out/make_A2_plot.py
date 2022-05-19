
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from os.path import join
from . import mp_cent_dens
from . import add_version_plot
from . import prep_plots
from . prep_plots import grid_x, grid_y, figsize_x, figsize_y


def main(npd, cld_i, pd, clp):
    """
    Make A2 block plots.
    """
    fig = plt.figure(figsize=(figsize_x, figsize_y))
    gs = gridspec.GridSpec(grid_y, grid_x)
    add_version_plot.main()

    # Obtain plotting parameters and data.
    x_min, x_max, y_min, y_max = prep_plots.frame_max_min(
        cld_i['x'], cld_i['y'])
    coord = "deg"
    if pd['xy_frame'] == 'equatorial':
        x_name, y_name = "ra", "dec"
        x_min, x_max = x_max, x_min
    else:
        x_name, y_name = "lon", "lat"
    asp_ratio = prep_plots.aspect_ratio(x_min, x_max, y_min, y_max)

    st_sizes_arr = prep_plots.star_size(cld_i['mags'][0])
    _, y_ax = prep_plots.ax_names(pd['colors'][0], pd['filters'][0], 'mag')

    # Structure plots.
    arglist = [
        # pl_full_frame: x,y finding chart of full frame.
        [gs, fig, x_name, y_name, coord, x_min, x_max, y_min, y_max, asp_ratio,
         clp['kde_cent'], cld_i['x'], cld_i['y'], st_sizes_arr,
         clp['clust_rad']],
        # pl_densmap: 2D Gaussian convolved histogram.
        [gs, fig, asp_ratio, x_name, y_name, coord, clp['bw_list'],
         clp['kde_cent'], clp['frame_kde_cent'], clp['pts_dens'],
         clp['clust_rad']],
        # # pl_knn_dens
        # [gs, fig, pd['plot_style'], asp_ratio, x_min, x_max, y_min, y_max,
        #  x_name, y_name, coord, clp['xy_filtered'],
        #  clp['pts_dens'], clp['NN_dist'], clp['kde_cent'], clp['clust_rad']],
        # pl_field_dens
        [gs, pd['plot_style'], coord, pd['fdens_method'], clp['xy_cent_dist'],
         clp['pts_dens'], clp['fdens_min_d'], clp['fdens_lst'],
         clp['fdens_std_lst'], clp['field_dens']],
        # pl_centdist_vs_mag
        [gs, fig, pd['plot_style'], y_ax, coord, cld_i['x'], cld_i['y'],
         cld_i['mags'][0], clp['kde_cent'], clp['clust_rad'],
         clp['integ_dists'], clp['integ_mags']]
    ]
    for n, args in enumerate(arglist):
        mp_cent_dens.plot(n, *args)

    fig.tight_layout()
    fname = join(npd['output_subdir'], npd['clust_name'] + '_A2' + npd['ext'])
    plt.savefig(fname)
    # Close to release memory.
    plt.clf()
    plt.close("all")


import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from os.path import join
from . import mp_cent_dens
from . import add_version_plot
from . import prep_plots


def main(
    npd, cld_i, pd, x_offset, y_offset, bw_list, kde_cent, frame_kde_cent,
    rdp_radii, integ_dists, integ_mags, xy_filtered, xy_cent_dist, NN_dist,
    fr_dens, fdens_min_d, fdens_lst, fdens_std_lst, field_dens_d, field_dens,
        field_dens_std, **kwargs):
    """
    Make A2 block plots.
    """
    if 'A2' in pd['flag_make_plot']:
        # figsize(x1, y1), GridSpec(y2, x2)
        fig = plt.figure(figsize=(30, 25))
        gs = gridspec.GridSpec(10, 12)
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
            # pl_densmap: 2D Gaussian convolved histogram.
            [gs, fig, asp_ratio, x_name, y_name, coord, bw_list, kde_cent,
             frame_kde_cent, fr_dens],
            # pl_knn_dens
            [gs, fig, asp_ratio, x_min, x_max, y_min, y_max, x_name, y_name,
             coord, pd['NN_dd'], xy_filtered, fr_dens, NN_dist, kde_cent],
            # pl_full_frame: x,y finding chart of full frame.
            [gs, fig, pd['project'], x_offset, y_offset, x_name, y_name, coord,
             x_min, x_max, y_min, y_max, asp_ratio, kde_cent, cld_i['x'],
             cld_i['y'], st_sizes_arr],
            # pl_field_dens
            [gs, coord, pd['fdens_method'], xy_cent_dist, fr_dens, fdens_min_d,
             fdens_lst, fdens_std_lst, field_dens_d, field_dens,
             field_dens_std],
            # pl_mag_cent
            [gs, coord, y_ax, integ_dists, integ_mags],
            # pl_rdp_rings
            [gs, fig, asp_ratio, x_min, x_max, y_min, y_max, x_name, y_name,
             coord, kde_cent, rdp_radii]
        ]
        for n, args in enumerate(arglist):
            mp_cent_dens.plot(n, *args)

        fig.tight_layout()
        plt.savefig(
            join(npd['output_subdir'], str(npd['clust_name']) +
                 '_A2.' + pd['plot_frmt']), dpi=pd['plot_dpi'],
            bbox_inches='tight')
        # Close to release memory.
        plt.clf()
        plt.close("all")

        print("<<Plots for A2 block created>>")
    else:
        print("<<Skip A2 block plot>>")

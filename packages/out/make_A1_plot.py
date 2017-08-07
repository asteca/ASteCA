
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from os.path import join
import mp_centers
import add_version_plot
import prep_plots


def main(
        npd, cld, pd, xy_mag_ranges, hist_2d_g, cents_bin_2d, st_dev_lst,
        **kwargs):
    '''
    Make A1 block plots.
    '''
    if 'A1' in pd['flag_make_plot']:
        # figsize(x1, y1), GridSpec(y2, x2) --> To have square plots: x1/x2 =
        # y1/y2 = 2.5
        fig = plt.figure(figsize=(30, 25))
        gs = gridspec.GridSpec(10, 12)
        add_version_plot.main()

        # Obtain plotting parameters and data.
        x_min, x_max, y_min, y_max = prep_plots.frame_max_min(
            cld['x'], cld['y'])
        asp_ratio = prep_plots.aspect_ratio(x_min, x_max, y_min, y_max)
        coord, x_name, y_name = prep_plots.coord_syst(pd['coords'])
        x_ax, y_ax, y_axis = prep_plots.ax_names(pd['filters'], pd['colors'])

        N_all, min_mag_all = len(cld['mags'][0]), min(cld['mags'][0])
        # Structure plots.
        arglist = []
        for mag_rng in xy_mag_ranges:
            x, y, m = zip(*mag_rng.values()[0])
            mag_range = mag_rng.keys()[0]
            st_sizes_arr = prep_plots.star_size(m, N=N_all, min_m=min_mag_all)
            arglist.append(
                # pl_full_frame: x,y finding chart of full frame.
                [gs, fig, x_name, y_name, coord, x_min, x_max, y_min, y_max,
                 asp_ratio, x, y, st_sizes_arr, mag_range, y_ax]
            )
        for n, args in enumerate(arglist):
            # with timeblock("{}".format(n)):
            mp_centers.plot(n, *args)

        # 2D Gaussian convolved histogram.
        arglist = []
        for h2d, cb2d in zip(*[hist_2d_g, cents_bin_2d]):
            arglist.append(
                # pl_center: 2D Gaussian convolved histogram.
                [gs, fig, asp_ratio, x_name, y_name, coord, st_dev_lst, h2d,
                 cb2d],
            )
        for n, args in enumerate(arglist, 5):
            # with timeblock("{}".format(n)):
            mp_centers.plot(n, *args)

        # Generate output file.
        fig.tight_layout()
        plt.savefig(
            join(npd['output_subdir'], str(npd['clust_name']) +
                 '_A1.' + pd['plot_frmt']), dpi=pd['plot_dpi'],
            bbox_inches='tight')
        # Close to release memory.
        plt.clf()
        plt.close()

        print("<<Plots for A1 block created>>")
    else:
        print("<<Skip A1 block plot>>")

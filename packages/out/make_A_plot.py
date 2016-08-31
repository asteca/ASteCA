
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from os.path import join
import warnings
from .._version import __version__
import mp_structure
import prep_plots


#############################################################
# # Timer function: http://stackoverflow.com/a/21860100/1391441
# from contextlib import contextmanager
# import time
# @contextmanager
# def timeblock(label):
#     start = time.clock()
#     try:
#         yield
#     finally:
#         end = time.clock()
#         print ('{} elapsed: {}'.format(label, end - start))
#############################################################

def main(
        npd, cld, pd, clust_cent, e_cent, K_cent_dens, clust_rad, e_rad,
        poisson_error, stars_out_rjct, field_regions, cent_bin, bin_width,
        hist_2d_g, approx_cents, st_dev_lst, core_rad, e_core, tidal_rad,
        e_tidal, K_conct_par, flag_2pk_conver, flag_3pk_conver, radii,
        rdp_points, field_dens, cont_index, kde_plot, cl_region,
        flag_no_fl_regs, **kwargs):
    '''
    Make A block plots.
    '''

    # flag_make_plot = pd['pl_params'][0]
    if pd['pl_params'][0]:

        # figsize(x1, y1), GridSpec(y2, x2) --> To have square plots: x1/x2 =
        # y1/y2 = 2.5
        fig = plt.figure(figsize=(30, 25))  # create the top-level container
        gs = gridspec.GridSpec(10, 12)      # create a GridSpec object
        # Add version number to top left.
        ver = '[ASteCA ' + __version__ + ']'
        x_coord = 0.957 - (len(__version__) - 6) * 0.001
        plt.figtext(x_coord, .988, ver, fontsize=9, color='#585858')

        # Obtain plotting parameters and data.
        x_min, x_max, y_min, y_max = prep_plots.frame_max_min(
            cld['x'], cld['y'])
        asp_ratio = prep_plots.aspect_ratio(x_min, x_max, y_min, y_max)
        coord, x_name, y_name = prep_plots.coord_syst(pd['coords'])
        x_zmin, x_zmax, y_zmin, y_zmax = prep_plots.frame_zoomed(
            x_min, x_max, y_min, y_max, clust_cent, clust_rad)
        x_data_z, y_data_z, mag_data_z = prep_plots.zoomed_frame(
            cld['x'], cld['y'], cld['mags'], x_zmin, x_zmax, y_zmin, y_zmax)
        st_sizes_arr = prep_plots.star_size(cld['mags'][0])
        st_sizes_arr_z = prep_plots.star_size(mag_data_z)

        # Structure plots.
        arglist = [
            # pl_hist_g: 2D Gaussian convolved histogram.
            [gs, fig, asp_ratio, x_name, y_name, coord, cent_bin, clust_rad,
                bin_width, hist_2d_g],
            # pl_centers: 2D Gaussian histograms' centers.
            [gs, x_name, y_name, coord, x_min, x_max, y_min, y_max, asp_ratio,
                approx_cents, bin_width, st_dev_lst],
            # pl_full_frame: x,y finding chart of full frame.
            [gs, fig, x_name, y_name, coord, x_min, x_max, y_min, y_max,
             asp_ratio, clust_cent, clust_rad, e_cent, cld['x'], cld['y'],
             st_sizes_arr, core_rad, e_core, tidal_rad, e_tidal, K_conct_par,
             flag_2pk_conver, flag_3pk_conver],
            # pl_rad_dens: Radial density plot.
            [gs, pd['mode'], radii, rdp_points, field_dens, coord,
             npd['clust_name'], clust_rad, e_rad, poisson_error, bin_width,
             core_rad, e_core, tidal_rad, e_tidal, K_cent_dens,
             flag_2pk_conver, flag_3pk_conver],
            # pl_zoom_frame: Zoom on x,y finding chart.
            [gs, fig, x_name, y_name, coord, x_zmin, x_zmax, y_zmin, y_zmax,
                cont_index, kde_plot, x_data_z, y_data_z,
                st_sizes_arr_z, clust_cent, clust_rad],
            # pl_cl_fl_regions: Cluster and field regions defined.
            [gs, fig, x_name, y_name, coord, x_min, x_max, y_min, y_max,
                asp_ratio, clust_cent, clust_rad, field_regions, cl_region,
                flag_no_fl_regs]
        ]
        for n, args in enumerate(arglist):
            # with timeblock("{}".format(n)):
            mp_structure.plot(n, *args)

        # Ignore warning issued by colorbar plotted in photometric diagram with
        # membership probabilities.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fig.tight_layout()

        # Generate output file for each data file.
        pl_fmt, pl_dpi = pd['pl_params'][1:3]
        plt.savefig(
            join(npd['output_subdir'], str('A_' + npd['clust_name']) +
                 '.' + pl_fmt), dpi=pl_dpi, bbox_inches='tight')

        # Close to release memory.
        plt.clf()
        plt.close()

        print("Plots from block 'A' created.")

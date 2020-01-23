
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from os.path import join
import warnings
from . import mp_structure
from . import add_version_plot
from . import prep_plots


def main(
    npd, cld_i, pd, x_offset, y_offset, bw_list, kde_cent, kde_plot,
    K_cent_dens, clust_rad, e_rad, rad_rads, rad_N_membs, rad_N_field, rad_CI,
    bin_cent, bin_width, kde_approx_cent, frame_kde_cent, core_rad, e_core,
    tidal_rad, e_tidal, K_conct_par, flag_2pk_conver, flag_3pk_conver,
    rdp_radii, rdp_points, rdp_stddev, xy_filtered, xy_cent_dist, NN_dist,
    fr_dens, fdens_min_d, fdens_lst, fdens_std_lst, field_dens_d, field_dens,
    field_dens_std, cont_index, n_memb_i, cl_region_i, frac_cl_area,
    cl_region_rjct_i, field_regions_rjct_i, field_regions_i, flag_no_fl_regs_i,
        **kwargs):
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
        x_zmin, x_zmax, y_zmin, y_zmax = prep_plots.frame_zoomed(
            x_min, x_max, y_min, y_max, kde_cent, clust_rad)
        x_data_z, y_data_z, mag_data_z = prep_plots.zoomed_frame(
            cld_i['x'], cld_i['y'], cld_i['mags'], x_zmin, x_zmax, y_zmin,
            y_zmax)
        st_sizes_arr = prep_plots.star_size(cld_i['mags'][0])
        st_sizes_arr_z = prep_plots.star_size(mag_data_z)

        # Structure plots.
        arglist = [
            # pl_center: 2D Gaussian convolved histogram.
            [gs, fig, asp_ratio, x_name, y_name, coord, bw_list, kde_cent,
             frame_kde_cent, fr_dens, clust_rad],
            # pl_knn_dens
            [gs, fig, asp_ratio, x_min, x_max, y_min, y_max, x_name, y_name,
             coord, pd['NN_dd'], xy_filtered, fr_dens, NN_dist, kde_cent,
             rdp_radii],
            # pl_field_dens
            [gs, coord, pd['fdens_method'], xy_cent_dist, fr_dens, fdens_min_d,
             fdens_lst, fdens_std_lst, field_dens_d, field_dens,
             field_dens_std],
            # pl_rad_find
            [gs, coord, clust_rad, e_rad, rad_rads, rad_N_membs, rad_N_field,
             rad_CI],
            # pl_rad_dens: Radial density plot.
            [gs, coord, rdp_radii, rdp_points, rdp_stddev, field_dens,
             field_dens_std, clust_rad, e_rad, core_rad,
             e_core, tidal_rad, e_tidal, K_cent_dens, flag_2pk_conver,
             flag_3pk_conver],
            # pl_full_frame: x,y finding chart of full frame.
            [gs, fig, pd['project'], x_offset, y_offset, x_name, y_name,
             coord, x_min, x_max, y_min, y_max, asp_ratio, kde_cent, clust_rad,
             frac_cl_area, cld_i['x'], cld_i['y'], st_sizes_arr, core_rad,
             e_core, tidal_rad, e_tidal, K_conct_par, flag_2pk_conver,
             flag_3pk_conver],
            # pl_zoom_frame: Zoom on x,y finding chart.
            [gs, fig, x_name, y_name, coord, x_zmin, x_zmax, y_zmin, y_zmax,
                cont_index, n_memb_i, kde_plot, x_data_z, y_data_z,
                st_sizes_arr_z, kde_cent, clust_rad],
            # pl_cl_fl_regions: Cluster and field regions defined.
            [gs, fig, x_name, y_name, coord, x_min, x_max, y_min, y_max,
                asp_ratio, kde_cent, clust_rad, field_regions_i,
                field_regions_rjct_i, cl_region_i, cl_region_rjct_i,
                flag_no_fl_regs_i]
        ]
        for n, args in enumerate(arglist):
            mp_structure.plot(n, *args)

        # Generate output file.
        # Ignore warning issued by log-log inset.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
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

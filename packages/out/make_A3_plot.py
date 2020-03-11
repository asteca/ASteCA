
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from os.path import join
import warnings
from . import mp_radius
from . import add_version_plot
from . import prep_plots
from . prep_plots import figsize_x, figsize_y, grid_x, grid_y


def main(
    npd, cld_i, pd, kde_cent, clust_rad, e_rad, rad_rads, rad_N_membs,
    rad_N_field, rad_CI, KP_cent_dens, membvsmag,
    KP_steps, KP_mean_afs, KP_tau_autocorr, KP_ESS, KP_samples, KP_Bys_rc,
    KP_Bys_rt, KP_memb_num, KP_conct_par, rdp_radii, rdp_points, rdp_stddev,
    field_dens, field_dens_std, cont_index, cl_region_i, frac_cl_area,
    cl_region_rjct_i, field_regions_rjct_i, field_regions_i,
        flag_no_fl_regs_i, **kwargs):
    """
    Make A3 block plots.
    """
    if 'A3' in pd['flag_make_plot']:
        # fig = plt.figure(figsize=(figsize_x, figsize_y))
        # gs = gridspec.GridSpec(grid_y, grid_x)
        fig = plt.figure(figsize=(25, 25))
        gs = gridspec.GridSpec(10, 12)

        add_version_plot.main()

        # Obtain plotting parameters and data.
        x_min, x_max, y_min, y_max = prep_plots.frame_max_min(
            cld_i['x'], cld_i['y'])
        asp_ratio = prep_plots.aspect_ratio(x_min, x_max, y_min, y_max)
        coord, x_name, y_name = prep_plots.coord_syst(pd['coords'])
        x_zmin, x_zmax, y_zmin, y_zmax = prep_plots.frame_zoomed(
            x_min, x_max, y_min, y_max, kde_cent, clust_rad, pd['kp_flag'],
            KP_Bys_rt)
        x_data_z, y_data_z, mag_data_z = prep_plots.zoomed_frame(
            cld_i['x'], cld_i['y'], cld_i['mags'], x_zmin, x_zmax, y_zmin,
            y_zmax)
        st_sizes_arr_z = prep_plots.star_size(mag_data_z)
        _, y_ax = prep_plots.ax_names(pd['colors'][0], pd['filters'][0], 'mag')

        # Structure plots.
        arglist = [
            # pl_rad_find
            [gs, coord, clust_rad, e_rad, rad_rads, rad_N_membs, rad_N_field,
             rad_CI],
            # pl_rad_dens: Radial density plot.
            [gs, coord, rdp_radii, rdp_points, rdp_stddev, field_dens,
             field_dens_std, clust_rad, e_rad, pd['kp_flag'], KP_Bys_rc,
             KP_Bys_rt, KP_cent_dens, KP_conct_par],
            # pl_KP_Bys
            [gs, coord, pd['kp_flag'], pd['kp_nburn'], KP_steps, KP_mean_afs,
             KP_tau_autocorr, KP_ESS, KP_samples, KP_Bys_rc, KP_Bys_rt],
            # pl_zoom_frame: Zoom on x,y finding chart.
            [gs, fig, x_name, y_name, coord, x_zmin, x_zmax, y_zmin, y_zmax,
             cont_index, x_data_z, y_data_z, st_sizes_arr_z,
             kde_cent, clust_rad, KP_Bys_rc, KP_Bys_rt, frac_cl_area,
             pd['kp_flag']],
            # pl_mag_membs
            [gs, y_ax, pd['kp_flag'], membvsmag],
            # pl_cl_fl_regions: Cluster and field regions defined.
            [gs, fig, x_name, y_name, coord, x_min, x_max, y_min, y_max,
             asp_ratio, kde_cent, clust_rad, field_regions_i,
             field_regions_rjct_i, cl_region_i, cl_region_rjct_i,
             flag_no_fl_regs_i, pd['kp_flag']]
        ]
        for n, args in enumerate(arglist):
            mp_radius.plot(n, *args)

        # Generate output file.
        # Ignore warning issued by log-log inset.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fig.tight_layout()
        plt.savefig(
            join(npd['output_subdir'], str(npd['clust_name']) +
                 '_A3.' + pd['plot_frmt']), dpi=pd['plot_dpi'],
            bbox_inches='tight')
        # Close to release memory.
        plt.clf()
        plt.close("all")

        print("<<Plots for A3 block created>>")
    else:
        print("<<Skip A3 block plot>>")

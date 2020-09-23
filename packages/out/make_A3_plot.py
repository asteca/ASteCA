
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from os.path import join
import warnings
from . import mp_radius, mp_KP_bayes
from . import add_version_plot
from . import prep_plots
from . prep_plots import figsize_x, figsize_y, grid_x, grid_y


def main(
    npd, cld_i, pd, kde_cent, frame_kde_cent, clust_rad, e_rad, rad_rads,
    rad_N_membs, rad_N_field, rad_CI, membvsmag, xy_filtered, xy_cent_dist,
    N_MC, rand_01_MC, cos_t, sin_t, KP_Bys_rc, KP_Bys_rt, KP_Bys_ecc,
    KP_Bys_theta, KP_plot, KP_conct_par, KP_memb_num, field_dens,
    field_dens_std, cont_index, cl_region_i, frac_cl_area, cl_region_rjct_i,
        field_regions_rjct_i, field_regions_i, flag_no_fl_regs_i, **kwargs):
    """
    Make A3 block plots.
    """
    if 'A3' in pd['flag_make_plot']:

        fig = plt.figure(figsize=(figsize_x, figsize_y))
        gs = gridspec.GridSpec(grid_y, grid_x)
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
        #rdp_radii, rdp_points, rdp_stddev = prep_plots.RDPCurve(
        #    xy_filtered, xy_cent_dist, kde_cent, N_MC, rand_01_MC, cos_t,
        #    sin_t)
        rdp_radii, rdp_points, rdp_stddev = prep_plots.RDPellipse(
            xy_filtered, xy_cent_dist, kde_cent, KP_Bys_ecc[3], KP_Bys_theta[3], 
            N_MC, rand_01_MC, cos_t, sin_t)



        # Structure plots.
        arglist = [
            # pl_rad_find
            [gs, pd['plot_style'], coord, clust_rad, e_rad, rad_rads,
             rad_N_membs, rad_N_field, rad_CI],
            # pl_rad_dens: Radial density plot.
            [gs, pd['plot_style'], coord, rdp_radii, rdp_points, rdp_stddev,
             field_dens, field_dens_std, clust_rad, e_rad, pd['kp_flag'],
             KP_Bys_rc, KP_Bys_rt, KP_plot, KP_conct_par],
            # pl_zoom_frame: Zoom on x,y finding chart.
            [gs, fig, x_name, y_name, coord, x_zmin, x_zmax, y_zmin, y_zmax,
             cont_index, x_data_z, y_data_z, st_sizes_arr_z,
             kde_cent, clust_rad, KP_Bys_rc, KP_Bys_rt, KP_Bys_ecc,
             KP_Bys_theta, frac_cl_area, pd['kp_flag']],
            # pl_dens_map
            [gs, fig, asp_ratio, x_name, y_name, coord, x_zmin, x_zmax,
             y_zmin, y_zmax, kde_cent, frame_kde_cent, clust_rad,
             pd['kp_flag'], KP_Bys_rc, KP_Bys_rt, KP_Bys_ecc, KP_Bys_theta],
            # pl_mag_membs
            [gs, pd['plot_style'], y_ax, membvsmag],
            # pl_cl_fl_regions: Cluster and field regions defined.
            [gs, fig, pd['plot_style'], x_name, y_name, coord, x_min, x_max,
             y_min, y_max, asp_ratio, kde_cent, clust_rad, field_regions_i,
             field_regions_rjct_i, cl_region_i, cl_region_rjct_i,
             flag_no_fl_regs_i]
        ]
        for n, args in enumerate(arglist):
            mp_radius.plot(n, *args)

        # Generate output file.
        # Ignore warning issued by log-log inset.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fig.tight_layout()
        plt.savefig(
            join(npd['output_subdir'], str(npd['clust_name']) + '_A3'),
            bbox_inches='tight')
        # Close to release memory.
        plt.clf()
        plt.close("all")

        # Bayesian convergence plots
        if pd['kp_flag']:

            fig = plt.figure(figsize=(figsize_x, figsize_y))
            gs = gridspec.GridSpec(grid_y, grid_x)
            add_version_plot.main(y_fix=.998)

            arglist = [
                # pl_KP_Bys
                [gs, coord, pd['kp_nburn'], KP_plot, KP_Bys_rc, KP_Bys_rt,
                 KP_Bys_ecc, KP_Bys_theta]
            ]
            for n, args in enumerate(arglist):
                mp_KP_bayes.plot(n, *args)

            fig.tight_layout()
            plt.savefig(
                join(npd['output_subdir'], str(npd['clust_name']) + '_A3_KP'),
                bbox_inches='tight')
            # Close to release memory.
            plt.clf()
            plt.close("all")

        print("<<Plots for A3 block created>>")
    else:
        print("<<Skip A3 block plot>>")

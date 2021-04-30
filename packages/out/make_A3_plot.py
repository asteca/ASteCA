
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from os.path import join
import warnings
from . import mp_radius, mp_KP_bayes
from . import add_version_plot
from . import prep_plots
from . prep_plots import figsize_x, figsize_y, grid_x, grid_y


def main(npd, cld_i, pd, clp):
    """
    Make A3 block plots.
    """
    fig = plt.figure(figsize=(figsize_x, figsize_y))
    gs = gridspec.GridSpec(grid_y, grid_x)
    add_version_plot.main(y_fix=.998)

    # Obtain plotting parameters and data.
    x_min, x_max, y_min, y_max = prep_plots.frame_max_min(
        cld_i['x'], cld_i['y'])
    asp_ratio = prep_plots.aspect_ratio(x_min, x_max, y_min, y_max)
    coord, x_name, y_name = prep_plots.coord_syst(pd['coords'])
    x_zmin, x_zmax, y_zmin, y_zmax = prep_plots.frame_zoomed(
        x_min, x_max, y_min, y_max, clp['kde_cent'], clp['clust_rad'],
        pd['kp_ndim'], clp['KP_Bys_rt'])
    x_data_z, y_data_z, mag_data_z = prep_plots.zoomed_frame(
        cld_i['x'], cld_i['y'], cld_i['mags'], x_zmin, x_zmax, y_zmin,
        y_zmax)
    st_sizes_arr_z = prep_plots.star_size(mag_data_z)
    _, y_ax = prep_plots.ax_names(pd['colors'][0], pd['filters'][0], 'mag')
    rdp_radii, rdp_points, rdp_stddev, rad_max = prep_plots.RDPCurve(
        pd['kp_ndim'], clp['xy_filtered'], clp['xy_cent_dist'],
        clp['kde_cent'], clp['clust_rad'], clp['KP_Bys_ecc'][3],
        clp['KP_Bys_theta'][3])
    membvsmag = prep_plots.NmembVsMag(
        cld_i['x'], cld_i['y'], cld_i['mags'], clp['kde_cent'],
        clp['clust_rad'], clp['cl_area'])
    CI_vals, N_membs, N_membs_16, N_membs_84, membs_ratio,\
        membs_ratio_smooth = prep_plots.membVSrad(
            clp['field_dens'], clp['field_dens_std'], clp['rad_radii'],
            clp['rad_areas'], clp['N_in_cl_rad'], clp['N_in_ring'])

    # Structure plots.
    arglist = [
        # pl_rad_find
        [gs, pd['plot_style'], coord, clp['clust_rad'], membs_ratio,
         membs_ratio_smooth, CI_vals, clp['rad_radii']],
        # pl_mag_membs
        [gs, pd['plot_style'], y_ax, membvsmag],
        # pl_cl_fl_regions: Cluster and field regions defined.
        [gs, fig, pd['plot_style'], x_name, y_name, coord, x_min, x_max,
         y_min, y_max, asp_ratio, clp['kde_cent'], clp['clust_rad'],
         clp['field_regions_i'], clp['field_regions_rjct_i'],
         clp['cl_region_i'], clp['cl_region_rjct_i'],
         clp['flag_no_fl_regs_i']],
        # pl_rad_dens: Radial density plot.
        [gs, pd['plot_style'], coord, rdp_radii, rdp_points, rdp_stddev,
         rad_max, clp['field_dens'], clp['field_dens_std'], clp['clust_rad'],
         pd['kp_ndim'], clp['KP_Bys_rc'], clp['KP_Bys_rt'],
         clp['KP_plot'], clp['KP_conct_par']],
        # pl_zoom_frame: Zoom on x,y finding chart.
        [gs, fig, x_name, y_name, coord, x_zmin, x_zmax, y_zmin, y_zmax,
         clp['cont_index'], x_data_z, y_data_z, st_sizes_arr_z,
         clp['kde_cent'], clp['clust_rad'], clp['KP_Bys_rc'], clp['KP_Bys_rt'],
         clp['KP_Bys_ecc'], clp['KP_Bys_theta'], clp['frac_cl_area'],
         pd['kp_ndim']],
        # pl_memb_vs_rad
        [gs, pd['plot_style'], coord, clp['clust_rad'], clp['rad_radii'],
         N_membs, N_membs_16, N_membs_84, clp['KP_Bys_rc'][1],
         clp['KP_Bys_rt'][1], pd['kp_ndim'], clp['KP_plot']],
        # pl_membs_dist
        [gs, fig, clp['members_dist'], clp['n_memb_i']]
    ]
    for n, args in enumerate(arglist):
        mp_radius.plot(n, *args)

    # Generate output file.
    # Ignore warning issued by log-log inset.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fig.tight_layout()
    plt.savefig(
        join(npd['output_subdir'], str(npd['clust_name']) + '_A3'
             + npd['ext']))
    # Close to release memory.
    plt.clf()
    plt.close("all")

    # Bayesian convergence plots
    if pd['kp_ndim'] in (2, 4):

        fig = plt.figure(figsize=(figsize_x, figsize_y))
        gs = gridspec.GridSpec(grid_y, grid_x)
        add_version_plot.main(y_fix=.998)

        arglist = [
            # pl_KP_Bys
            [gs, coord, pd['kp_nburn'], clp['KP_plot'], clp['KP_Bys_rc'],
             clp['KP_Bys_rt'], clp['KP_Bys_ecc'], clp['KP_Bys_theta']]
        ]
        for n, args in enumerate(arglist):
            mp_KP_bayes.plot(n, *args)

        fig.tight_layout()
        plt.savefig(
            join(npd['output_subdir'], str(npd['clust_name']) + '_A3_KP'
                 + npd['ext']))
        # Close to release memory.
        plt.clf()
        plt.close("all")

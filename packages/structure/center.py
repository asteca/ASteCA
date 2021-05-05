
import numpy as np
from scipy import stats
from .xy_density import cent_bin as center_bin
from ..inp.get_data import coordsProject


def main(cld_i, clp, project, cent_method, **kwargs):
    """
    Obtains the center of the putative cluster.
    """

    # Restrict the KDE to a smaller area (to improve performance).
    radius = 0.25 * min(
        np.nanmax(cld_i['x']) - np.nanmin(cld_i['x']),
        np.nanmax(cld_i['y']) - np.nanmin(cld_i['y']))

    if cent_method[0] == 'a':
        # Obtain center coordinates as the maximum KDE value. Use the
        # approximate center obtained with the full frame and the given
        # bandwidth.
        cent = clp['kde_approx_cent']
        kde_cent, kde_plot = kde_center_zoom(
            cld_i['x'], cld_i['y'], cent, radius)

        # Find bin where the center xy coordinates are located.
        bin_cent = center_bin(clp['xedges'], clp['yedges'], kde_cent)

        print("Auto center found (bw={:g}): ({:g}, {:g}) deg".format(
            clp['bw_list'][1], kde_cent[0], kde_cent[1]))

    else:
        # De-project center coordinates if needed.
        x0, y0, _, _ = coordsProject(
            cent_method[0], cent_method[1], project,
            clp['x_offset'], clp['y_offset'])

        # Obtain KDE plot.
        _, kde_plot = kde_center_zoom(cld_i['x'], cld_i['y'], (x0, y0), radius)

        kde_cent = (x0, y0)
        print("Manual center fixed: ({:g}, {:g}) deg".format(*cent_method))

        # Find bin where the center xy coordinates are located.
        bin_cent = center_bin(clp['xedges'], clp['yedges'], kde_cent)

    # Add data to dictionary.
    center_params = {
        'kde_cent': kde_cent, 'kde_plot': kde_plot, 'bin_cent': bin_cent}
    clp.update(center_params)

    return clp


def kde_center_zoom(x_data, y_data, kde_approx_cent, radius):
    """
    Find the KDE maximum value pointing to the center coordinates.
    """

    # Unpack approximate center values.
    x_cent_pix, y_cent_pix = kde_approx_cent

    # Generate zoom around approx center value to speed things up.
    xmin_z, xmax_z = x_cent_pix - radius, x_cent_pix + radius
    ymin_z, ymax_z = y_cent_pix - radius, y_cent_pix + radius
    ext_range = [xmin_z, xmax_z, ymin_z, ymax_z]
    # Use reduced region around the center.
    x_zoom, y_zoom = [], []
    for indx, star_x in enumerate(x_data):
        if xmin_z < star_x < xmax_z and ymin_z < y_data[indx] < ymax_z:
            x_zoom.append(star_x)
            y_zoom.append(y_data[indx])
    values = np.vstack([x_zoom, y_zoom])

    # Check if there is at least one star selected here.
    if not values.any():
        raise ValueError(
            "ERROR: cluster region is empty and no center value\n"
            "could be estimated. Check that x,y columns are correct\n"
            "in 'params_input.dat' file.")

    # Obtain Gaussian KDE.
    kernel = stats.gaussian_kde(values)
    # Grid density (number of points).
    gd = 100
    gd_c = complex(0, gd)
    # Define x,y grid.
    x_grid, y_grid = np.mgrid[xmin_z:xmax_z:gd_c, ymin_z:ymax_z:gd_c]
    positions = np.vstack([x_grid.ravel(), y_grid.ravel()])
    # Evaluate kernel in grid positions.
    k_pos = kernel(positions)
    # Coordinates of max value in x,y grid (ie: center position).
    kde_cent = positions.T[np.argmax(k_pos)]

    # Pass for plotting.
    kde_plot = [ext_range, x_grid, y_grid, k_pos]

    return kde_cent, kde_plot

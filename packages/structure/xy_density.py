
import numpy as np
# from scipy.ndimage.filters import gaussian_filter
from scipy import stats
import bisect
from ..update_progress import updt


def main(clp, cld_i, center_bw, flag_make_plot, **kwargs):
    """
    Obtain Gaussian filtered 2D x,y histograms and the maximum values in them
    as centers.
    """
    print("Obtaining KDEs for the frame's coordinates.")

    # Filter possible nan values in (x, y)
    mskx, msky = np.isnan(cld_i['x']), np.isnan(cld_i['y'])
    msk = ~mskx & ~msky
    x_data, y_data = cld_i['x'][msk], cld_i['y'][msk]

    if center_bw == 0.:
        # Use Scotts factor (scipy's default).
        values = np.vstack([x_data, y_data])
        kernel = stats.gaussian_kde(values)
        c_bw = kernel.covariance_factor() * np.max(values.std(axis=1))
    else:
        c_bw = center_bw

    # KDE factor values for the KDE filter.
    bw_list = (c_bw * .5, c_bw, c_bw * 2.)

    cents_xy, frame_kdes = [], []
    if 'A1' in flag_make_plot:
        N = 1.
        for xym_rang in clp['xy_mag_ranges']:
            xr, yr, dummy = list(zip(*list(xym_rang.values())[0]))
            for bdw in bw_list:
                a, b = kde_center(xr, yr, bdw)
                cents_xy.append(a)
                frame_kdes.append(b)
                updt(15., N)
                N += 1.
        kde_approx_cent, frame_kde_cent = cents_xy[1], frame_kdes[1]
    else:
        # Only obtain the KDE for the full magnitude range.
        kde_approx_cent, frame_kde_cent = kde_center(
            x_data, y_data, bw_list[1])

    # Run once more for plotting.
    kernel, x_grid, y_grid, positions, k_pos = kde_center(
        x_data, y_data, bw_list[1], True)
    kde_dens_max, kde_dens_min = coordsDens(
        len(x_data), x_grid, y_grid, kernel, positions, k_pos)

    clp['cents_xy'], clp['kde_approx_cent'], clp['bw_list'],\
        clp['frame_kdes'], clp['frame_kde_cent'], clp['kde_dens_max'],\
        clp['kde_dens_min'] = cents_xy, kde_approx_cent, bw_list, frame_kdes,\
        frame_kde_cent, kde_dens_max, kde_dens_min

    return clp


def kde_center(x_data, y_data, bdw, plotFlag=False):
    '''
    Find the KDE maximum value pointing to the center coordinates.
    '''
    values = np.vstack([x_data, y_data])
    # Obtain Gaussian KDE.
    kernel = stats.gaussian_kde(
        values, bw_method=bdw / np.max(values.std(axis=1)))
    # Grid density (number of points).
    gd = 50
    gd_c = complex(0, gd)
    # Define x,y grid.
    xmin, xmax = np.min(x_data), np.max(x_data)
    ymin, ymax = np.min(y_data), np.max(y_data)
    x_grid, y_grid = np.mgrid[xmin:xmax:gd_c, ymin:ymax:gd_c]
    positions = np.vstack([x_grid.ravel(), y_grid.ravel()])
    # Evaluate kernel in grid positions.
    k_pos = kernel(positions)
    if plotFlag:
        return kernel, x_grid, y_grid, positions, k_pos

    # Coordinates of max value in x,y grid (ie: center position).
    kde_cent = positions.T[np.argmax(k_pos)]

    # Pass for plotting.
    ext_range = [xmin, xmax, ymin, ymax]
    kde_plot = [ext_range, x_grid, y_grid, k_pos]

    return kde_cent, kde_plot


def cent_bin(xedges, yedges, xy_cent):
    '''
    Take x,y center coordinates and return the 2D bin where they are located.
    '''
    x_cent_bin = bisect.bisect_left(xedges, xy_cent[0])
    y_cent_bin = bisect.bisect_left(yedges, xy_cent[1])
    # Store center bin coords for the filtered hist.
    cent_bin = [(x_cent_bin - 1), (y_cent_bin - 1)]

    return cent_bin


def coordsDens(N_stars, x_grid, y_grid, kernel, positions, k_pos):
    """
    Values used fort plotting the colorbar in the coordinates density map
    of the A2 block.
    """
    # Use a "bin width" (an area) that is dependent on the coordinates used
    # (pixels or celestials), but that it is small enough to converge to the
    # actual density value.
    bw = np.mean([
        x_grid[:, 0][1] - x_grid[:, 0][0], y_grid[0, :][1] - y_grid[0, :][0]])
    bw = bw / 10.

    # Coordinates for maximum KDE value
    kde_cent = positions.T[np.argmax(k_pos)]
    intg = kernel.integrate_box((kde_cent - bw), (kde_cent + bw))
    kde_dens_max = intg * N_stars / bw**2
    # Coordinates for minimum KDE value
    kde_min = positions.T[np.argmin(k_pos)]
    intg = kernel.integrate_box((kde_min - bw), (kde_min + bw))
    kde_dens_min = intg * N_stars / bw**2

    return kde_dens_max, kde_dens_min

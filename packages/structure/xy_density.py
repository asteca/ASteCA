
import numpy as np
# from scipy.ndimage.filters import gaussian_filter
from scipy import stats
import bisect
from ..update_progress import updt


def main(
    clp, cld, flag_make_plot, mirror_flag=True, Nmax=25000, center_bw=0.,
        **kwargs):
    """
    Obtain Gaussian filtered 2D x,y histograms and the maximum values in them
    as centers.

    HARDCODED
    Nmax: maximum number of stars used in the center estimation process
    """
    print("Obtaining KDEs for the frame's coordinates")

    x_data, y_data = cld['x'], cld['y']
    if x_data.size > Nmax:
        print(("  WARNING: too many stars. Selecting {} random\n"
               "  stars for center estimation").format(Nmax))
        step = int(x_data.size / Nmax)
        x_data, y_data = x_data[::step], y_data[::step]

    if center_bw == 0.:
        # Use half of Scotts factor (scipy's default).
        values = np.vstack([x_data, y_data])
        kernel = stats.gaussian_kde(values)
        c_bw = kernel.covariance_factor() * np.max(values.std(axis=1)) * .5
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
                a, b = kde_center(xr, yr, bdw, mirror_flag)
                cents_xy.append(a)
                frame_kdes.append(b)
                updt(15., N)
                N += 1.
        kde_approx_cent, frame_kde_cent = cents_xy[1], frame_kdes[1]
    else:
        # Only obtain the KDE for the full magnitude range.
        kde_approx_cent, frame_kde_cent = kde_center(
            x_data, y_data, bw_list[1], mirror_flag)

    clp['cents_xy'], clp['kde_approx_cent'], clp['bw_list'],\
        clp['frame_kdes'], clp['frame_kde_cent'] = cents_xy, kde_approx_cent,\
        bw_list, frame_kdes, frame_kde_cent

    return clp


def kde_center(x_data, y_data, bdw, mirror_flag):
    """
    Find the KDE maximum value pointing to the center coordinates.
    """

    xmin, xmax = np.min(x_data), np.max(x_data)
    ymin, ymax = np.min(y_data), np.max(y_data)
    values = np.vstack([x_data, y_data])

    # Fix KDE artifact around the borders
    if mirror_flag:
        values = dataMirror(values.T)

    # Obtain Gaussian KDE.
    try:
        kernel = stats.gaussian_kde(
            values, bw_method=bdw / np.max(values.std(axis=1)))
    except ValueError:
        print("  WARNING: could not generate coordinates KDE")
        return [np.nan, np.nan], []

    # Grid density (number of points).
    gd = 50
    gd_c = complex(0, gd)
    # Define x,y grid.
    x_grid, y_grid = np.mgrid[xmin:xmax:gd_c, ymin:ymax:gd_c]
    positions = np.vstack([x_grid.ravel(), y_grid.ravel()])
    # Evaluate kernel in grid positions.

    try:
        k_pos = kernel(positions)
    except np.linalg.LinAlgError:
        print("  WARNING: could not generate coordinates KDE")
        return [np.nan, np.nan], []

    # Coordinates of max value in x,y grid (ie: center position).
    kde_cent = positions.T[np.argmax(k_pos)]

    # Pass for plotting.
    ext_range = [xmin, xmax, ymin, ymax]
    kde_plot = [ext_range, x_grid, y_grid, k_pos]

    return kde_cent, kde_plot


def dataMirror(data, perc=.05):
    """
    Mirror a small percentage of data in all borders to remove the KDE artifact
    that lowers the density close to the edges.

    Source: https://stackoverflow.com/a/33602171/1391441
    """

    # Box
    xmin, ymin = np.min(data, 0)
    xmax, ymax = np.max(data, 0)
    midy = (ymax - ymin) * .5

    points_left_up = np.copy(data)
    points_left_up[:, 0] = xmin - (points_left_up[:, 0] - xmin)
    points_left_up[:, 1] += midy

    points_left_down = np.copy(data)
    points_left_down[:, 0] = xmin - (points_left_down[:, 0] - xmin)
    points_left_down[:, 1] -= midy

    points_right_up = np.copy(data)
    points_right_up[:, 0] = xmax + (xmax - points_right_up[:, 0])
    points_right_up[:, 1] += midy

    points_right_down = np.copy(data)
    points_right_down[:, 0] = xmax + (xmax - points_right_down[:, 0])
    points_right_down[:, 1] -= midy

    points_down = np.copy(data)
    points_down[:, 1] = ymin - (points_down[:, 1] - ymin)
    points_up = np.copy(data)
    points_up[:, 1] = ymax + (ymax - points_up[:, 1])

    # plt.scatter(*data.T, c='g', alpha=.3)
    # plt.scatter(*points_left_up.T, c='b', alpha=.2)
    # plt.scatter(*points_left_down.T, c='r', alpha=.2)
    # plt.scatter(*points_right_up.T, c='b', alpha=.2)
    # plt.scatter(*points_right_down.T, c='r', alpha=.2)
    # plt.scatter(*points_up.T, c='cyan', alpha=.2)
    # plt.scatter(*points_down.T, c='orange', alpha=.2)
    # plt.show()

    # Mirrored data
    mirror_data = np.append(
        np.append(points_left_down, points_left_up, axis=0),
        np.append(points_right_down, points_right_up, axis=0), axis=0)
    mirror_data = np.append(
        mirror_data, np.append(points_down, points_up, axis=0), axis=0)

    # Combine all data
    points = np.append(data, mirror_data, axis=0)

    # Trim mirrored frame to within a 'perc' pad of the full (x, y) range
    xr, yr = np.ptp(data.T[0]) * perc, np.ptp(data.T[1]) * perc
    xmin, xmax = xmin - xr, xmax + xr
    ymin, ymax = ymin - yr, ymax + yr
    msk = (points[:, 0] > xmin) & (points[:, 0] < xmax) &\
        (points[:, 1] > ymin) & (points[:, 1] < ymax)

    return points[msk].T


def cent_bin(xedges, yedges, xy_cent):
    """
    Take x,y center coordinates and return the 2D bin where they are located.
    """
    x_cent_bin = bisect.bisect_left(xedges, xy_cent[0])
    y_cent_bin = bisect.bisect_left(yedges, xy_cent[1])
    # Store center bin coords for the filtered hist.
    cent_bin = [(x_cent_bin - 1), (y_cent_bin - 1)]

    return cent_bin

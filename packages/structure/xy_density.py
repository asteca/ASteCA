
import numpy as np
# from scipy.ndimage.filters import gaussian_filter
from scipy import stats
import bisect
from ..update_progress import updt


def main(clp, cld_i, center_kf, flag_make_plot, **kwargs):
    """
    Obtain Gaussian filtered 2D x,y histograms and the maximum values in them
    as centers.
    """
    print("Obtaining KDEs for the frame's coordinates.")

    if center_kf == 0.:
        # Scotts factor (scipy's default).
        values = np.vstack([cld_i['x'], cld_i['y']])
        kernel = stats.gaussian_kde(values)
        kf = kernel.covariance_factor()
    else:
        kf = center_kf

    # KDE factor values for the KDE filter.
    kf_list = (kf * .5, kf, kf * 2.)

    cents_xy, frame_kdes = [], []
    if 'A1' in flag_make_plot:
        N = 1.
        for xym_rang in clp['xy_mag_ranges']:
            xr, yr, dummy = list(zip(*list(xym_rang.values())[0]))
            for k_f in kf_list:
                a, b = kde_center(xr, yr, k_f)
                cents_xy.append(a)
                frame_kdes.append(b)
                updt(15., N)
                N += 1.
        kde_approx_cent, frame_kde_cent = cents_xy[1], frame_kdes[1]
    else:
        # Only obtain the KDE for the full magnitude range.
        kde_approx_cent, frame_kde_cent = kde_center(
            cld_i['x'], cld_i['y'], kf_list[1])

    clp['cents_xy'], clp['kde_approx_cent'], clp['kf_list'],\
        clp['frame_kdes'], clp['frame_kde_cent'] = cents_xy, kde_approx_cent,\
        kf_list, frame_kdes, frame_kde_cent

    return clp


def kde_center(x_data, y_data, k_f):
    '''
    Find the KDE maximum value pointing to the center coordinates.
    '''
    values = np.vstack([x_data, y_data])
    # Obtain Gaussian KDE.
    kernel = stats.gaussian_kde(values, bw_method=k_f)
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

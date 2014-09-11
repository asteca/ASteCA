
import numpy as np
import bisect
import matplotlib.pyplot as plt
from scipy import stats
from scipy.ndimage.filters import gaussian_filter
from display_cent import disp_cent as d_c


def center_approx(hist_lst, st_dev_lst):
    '''
    Function that returns a filtered 2D histogram and approximate center
    coordinates obtained using different standard deviations.
    '''

    hist, xedges, yedges = hist_lst

    approx_cents = []
    for indx, st_dev in enumerate(st_dev_lst):

        # 2D histogram with a gaussian filter applied.
        h_g = gaussian_filter(hist, st_dev, mode='constant')

        # x,y coordinates of the bin with the maximum value in the 2D
        # filtered histogram.
        x_cent_bin, y_cent_bin = np.unravel_index(h_g.argmax(), h_g.shape)
        # x,y coords of the center of the above bin.
        x_cent_pix, y_cent_pix = np.average(xedges[x_cent_bin:x_cent_bin + 2]),\
        np.average(yedges[y_cent_bin:y_cent_bin + 2])

        # Only pass this one for plotting.
        if indx == 0:
            hist_2d_g = h_g

        approx_cents.append([x_cent_pix, y_cent_pix])

    # Find bin where the center xy coordinates are located.
    x_cent_bin = bisect.bisect_left(xedges, approx_cents[0][0])
    y_cent_bin = bisect.bisect_left(yedges, approx_cents[0][1])
    # Store center bin coords for the filtered hist.
    cent_bin = [(x_cent_bin - 1), (y_cent_bin - 1)]

    return hist_2d_g, cent_bin, approx_cents


def kde_center(x_data, y_data, approx_cents, radius, gc_params):
    '''
    Find the KDE maximum value wich points to the center coordinates.
    '''

    kde_centers = []
    for indx_c, center in enumerate(approx_cents):

        # Unpack approximate center values.
        x_cent_pix, y_cent_pix = center

        # Generate zoom around approx center value to spped things up.
        xmin_z, xmax_z = x_cent_pix - radius, x_cent_pix + radius
        ymin_z, ymax_z = y_cent_pix - radius, y_cent_pix + radius
        # Use reduced region around the center.
        x_zoom, y_zoom = [], []
        for indx, star_x in enumerate(x_data):
            if xmin_z < star_x < xmax_z and ymin_z < y_data[indx] < ymax_z:
                x_zoom.append(star_x)
                y_zoom.append(y_data[indx])
        values = np.vstack([x_zoom, y_zoom])

        # Obtain Gaussian KDE.
        kernel = stats.gaussian_kde(values)
        # Define bandwidth value.
        if gc_params[0] != 'auto':
            # Change badwidth to manual value.
            kernel.set_bandwidth(bw_method=gc_params[1])

        # Define x,y grid.
        # Grid density (number of points).
        gd = 100
        gd_c = complex(0, gd)
        x, y = np.mgrid[xmin_z:xmax_z:gd_c, ymin_z:ymax_z:gd_c]
        positions = np.vstack([x.ravel(), y.ravel()])

        # Evaluate kernel in grid positions.
        k_pos = kernel(positions)
        # Usa values obtained with the approx center derived via the minimum
        # standard deviation value applied to the 2D histogram.
        if indx_c == 0:
            ext_range = [xmin_z, xmax_z, ymin_z, ymax_z]
            x_grid, y_grid = x, y
            k_pos_plot = k_pos
            # The error is associated with the grid density used and the
            # zoomed area defined.
            x_range, y_range = max(x_zoom) - min(x_zoom), max(y_zoom) - \
            min(y_zoom)
            e_cent = [x_range / gd, y_range / gd]
        # Coordinates of max value in x,y grid (ie: center position).
        x_cent_kde, y_cent_kde = positions.T[np.argmax(k_pos)]
        # Append values to list.
        kde_centers.append([x_cent_kde, y_cent_kde])

    # Pass for plotting.
    kde_plot = [ext_range, x_grid, y_grid, k_pos_plot]

    return kde_centers, e_cent, kde_plot


def get_center(x_data, y_data, mag_data, hist_lst, gc_params, mode,
    semi_return, coord_lst):
    """
    Obtains the center of the putative cluster. Returns the center values
    along with its errors and several arrays related to histograms, mainly for
    plotting purposes.
    """

    coord = coord_lst[0]

    st_dev_lst = [2., 2.5, 3., 3.5, 4.]
    # Set flags.
    flag_center_med, flag_center_std = False, False
    flag_center_manual = False

    mode_semi = True
    if mode == 'semi':
        # Unpack semi values.
        cent_cl_semi, cl_rad_semi, cent_flag_semi = semi_return[:3]

        # Only apply if flag is true, else skip semi-center assignment for
        # this cluster.
        if cent_flag_semi == 1:
            # Search for new center values using the center coordinates
            # and radius given as initial values.

            # Call funct to obtain the pixel coords of the maximum KDE value.
            approx_cent_semi = [cent_cl_semi]
            kde_centers, e_cent, kde_plot = kde_center(x_data, y_data,
                approx_cent_semi, cl_rad_semi, gc_params)

            # Find bin where the center xy coordinates are located.
            hist, xedges, yedges = hist_lst[:-1]
            x_cent_bin = bisect.bisect_left(xedges, kde_centers[0][0])
            y_cent_bin = bisect.bisect_left(yedges, kde_centers[0][1])
            # Store center bin coords for the filtered hist.
            cent_bin = [(x_cent_bin - 1), (y_cent_bin - 1)]

            # For plotting.
            # 2D histogram with a gaussian filter applied.
            hist_2d_g = gaussian_filter(hist, st_dev_lst[0], mode='constant')

            print 'Semi center found: ({:g}, {:g}) {c}.'.format(
                *kde_centers[0], c=coord)
        else:
            # Use 'auto' mode.
            mode_semi = False

    if mode == 'auto' or mode_semi is False:

        # Obtain approximate values for center coordinates using several
        # Gaussian filters with different standard deviation values, on the
        # 2D histogram.
        hist_2d_g, cent_bin, approx_cents = center_approx(hist_lst[:-1],
            st_dev_lst)

        # This is the radius used in auto mode to restrict the search of the
        # KDE center coordinates.
        x_span, y_span = max(x_data) - min(x_data), max(y_data) - min(y_data)
        radius = 0.25 * min(x_span, y_span)

        # Call funct to obtain the pixel coords of the maximum KDE value.
        kde_centers, e_cent, kde_plot = kde_center(x_data, y_data, approx_cents,
            radius, gc_params)

        # Calculate the median value for the cluster's center (use median
        # instead of mean to reject possible outliers) and the standard
        # deviation using all the coordinates obtained.
        cent_median, cent_std_dev = np.mean(np.array(kde_centers), axis=0), \
        np.std(np.array(kde_centers), axis=0)

        # Raise a flag if either median cluster's central coordinate is
        # more than 10% away from the ones assigned as the cluster's center.
        if abs(cent_median[0] - kde_centers[0][0]) > 0.1 * kde_centers[0][0] \
        or abs(cent_median[1] - kde_centers[0][1]) > 0.1 * kde_centers[0][1]:
            flag_center_med = True
        # Raise a flag if the standard deviation for either coord is larger than
        # 10% of the center coord values.
        if cent_std_dev[0] > 0.1 * kde_centers[0][0] or \
            cent_std_dev[1] > 0.1 * kde_centers[0][1]:
            flag_center_std = True

        print 'Auto center found: ({:g}, {:g}) {c}.'.format(*kde_centers[0],
            c=coord)

    # If Manual mode is set, display center and ask the user to accept it or
    # input new one.
    elif mode == 'manual':

        # Obtain approximate values for center coordinates using several
        # Gaussian filters with different standard deviation values, on the
        # 2D histogram.
        hist_2d_g, cent_bin, approx_cents = center_approx(hist_lst[:-1],
            [st_dev_lst[0]])

        # Show plot with center obtained.
        d_c(x_data, y_data, mag_data, approx_cents[0], cent_bin, hist_2d_g)
        plt.show()

        kde_centers = [approx_cents[0]]
        e_cent = [0., 0.]
        kde_plot = []

        wrong_answer = True
        while wrong_answer:
            answer_cen = raw_input('Input new center values? (y/n) ')
            if answer_cen == 'n':
                print 'Value accepted.'
                wrong_answer = False
            elif answer_cen == 'y':
                print 'Input new center values.'
                kde_centers[0][0] = float(raw_input('x: '))
                kde_centers[0][1] = float(raw_input('y: '))

                # Find bin where the center xy coordinates are located.
                hist, xedges, yedges = hist_lst[:-1]
                x_cent_bin = bisect.bisect_left(xedges, kde_centers[0][0])
                y_cent_bin = bisect.bisect_left(yedges, kde_centers[0][1])
                # Store center bin coords for the filtered hist.
                cent_bin = [(x_cent_bin - 1), (y_cent_bin - 1)]
                wrong_answer = False
                flag_center_manual = True
            else:
                print 'Wrong input. Try again.\n'

    center_params = [cent_bin, kde_centers, e_cent, st_dev_lst, hist_2d_g,
        kde_plot, flag_center_med, flag_center_std, flag_center_manual]

    return center_params
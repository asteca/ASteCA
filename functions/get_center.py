
import numpy as np
import bisect
import matplotlib.pyplot as plt
from scipy import stats
from display_cent import disp_cent as d_c


def kde_center(indx_b, x_data, y_data, x_cent_pix, y_cent_pix, radius,
    gc_params):
    '''
    Find the KDE maximum value wich points to the center coordinates.
    '''
    # Generate zoom around initial center value.
    xmin_z, xmax_z = x_cent_pix - radius, x_cent_pix + radius
    ymin_z, ymax_z = y_cent_pix - radius, y_cent_pix + radius
    # Use reduced region around the center.
    x_zoom, y_zoom = [], []
    for indx, star_x in enumerate(x_data):
        if xmin_z < star_x < xmax_z and ymin_z < y_data[indx] < ymax_z:
            x_zoom.append(star_x)
            y_zoom.append(y_data[indx])
    values = np.vstack([x_zoom, y_zoom])
    # Define x,y grid.
    x, y = np.mgrid[xmin_z:xmax_z:100j, ymin_z:ymax_z:100j]
    positions = np.vstack([x.ravel(), y.ravel()])
    # Obtain Gaussian KDE.
    kernel = stats.gaussian_kde(values)
    if gc_params[0] == 'auto':
        # Get default bandwidth value.
        bw = kernel.covariance_factor()
    else:
        # Use manual value.
        bw = gc_params[1]
    kde_centers = []
    indx_r = 4 if indx_b == 0 else 1
    for i in range(indx_r):
        kernel = stats.gaussian_kde(values, bw_method=(bw * (i + 1)))
        # Evaluate kernel in this positions.
        k_pos = kernel(positions)
        # Save this one for plotting.
        if i == 0:
            k_pos_plot = k_pos
        # Coordinates of max value in x,y grid.
        x_cent_kde, y_cent_kde = positions.T[np.argmax(k_pos)]
        # Append values to list.
        kde_centers.append([x_cent_kde, y_cent_kde])
    if indx_b == 0:
        x_cent_kde, y_cent_kde = np.mean(kde_centers, axis=0)
        e_cent = np.std(kde_centers, axis=0)
    else:
        e_cent = 0.

    # Pass for plotting.
    ext_range = [xmin_z, xmax_z, ymin_z, ymax_z]
    kde_plot = [ext_range, x, y, k_pos_plot]

    return x_cent_kde, y_cent_kde, e_cent, kde_plot


def get_center(x_data, y_data, mag_data, bin_list, gc_params, mode,
    semi_return):
    """
    Obtains the center of the putative cluster. Returns the center values
    along with its errors and several arrays related to histograms, mainly for
    plotting purposes.
    """

    # This is the radius used in auto mode to restrict the search of the
    # KDE center coordinates.
    x_span, y_span = max(x_data) - min(x_data), max(y_data) - min(y_data)
    radius = 0.15 * min(x_span, y_span)

    # Arrays that store the cluster's center values calculated varying
    # the bin size 'd_b'.
    centers_kde = []

    # Iterate for the defined bin widths.
    for indx, d_b in enumerate(bin_list):

        hist, xedges, yedges, h_g = center_fun(x_data, y_data, d_b)
        # x,y coordinates of the bin with the maximum value in the 2D
        # filtered histogram.
        x_cent_bin, y_cent_bin = np.unravel_index(h_g.argmax(), h_g.shape)
        # x,y coords of the center of the above bin.
        x_cent_pix, y_cent_pix = np.average(xedges[x_cent_bin:x_cent_bin + 2]),\
        np.average(yedges[y_cent_bin:y_cent_bin + 2])

        # Call funct to obtain the pixel coords of the maximum KDE value.
        x_cent_kde, y_cent_kde, e_cent, kde_plot = kde_center(indx, x_data,
            y_data, x_cent_pix, y_cent_pix, radius, gc_params)

        # Store center coords in pixel coordinates.
        centers_kde.append([x_cent_kde, y_cent_kde])

        # Only store for the smallest value of 'd_b'.
        if indx == 0:
            # Store center error.
            e_center = round(max(e_cent), 1)
            # Store min width bin edges.
            hist_xyedges = [xedges, yedges]
            # Find bin where the center xy coordinates are located.
            x_cent_bin = bisect.bisect_left(xedges, x_cent_kde)
            y_cent_bin = bisect.bisect_left(yedges, y_cent_kde)
            # Store center bin coords for the filtered hist.
            bin_center = [(x_cent_bin - 1), (y_cent_bin - 1)]
            # Store not-filtered 2D hist.
            h_not_filt = hist
            # Store filtered 2D hist.
            h_filter = h_g
            # Store KDE for plotting.
            kde_pl = kde_plot

    # Store auto found center with min bin width.
    center_cl = [centers_kde[0][0], centers_kde[0][1]]

    # Store in this array the values for the cluster's center x,y coordinates
    # obtained with the different bin sizes.
    arr_g = np.array(centers_kde)
    # Calculate the median value for the cluster's center (use median instead
    # of mean to reject possible outliers) and the standard deviation using all
    # the coordinates obtained with different bin widths.
    median_coords, std_dev = np.median(arr_g, axis=0), np.std(arr_g, axis=0)
    # Store stats values.
    cent_stats = [e_center, median_coords, std_dev]
    # Set flags.
    flag_center_med, flag_center_std = False, False
    # Raise a flag if either median cluster's central coordinate is more than
    # 10% away from the ones obtained with the min bin width.
    if abs(median_coords[0] - center_cl[0]) > 0.1 * center_cl[0] or \
    abs(median_coords[1] - center_cl[1]) > 0.1 * center_cl[1]:
        flag_center_med = True
    # Raise a flag if the standard deviation for either coord is larger than
    # 10% of the center coord values.
    if std_dev[0] > 0.1 * center_cl[0] or std_dev[1] > 0.1 * center_cl[1]:
        flag_center_std = True

    flag_center_manual = False
    if mode == 'auto':
        print 'Auto center found: (%0.2f, %0.2f) px.' % (center_cl[0],
        center_cl[1])

    elif mode == 'semi':
        # Unpack semi values.
        cent_cl_semi, cl_rad_semi, cent_flag_semi, rad_flag_semi, \
        err_flag_semi = semi_return

        # Only apply if flag is true, else skip semi-center assignment for
        # this cluster.
        if cent_flag_semi == 1:
            # Search for new center values using the center coordinates
            # and radius given as initial values.

            # Call funct to obtain the pixel coords of the maximum KDE value.
            x_cent_kde, y_cent_kde, e_cent, kde_plot = kde_center(0, x_data,
                y_data, cent_cl_semi[0], cent_cl_semi[1], cl_rad_semi,
                gc_params)

            # Update KDE for plotting.
            kde_pl = kde_plot

            # Update list with new center coords.
            centers_kde[0][0], centers_kde[0][1] = x_cent_kde, y_cent_kde

            # Find bin where the center xy coordinates are located.
            x_cent_bin = bisect.bisect_left(hist_xyedges[0], x_cent_kde)
            y_cent_bin = bisect.bisect_left(hist_xyedges[1], y_cent_kde)
            # Update center bin coords for the filtered hist.
            bin_center = [(x_cent_bin - 1), (y_cent_bin - 1)]

            print 'Semi center found: (%0.2f, %0.2f) px.' % (centers_kde[0][0],
                centers_kde[0][1])
        else:
            print 'Auto center found: (%0.2f, %0.2f) px.' % (center_cl[0],
            center_cl[1])

    # If Manual mode is set, display center and ask the user to accept it or
    # input new one.
    elif mode == 'manual':
        # Show plot with center obtained.
        d_c(x_data, y_data, mag_data, center_cl, bin_center, h_filter,
            bin_list[0])
        plt.show()

        wrong_answer = True
        while wrong_answer:
            answer_cen = raw_input('Input new center values? (y/n) ')
            if answer_cen == 'n':
                print 'Value accepted.'
                wrong_answer = False
            elif answer_cen == 'y':
                print 'Input new center values (in px).'
                center_cl[0] = float(raw_input('x: '))
                center_cl[1] = float(raw_input('y: '))
                # Store center values in bin coordinates. We substract
                # the min (x,y) coordinate values otherwise the bin
                # coordinates won't be aligned.
                bin_center = [int(round((center_cl[0] - min(x_data))
                 / bin_list[0])), int(round((center_cl[1] - min(y_data))
                  / bin_list[0]))]
                wrong_answer = False
                flag_center_manual = True
            else:
                print 'Wrong input. Try again.\n'

    center_params = [bin_list, h_not_filt, hist_xyedges, h_filter, bin_center,
        centers_kde, cent_stats, kde_pl, flag_center_med, flag_center_std,
        flag_center_manual]

    return center_params
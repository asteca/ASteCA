
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter
from display_cent import disp_cent as d_c


def center_fun(x_data, y_data, mag_data, d_b):
    '''
    Function that returns center coordinates for the data given.
    '''

    xmin, xmax = min(x_data), max(x_data)
    ymin, ymax = min(y_data), max(y_data)
    rang = [[xmin, xmax], [ymin, ymax]]
    # Obtain relative intensities.
    rel_int = 10 ** ((max(mag_data) - mag_data) / 2.5)

    # Number of bins in x,y given the bin width 'd_b'
    binsxy = [int((xmax - xmin) / d_b), int((ymax - ymin) / d_b)]

    # hist is the 2D histogran, *edges store the edges of the bins.
    hist, xedges, yedges = np.histogram2d(x_data, y_data,
        range=rang, bins=binsxy)
    hist_w, xedges_w, yedges_w = np.histogram2d(x_data, y_data,
        range=rang, bins=binsxy, weights=rel_int)

    # H_g is the 2D histogram with a gaussian filter applied.
    h_g = gaussian_filter(hist, 2, mode='constant')
    h_g_w = gaussian_filter(hist_w, 2, mode='constant')

    return hist, hist_w, h_g, h_g_w, xedges, yedges, xedges_w, yedges_w


def get_center(x_data, y_data, mag_data, gc_params, mode, semi_return):
    """
    Obtains the center of the putative cluster. Returns the center values
    along with its errors and several arrays related to histograms, mainly for
    plotting purposes.
    """

    xmin, xmax = min(x_data), max(x_data)
    ymin, ymax = min(y_data), max(y_data)
    rang = [[xmin, xmax], [ymin, ymax]]

    # Either calculate or read the number of bins used.
    if gc_params[0] == 'auto':
        min_rang = min((rang[0][1] - rang[0][0]), (rang[1][1] - rang[1][0]))
        # Number of bins given by 1%, 2% and 3% of the minimum axis range.
        bin_list = [int(i * min_rang / 100.) for i in range(1, 5)]
    else:
        bin_list = gc_params[1:]
        bin_list.sort()

    # Arrays that store the cluster's center values calculated varying
    # the bin size 'd_b'.
    centers_h, centers_hw = [], []

    # Iterate for the defined bin widths.
    for indx, d_b in enumerate(bin_list):

        hist, hist_w, h_g, h_g_w, xedges, yedges, xedges_w, yedges_w = \
        center_fun(x_data, y_data, mag_data, d_b)

        # x_cent_g & y_cent_g store the x,y coordinates of the bin with the
        # maximum value in the 2D filtered histogram.
        x_cent_g, y_cent_g = np.unravel_index(h_g.argmax(), h_g.shape)
        x_cent_g_w, y_cent_g_w = np.unravel_index(h_g_w.argmax(), h_g_w.shape)

        # Only store for the smallest value of 'd_b'.
        if indx == 0:
            # Store min width bin edges.
            xedges_min_db, yedges_min_db = xedges, yedges
            # Store not-filtered 2D hist.
            h_not_filt = hist
            # Store filtered 2D hist.
            h_filter = [h_g, h_g_w]
            # Store center bins for the filtered hist.
            x_center_bin, y_center_bin = x_cent_g, y_cent_g

        # Store center coords in pixel coordinates.
        centers_h.append([np.average(xedges[x_cent_g:x_cent_g + 2]),
            np.average(yedges[y_cent_g:y_cent_g + 2])])
        centers_hw.append([np.average(xedges_w[x_cent_g_w:x_cent_g_w + 2]),
            np.average(yedges_w[y_cent_g_w:y_cent_g_w + 2])])

    # Store in this array the values for the cluster's center x,y coordinates
    # obtained with the different bin sizes.
    arr_g = np.array(centers_h)
    # Calculate the median value for the cluster's center (use median instead
    # of mean to reject possible outliers) and the standard deviation using all
    # the coordinates obtained with different bin widths.
    median_coords, std_dev = np.median(arr_g, axis=0), np.std(arr_g, axis=0)
    # Set flags.
    flag_center_med, flag_center_std = False, False
    # Raise a flag if either median cluster's central coordinate is more than
    # 10% away from the ones obtained with the min bin width.
    if abs(median_coords[0] - centers_h[0][0]) > 0.1 * centers_h[0][0] or \
    abs(median_coords[1] - centers_h[0][1]) > 0.1 * centers_h[0][1]:
        flag_center_med = True
    # Raise a flag if the standard deviation for either coord is larger than
    # 10% the coord value.
    if std_dev[0] > 0.1 * centers_h[0][0] or \
    std_dev[1] > 0.1 * centers_h[0][1]:
        flag_center_std = True

    # Pass the center coordinates obtained with the smallest bin.
    center_coords = [centers_h, centers_hw]
    # Store auto found center with min bin width.
    center_cl = [center_coords[0][0][0], center_coords[0][0][1]]
    # Pass the errors as the width of the bin.
    cent_cl_err = [bin_list[0], bin_list[0]]

    flag_center_manual = False
    if mode == 'a':
        print 'Auto center found: (%0.2f, %0.2f) px.' % (center_cl[0],
        center_cl[1])

    elif mode == 's':
        # Unpack semi values.
        cent_cl_semi, cl_rad_semi, cent_flag_semi, rad_flag_semi, \
        err_flag_semi = semi_return

        if cent_flag_semi == 1:
            # Search for new center values using the initial center coordinates
            # and radius given.

            # Define maximum x,y values for the search.
            x_min, x_max = cent_cl_semi[0] - cl_rad_semi, \
            cent_cl_semi[0] + cl_rad_semi
            y_min, y_max = cent_cl_semi[1] - cl_rad_semi, \
            cent_cl_semi[1] + cl_rad_semi
            # Convert range in pixel coords to bin coords.
            x_bin_max, x_bin_min = int(x_max / bin_list[0]), \
            int(x_min / bin_list[0])
            y_bin_max, y_bin_min = int(y_max / bin_list[0]), \
            int(y_min / bin_list[0])

            # Create new array filled with inf.
            A = np.full_like(h_filter[0], -np.inf)
            # Slice new array.
            A[x_bin_min:x_bin_max, y_bin_min:y_bin_max] = \
            h_filter[0][x_bin_min:x_bin_max, y_bin_min:y_bin_max]
            # Search for maximum value in this region in bin coords.
            x_center_bin, y_center_bin = np.unravel_index(A.argmax(),
                A.shape)
            # Convert to pixel coordinates.
            B = [np.average(xedges_min_db[x_center_bin:x_center_bin + 2]),
                np.average(yedges_min_db[y_center_bin:y_center_bin + 2])]
            # Update list with new center coords.
            center_cl = B

            print 'Semi center found: (%0.2f, %0.2f) px.' % (center_cl[0],
                center_cl[1])

    # If Manual mode is set, display center and ask the user to accept it or
    # input new one.
    elif mode == 'm':
        # Show plot with center obtained.
        d_c(x_data, y_data, mag_data, center_cl, cent_cl_err, x_center_bin,
            y_center_bin, h_filter)
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
                # Update values.
                cent_cl_err[0], cent_cl_err[1] = 0., 0.
                # Store center values in bin coordinates. We substract
                # the min (x,y) coordinate values otherwise the bin
                # coordinates won't be aligned.
                x_center_bin, y_center_bin = int(round((center_cl[0] -
                min(x_data)) / bin_list[0])), int(round((center_cl[1] -
                min(y_data)) / bin_list[0]))
                wrong_answer = False
                flag_center_manual = True
            else:
                print 'Wrong input. Try again.\n'

    center_params = [center_cl, h_not_filt, x_center_bin,
        y_center_bin, xedges_min_db, yedges_min_db, bin_list[0], h_filter,
        flag_center_med, flag_center_std, flag_center_manual]

    return center_params
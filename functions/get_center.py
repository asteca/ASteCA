"""
@author: gabriel
"""

import numpy as np
from scipy.ndimage.filters import gaussian_filter


def get_center(x_data, y_data, mag_data, gc_params):
    """
    Obtains the center of the putative cluster. Returns the center values
    along with its errors and several arrays related to histograms, mainly for
    plotting purposes.
    """

    xmin, xmax = min(x_data), max(x_data)
    ymin, ymax = min(y_data), max(y_data)
    rang = [[xmin, xmax], [ymin, ymax]]

    # Arrays that store the cluster's center values calculated varying the bin
    # size 'd_b' and the gaussian filtered 2D histogram.
    center_x_g, center_y_g = [], []
    # Store the bin centers for each value of 'd_b'
    x_center_bin, y_center_bin = [], []
    # 2D histograms
    h_not_filt, h_filter = [], []

    # Either calculate or read the number of bins used.
    if gc_params[0] == 'auto':
        min_rang = min((rang[0][1] - rang[0][0]), (rang[1][1] - rang[1][0]))
        # Number of bins given by 1%, 2% and 3% of the minimum axis range.
        bin_list = [int(i * min_rang / 100.) for i in range(1, 4)]
    else:
        bin_list = gc_params[1:]
        bin_list.sort()

    # Iterate for the defined bin widths.
    for indx, d_b in enumerate(bin_list):
        # Number of bins in x,y given the bin width 'd_b'
        binsxy = [int((xmax - xmin) / d_b), int((ymax - ymin) / d_b)]

        # hist is the 2D histogran, xedges & yedges store the edges of the bins
        hist, xedges, yedges = np.histogram2d(x_data, y_data, range=rang,
                                              bins=binsxy)
        # Store not-filtered 2D hist in 'h_not_filt'.
        h_not_filt.append(hist)

        # Only store the edges for the smallest value of 'd_b'.
        if indx == 0:
            xedges_min_db, yedges_min_db = xedges, yedges

        # H_g is the 2D histogram with a gaussian filter applied
        h_g = gaussian_filter(hist, 2, mode='constant')
        # Store all filtered arrays with different bin sizes in 'h_filter'
        h_filter.append(h_g)

        # x_cent_g & y_cent_g store the x,y coordinates of the bin with the
        # maximum value in the 2D filtered histogram, ie: the center of the
        # putative cluster.
        # 'argmax' returns the flattened histogram's bin index that stores the
        # biggest value.
        # 'shape' returns the number of rows and columns in the histogram.
        # 'unravel' takes the 'argmax' value and using the 'shape' output
        # returns the bin index that stores the max value in the histogram in a
        # non-flattened array.
        x_cent_g, y_cent_g = np.unravel_index(h_g.argmax(), h_g.shape)
        # Store center bins for the filtered hist.
        x_center_bin.append(x_cent_g)
        y_center_bin.append(y_cent_g)

        # Calculate center coords in pixel coordinates.
        center_x_g.append(np.average(xedges[x_cent_g:x_cent_g + 2]))
        center_y_g.append(np.average(yedges[y_cent_g:y_cent_g + 2]))

    # Run again for the minimum bin width, this time weighting according to
    # relative intensities.
    # Obtain relative intensities.
    rel_int = 10 ** ((max(mag_data) - mag_data) / 2.5)
    # Use minimum bin width.
    d_b = bin_list[0]
    bin_list.append(d_b)
    binsxy = [int((xmax - xmin) / d_b), int((ymax - ymin) / d_b)]
    hist, xedges, yedges = np.histogram2d(x_data, y_data, range=rang,
                                          bins=binsxy, weights=rel_int)
    # Store not-filtered 2D hist.
    h_not_filt.append(hist)
    h_g = gaussian_filter(hist, 2, mode='constant')
    # Store filtered 2D hist.
    h_filter.append(h_g)
    x_cent_g, y_cent_g = np.unravel_index(h_g.argmax(), h_g.shape)
    # Store center bins for the filtered hist.
    x_center_bin.append(x_cent_g)
    y_center_bin.append(y_cent_g)
    # Append center coords in pixel coordinates.
    center_x_g.append(np.average(xedges[x_cent_g:x_cent_g + 2]))
    center_y_g.append(np.average(yedges[y_cent_g:y_cent_g + 2]))

    # Store in this array the values for the cluster's center x,y coordinates
    # obtained with the different bin sizes.
    arr_g = np.array([center_x_g, center_y_g])

    # Calculate the median value for the cluster's center (use median instead of
    # mean to reject possible outliers) and the standard deviation using all
    # the coordinates obtained with different bin widths.
    median_coords, std_dev = np.median(arr_g, axis=1), np.std(arr_g, axis=1)

    # Raise a flag if either median cluster's central coordinate is more than
    # 2 sigmas away from the ones obtained with the min bin width.
    flag_center = False
    if abs(median_coords[0] - center_x_g[0]) > 2 * std_dev[0] or \
    abs(median_coords[1] - center_y_g[0]) > 2 * std_dev[1]:
        flag_center = True

    # Pass the center coordinates obtained with the smallest bin.
    center_coords = [center_x_g, center_y_g]
    # Pass the errors as the width of the bin.
    cent_coo_err = bin_list

    return center_coords, cent_coo_err, h_filter, h_not_filt, xedges_min_db, \
    yedges_min_db, x_center_bin, y_center_bin, bin_list, flag_center

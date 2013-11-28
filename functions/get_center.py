"""
@author: gabriel
"""

import numpy as np
from scipy.ndimage.filters import gaussian_filter

def get_center(x_data, y_data):
    """   
    Obtain the center of the putative cluster and return a list of two
    items: [x_center, y_center] and the array that contains the gausssian
    filtered 2D histogram and lots of other stuff.
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
    # Store the widths of the bins used
    width_bins = []

    # Iterate for bin widths of: 25, 50, 75 and 100 px
    for d_b in range(25, 105, 25):
        # Number of bins in x,y given the bin width 'd_b'
        binsxy = [int((xmax - xmin) / d_b), int((ymax - ymin) / d_b)]
        
        # Store the values of the widths of the used bins
        width_bins.append(d_b)               
        
        # hist is the 2D histogran, xedges & yedges store the edges of the bins
        hist, xedges, yedges = np.histogram2d(x_data, y_data, range=rang, 
                                              bins=binsxy)
        # Store not-filtered 2D hist with d_b=20 in 'h_not_filt' (used to get
        # background value)
        h_not_filt.append(hist)
        
        # Only store the edges for the smallest value of 'd_b' <- HARDCODED!!
        if d_b == 25:
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
        # Store center bins for the filtered hist
        x_center_bin.append(x_cent_g)
        y_center_bin.append(y_cent_g)

        # Calculate centers in pixel coordinates.
        center_x_g.append(np.average(xedges[x_cent_g:x_cent_g + 2]))
        center_y_g.append(np.average(yedges[y_cent_g:y_cent_g + 2]))

    # Store in this array the values for the cluster's center x,y coordinates
    # obtained with the different bin sizes.
    arr_g = np.array([center_x_g, center_y_g])
    
    # Calculate the median value for the cluster's center (use median instead of
    # mean to reject possible outliers) and the standard deviation that will be
    # used as a measure of the error in the assignation of the cluster's center.
    center_coords, std_dev = np.median(arr_g, axis=1), np.std(arr_g, axis=1)
    cent_coo_err = [int(std_dev[0]), int(std_dev[1])]
    
    # Raise a flag if the median cluster's central coordinates are more than
    # 50 px away from the ones obtained with the min bin width.
    flag_center = False
    if abs(center_coords[0]-center_x_g[0]) > 50 or \
    abs(center_coords[1]-center_y_g[0]) > 50:
        flag_center = True
        
    # If the standard deviation for the cluster's center is bigger
    # than 50 px for any axis -> raise a flag.
    flag_std_dev = False
    if cent_coo_err[0] > 50 or cent_coo_err[1] > 50:
        flag_std_dev = True
        
    # Pass the center coordinates obtained with the smallest bin.
    center_coords[0], center_coords[1] = center_x_g[0], center_y_g[0]
    # Pass the errors as half of the width of the smallest bin.
    cent_coo_err[0], cent_coo_err[1] = int(round(width_bins[0]/2.)), \
    int(round(width_bins[0]/2.))

    return center_coords, cent_coo_err, h_filter, h_not_filt, xedges_min_db, \
    yedges_min_db, x_center_bin, y_center_bin, width_bins, flag_center, \
    flag_std_dev

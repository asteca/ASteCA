
import numpy as np
from scipy.ndimage.filters import gaussian_filter
import bisect


def cent_bin(xedges, yedges, xy_cent):
    '''
    Take x,y center coordinates and return the 2D bin where they are located.
    '''
    x_cent_bin = bisect.bisect_left(xedges, xy_cent[0])
    y_cent_bin = bisect.bisect_left(yedges, xy_cent[1])
    # Store center bin coords for the filtered hist.
    cent_bin = [(x_cent_bin - 1), (y_cent_bin - 1)]

    return cent_bin


def center_xy(hist, xedges, yedges, st_dev_lst):
    '''
    Function that returns a filtered 2D histogram and approximate center
    coordinates obtained using different standard deviations.
    '''

    cents_xy, hist_2d_g, cents_bin = [], [], []
    for h2d in hist:
        for st_dev in st_dev_lst:
            # 2D histogram with a Gaussian filter applied.
            h_g = gaussian_filter(h2d, st_dev, mode='constant')

            # x,y coordinates of the bin with the maximum value in the 2D
            # filtered histogram.
            x_cent_bin, y_cent_bin = np.unravel_index(h_g.argmax(), h_g.shape)
            # x,y coords of the center of the above bin.
            x_cent_pix = np.average(xedges[x_cent_bin:x_cent_bin + 2])
            y_cent_pix = np.average(yedges[y_cent_bin:y_cent_bin + 2])

            # Store for plotting.
            hist_2d_g.append(h_g)

            cents_xy.append([x_cent_pix, y_cent_pix])
            cents_bin.append(
                cent_bin(xedges, yedges, [x_cent_pix, y_cent_pix]))

    return cents_xy, hist_2d_g, cents_bin


def main(clp, center_stddev, **kwargs):
    """
    Obtain Gaussian filtered 2D x,y histograms and the maximum values in them
    as centers.
    """

    # Standard deviation values for the Gaussian filter.
    st_dev_lst = (center_stddev * .5, center_stddev, center_stddev * 2.)

    # Obtain center coordinates using Gaussian filters with different
    # standard deviation values, applied on the 2D (x,y) histogram.
    cents_xy, hist_2d_g, cents_bin_2d = center_xy(
        clp['hist_2d'], clp['xedges'], clp['yedges'], st_dev_lst)

    # Raise a flag if the standard deviation for either coordinate is larger
    # than 10% of that axis range. Use the full x,y positions list to
    # calculate the STDDEV.
    flag_center_std = False
    stddev = np.std(list(zip(*cents_xy[:3])), axis=1)
    if stddev[0] > 0.1 * np.ptp(clp['xedges']) or \
            stddev[1] > 0.1 * np.ptp(clp['yedges']):
        flag_center_std = True

    clp['flag_center_std'], clp['cents_xy'], clp['hist_2d_g'],\
        clp['cents_bin_2d'], clp['st_dev_lst'] = flag_center_std, cents_xy,\
        hist_2d_g, cents_bin_2d, st_dev_lst

    return clp


import numpy as np
import get_in_params as gip


def get_2d_histo(id_coords):
    '''
    Return list of 2D histograms for the positional data.
    '''

    x_data, y_data = id_coords[1:]
    xmin, xmax = min(x_data), max(x_data)
    ymin, ymax = min(y_data), max(y_data)
    rang = [[xmin, xmax], [ymin, ymax]]

    # Calculate the number of bins used.
    x_rang, y_rang = rang[0][1] - rang[0][0], rang[1][1] - rang[1][0]
    min_rang = min(x_rang, y_rang)

    # Bin width to create the 2D histogram.
    if gip.gh_params[0] == 'auto':
        bin_width = 1. * min_rang / 100.
    else:
        bin_width = gip.gh_params[1]

    # Number of bins in x,y given the bin width 'd_b'
    binsxy = [int((xmax - xmin) / bin_width), int((ymax - ymin) / bin_width)]

    # hist is the 2D histogran, *edges store the edges of the bins.
    hist_2d, xedges, yedges = np.histogram2d(x_data, y_data,
        range=rang, bins=binsxy)

    hist_lst = [hist_2d, xedges, yedges, bin_width]

    return hist_lst

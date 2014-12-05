
import numpy as np


def get_2d_histo(x_data, y_data, gh_params):
    '''
    Return list of 2D histograms for the positional data.
    '''

    xmin, xmax = min(x_data), max(x_data)
    ymin, ymax = min(y_data), max(y_data)

    # Calculate the number of bins used.
    x_rang, y_rang = (xmax - xmin), (ymax - ymin)

    # Bin width to create the 2D histogram.
    if gh_params[0] == 'auto':
        bin_width = min(x_rang, y_rang) / 100.
    else:
        bin_width = gh_params[1]

    # Number of bins in x,y given the bin width 'd_b'
    binsxy = [int(x_rang / bin_width), int(y_rang / bin_width)]

    # hist is the 2D histogran, *edges store the edges of the bins.
    hist_2d, xedges, yedges = np.histogram2d(x_data, y_data,
        range=[[xmin, xmax], [ymin, ymax]], bins=binsxy)

    hist_lst = [hist_2d, xedges, yedges, bin_width]

    return hist_lst

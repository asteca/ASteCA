
import numpy as np
from ..inp import input_params as g


def main(x_data, y_data):
    '''
    Return 2D histogram for the positional data.
    '''

    xmin, xmax = min(x_data), max(x_data)
    ymin, ymax = min(y_data), max(y_data)

    # Calculate the number of bins used.
    x_rang, y_rang = (xmax - xmin), (ymax - ymin)

    # Bin width to create the 2D histogram.
    if g.gh_params[0] == 'auto':
        bin_width = min(x_rang, y_rang) / 100.
    else:
        bin_width = g.gh_params[1]

    # Number of bins in x,y given the bin width.
    binsxy = [int(x_rang / bin_width), int(y_rang / bin_width)]

    # hist_2d is the 2D histogram, *edges store the edges of the bins.
    hist_2d, xedges, yedges = np.histogram2d(
        x_data, y_data, range=[[xmin, xmax], [ymin, ymax]], bins=binsxy)

    hist_lst = [hist_2d, xedges, yedges, bin_width]

    print "Frame's 2D histogram obtained"

    return hist_lst


import sys
import traceback
import numpy as np


def main(pd, x, y, **kwargs):
    '''
    Return 2D histogram for the x,y positional data.
    '''

    xmin, xmax = min(x), max(x)
    ymin, ymax = min(y), max(y)

    # Calculate the number of bins used.
    x_rang, y_rang = (xmax - xmin), (ymax - ymin)

    # Bin width to create the 2D histogram.
    if pd['histo_struct_method'] == 'auto':
        bin_width = min(x_rang, y_rang) / 100.
    else:
        bin_width = pd['histo_struct_bw']

    # Number of bins in x,y given the bin width.
    binsxy = [int(x_rang / bin_width), int(y_rang / bin_width)]

    try:
        # hist_2d is the 2D histogram, *edges store the edges of the bins.
        hist_2d, xedges, yedges = np.histogram2d(
            x, y, range=[[xmin, xmax], [ymin, ymax]], bins=binsxy)
    except ValueError:
        print traceback.format_exc()
        sys.exit("\nERROR: Could not generate the positional 2D histogram.\n"
                 "If the 'bin_width' was manually given in the 'Structure'\n"
                 "block of the 'params_input.dat' file, check that it is"
                 "not too large.\n")

    # Create cluster's parameters dictionary.
    clp = {'hist_2d': hist_2d, 'xedges': xedges, 'yedges': yedges,
           'bin_width': bin_width}

    print("Frame's 2D histogram obtained.")

    return clp

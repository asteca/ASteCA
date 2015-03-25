"""
@author: gabriel
"""

import matplotlib.pyplot as plt
from ..out import prep_plots


def disp_frame(x_data, y_data, mag_data):
    '''
    Show full frame.
    '''

    coord, x_name, y_name = prep_plots.coord_syst()
    st_sizes_arr = prep_plots.star_size(mag_data)

    # Get max and min values in x,y
    x_min, x_max = min(x_data), max(x_data)
    y_min, y_max = min(y_data), max(y_data)
    #Set plot limits
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    if coord == 'deg':
        # If RA is used, invert axis.
        plt.gca().invert_xaxis()
    #Set axis labels
    plt.xlabel('{} ({})'.format(x_name, coord), fontsize=12)
    plt.ylabel('{} ({})'.format(y_name, coord), fontsize=12)
    # Set minor ticks
    plt.minorticks_on()
    # Set grid
    plt.grid(b=True, which='major', color='k', linestyle='-', zorder=1)
    plt.grid(b=True, which='minor', color='k', linestyle='-', zorder=1)
    plt.scatter(x_data, y_data, marker='o', c='black', s=st_sizes_arr)

    plt.draw()
    print 'Plot displayed, waiting for it to be closed.'

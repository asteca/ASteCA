"""
@author: gabriel
"""

import numpy as np
import matplotlib.pyplot as plt


def disp_frame(x_data, y_data, mag_data):
    '''
    Show full frame.
    '''

    # Get max and min values in x,y
    x_min, x_max = min(x_data), max(x_data)
    y_min, y_max = min(y_data), max(y_data)
    #Set plot limits
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    #Set axis labels
    plt.xlabel('x (pixel)', fontsize=12)
    plt.ylabel('y (pixel)', fontsize=12)
    # Set minor ticks
    plt.minorticks_on()
    # Set grid
    plt.grid(b=True, which='major', color='k', linestyle='-', zorder=1)
    plt.grid(b=True, which='minor', color='k', linestyle='-', zorder=1)
    st_sizes_arr = 0.1 + 100. * 10 ** ((np.array(mag_data[0]) -
        min(mag_data[0])) / -2.5)
    plt.scatter(x_data, y_data, marker='o', c='black', s=st_sizes_arr)

    plt.draw()
    print 'Plot displayed, waiting for it to be closed.'

"""
@author: gabriel
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


def disp_cent(x_data, y_data, mag_data, center_cl, cent_bin, h_filter):
    '''
    Show plot of cluster with value of center obtained.
    '''

    # Plot all outputs
    plt.figure(figsize=(18, 8))  # create the top-level container
    gs = gridspec.GridSpec(1, 2)  # create a GridSpec object

    # 1 subplot: 2D filtered histogram
    ax1 = plt.subplot(gs[0, 0])
    #Set axis labels
    plt.xlabel('x', fontsize=12)
    plt.ylabel('y', fontsize=12)
    # Set minor ticks
    ax1.minorticks_on()
    # Set grid
    ax1.grid(b=True, which='major', color='k', linestyle='--', zorder=3)
    ax1.grid(b=True, which='minor', color='k', linestyle='--', zorder=3)
    plt.axvline(x=cent_bin[0], linestyle='--', color='white')
    plt.axhline(y=cent_bin[1], linestyle='--', color='white')
    plt.imshow(h_filter.transpose(), origin='lower')

    # 2 subplot: x,y finding chart of full frame
    ax2 = plt.subplot(gs[0, 1])
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
    ax2.minorticks_on()
    # Set grid
    ax2.grid(b=True, which='major', color='k', linestyle='-', zorder=1)
    ax2.grid(b=True, which='minor', color='k', linestyle='-', zorder=1)
    # Add lines through the center of the cluster
    plt.axvline(x=center_cl[0], linestyle='--', color='red', lw=2.)
    plt.axhline(y=center_cl[1], linestyle='--', color='red', lw=2.)
    # Add text box
    text1 = r'$x_{cent} = %d$' % (center_cl[0])
    text2 = '\n'
    text3 = r'$y_{cent} = %d$' % (center_cl[1])
    text = text1 + text2 + text3
    plt.text(0.7, 0.9, text, transform=ax2.transAxes,
    bbox=dict(facecolor='white', alpha=0.85), fontsize=15)
    st_sizes_arr = 0.1 + 100. * 10 ** ((np.array(mag_data) -
        min(mag_data)) / -2.5)
    plt.scatter(x_data, y_data, marker='o', c='black', s=st_sizes_arr)

    plt.draw()
    print 'Plot displayed, waiting for it to be closed.'

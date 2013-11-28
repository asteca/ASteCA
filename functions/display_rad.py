"""
@author: gabriel
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle

def disp_rad(x_data, y_data, mag_data, center_cl, cent_cl_err, clust_rad, \
    x_center_bin, y_center_bin, h_filter, backg_value, radii, delta_backg,
    ring_density, clust_name, poisson_error, width_bins, delta_percentage,
    inner_ring, outer_ring):
    '''
    Plot cluster and its radius.
    '''
    
    # Plot all outputs
    fig = plt.figure(figsize=(12, 12)) # create the top-level container
    gs = gridspec.GridSpec(2, 2)  # create a GridSpec object
        
    # 1 subplot: 2D filtered histogram, d_b=25
    ax1 = plt.subplot(gs[0, 0])
    #Set axis labels
    plt.xlabel('x', fontsize=12)
    plt.ylabel('y', fontsize=12)
    # Set minor ticks
    ax1.minorticks_on()
    # Add lines through the center of the cluster
    plt.axvline(x=x_center_bin[0], linestyle='--', color='white')
    plt.axhline(y=y_center_bin[0], linestyle='--', color='white')
    # Radius
    circle = plt.Circle((x_center_bin[0], y_center_bin[0]), clust_rad[0]/25., 
                        color='w', fill=False)
    fig.gca().add_artist(circle)
    # Cluster's name in a text box
    text = 'Center bins (%d, %d)' % (x_center_bin[0], y_center_bin[0])
    plt.text(0.5, 0.95, text, transform = ax1.transAxes, \
    bbox=dict(facecolor='white', alpha=0.8), fontsize=12)
    plt.imshow(h_filter[0].transpose(), origin='lower')

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
    # Draw circle radius
    circle = plt.Circle((center_cl[0], center_cl[1]), clust_rad[0], color='r', 
                        fill=False)
    fig.gca().add_artist(circle)
    if backg_value[0] > 0.005:
        # High density field
        plt.scatter(x_data, y_data, marker='o', c='black', 
                    s=500*np.exp(-0.004*mag_data**2.5))
    else:
        # Low density field
        plt.scatter(x_data, y_data, marker='o', c='black', 
                    s=500*np.exp(-0.0035*mag_data**2.5))
    # Plot inner and outer rectangles used to calculate the background value. 
    # Inner ring
    boxes = plt.gca()
    boxes.add_patch(Rectangle(((center_cl[0]-inner_ring), 
                               (center_cl[1]-inner_ring)), inner_ring*2,\
                               inner_ring*2, facecolor='none', edgecolor='b',\
                               ls='dashed'))
    # Outer ring
    boxes.add_patch(Rectangle(((center_cl[0]-outer_ring), 
                               (center_cl[1]-outer_ring)), outer_ring*2,\
                               outer_ring*2, facecolor='none', edgecolor='b',\
                               ls='dashed'))


    # Density profile
    ax4 = plt.subplot(gs[1, 0:2])
    # Get max and min values in x,y
    x_max = max(radii[0])+10
    y_min, y_max = (backg_value[0]-delta_backg)-(max(ring_density[0])-\
    min(ring_density[0]))/10, max(ring_density[0])+(max(ring_density[0])-\
    min(ring_density[0]))/10
    # Set plot limits
    plt.xlim(-10, x_max)
    plt.ylim(y_min, y_max)
    # Set axes labels
    plt.xlabel('radius (px)', fontsize=12)
    plt.ylabel(r"stars/px$^{2}$", fontsize=12)
    # Cluster's name.
    text = str(clust_name)
    plt.text(0.5, 0.9, text, transform = ax4.transAxes, fontsize=14)
    # Plot poisson error bars
    plt.errorbar(radii[0], ring_density[0], yerr = poisson_error[0], fmt='ko', 
                 zorder=1)
    ## Plot the delta around the background value used to asses when the density
    # has become stable
    plt.hlines(y=(backg_value[0]-delta_backg), xmin=0, xmax=max(radii[0]), 
               color='k', linestyles='dashed', zorder=2)
    plt.hlines(y=(backg_value[0]+delta_backg), xmin=0, xmax=max(radii[0]), 
               color='k', linestyles='dashed', zorder=2)
    # Legend texts
    texts = ['Density Profile (%d px)' % width_bins[0], 'Density Profile \
(%d px)' % width_bins[1], \
    'Density Profile (%d px)' % width_bins[2], 'Density Profile (%d px)' % \
    width_bins[3]]
    # Plot density profile with the smallest bin size
    ax4.plot(radii[0], ring_density[0], 'ko-', zorder=3, label=texts[0])
    # Plot density profiles for the last 3 (bigger) bin widths
    ax4.plot(radii[1], ring_density[1], 'go--', zorder=6, label=texts[1])
    ax4.plot(radii[2], ring_density[2], 'bo--', zorder=7, label=texts[2])
    ax4.plot(radii[3], ring_density[3], 'ro--', zorder=8, label=texts[3])
    # Plot background and radius lines
    ax4.vlines(x=clust_rad[0], ymin=0, ymax=max(ring_density[0]), label='Radius\
 = %d px' % clust_rad[0], color='r', zorder=4)
    ax4.hlines(y=backg_value[0], xmin=0, xmax=max(radii[0]), 
               label=r'Background ($\Delta$=%d%%)' % delta_percentage, \
               color='k', zorder=5)
    # get handles
    handles, labels = ax4.get_legend_handles_labels()
    # use them in the legend
    ax4.legend(handles, labels, loc='upper right', numpoints=1, fontsize=10)
    

    plt.draw()
    print 'Plot displayed, waiting for it to be closed.'

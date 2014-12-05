"""
@author: gabriel
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


def disp_rad(phot_data, bin_width, center_params, clust_rad, backg_value,
    rdp_params):
    '''
    Plot cluster and its radius.
    '''

    # Unpack.
    x_data, y_data, mag_data, col1_data = phot_data[1], phot_data[2], \
    phot_data[3], phot_data[5]

    cent_bin, kde_cent, e_cent = center_params[:3]
    hist_2d_g = center_params[5]

    center_cl = [kde_cent[0], kde_cent[1]]
    #h_not_filt, hist_xyedges, h_filter = center_params[:5]
    #xedges_min_db, yedges_min_db = hist_xyedges
    x_center_bin, y_center_bin = cent_bin

    #bin_width, cent_cl_err = bin_list[0], bin_list[0]
    radii, ring_density, poisson_error = rdp_params[:3]

    # Plot all outputs
    fig = plt.figure(figsize=(12, 12))
    gs = gridspec.GridSpec(2, 2)  # create a GridSpec object

    # 2D filtered histogram, smallest bin width.
    ax1 = plt.subplot(gs[0, 0])
    plt.xlabel('x (bins)', fontsize=12)
    plt.ylabel('y (bins)', fontsize=12)
    ax1.minorticks_on()
    plt.axvline(x=x_center_bin, linestyle='--', color='white')
    plt.axhline(y=y_center_bin, linestyle='--', color='white')
    # Radius
    circle = plt.Circle((x_center_bin, y_center_bin),
        clust_rad / bin_width, color='w', fill=False)
    fig.gca().add_artist(circle)
    text = 'Bin: %d px' % (bin_width)
    plt.text(0.7, 0.92, text, transform=ax1.transAxes,
             bbox=dict(facecolor='white', alpha=0.8), fontsize=12)
    plt.imshow(hist_2d_g.transpose(), origin='lower')

    # Finding chart.
    ax2 = plt.subplot(gs[0, 1])
    # Get max and min values in x,y
    x_min, x_max = min(x_data), max(x_data)
    y_min, y_max = min(y_data), max(y_data)
    #Set plot limits
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    #Set axis labels
    plt.xlabel('x (px)', fontsize=12)
    plt.ylabel('y (px)', fontsize=12)
    # Set minor ticks
    ax2.minorticks_on()
    circle = plt.Circle((center_cl[0], center_cl[1]), clust_rad, color='r',
                        fill=False)
    fig.gca().add_artist(circle)
    # Add text box
    text1 = '$x_{cent} = %d \pm %d px$' % (center_cl[0], e_cent[0])
    text2 = '\n'
    text3 = '$y_{cent} = %d \pm %d px$' % (center_cl[1], e_cent[1])
    text4 = text1 + text2 + text3
    plt.text(0.5, 0.85, text4, transform=ax2.transAxes,
        bbox=dict(facecolor='white', alpha=0.85), fontsize=15)
    # Plot stars.
    a, b, c = 50, -0.003, 2.5
    plt.scatter(x_data, y_data, marker='o', c='black',
        s=a * np.exp(b * mag_data ** c))

    # Radial density plot.
    ax3 = plt.subplot(gs[1, 0:2])
    # Get max and min values in x,y
    x_min, x_max = min(radii) - (max(radii) / 10.), \
    max(radii) + (max(radii) / 10.)
    max(radii) + (max(radii) / 20.)
    delta_total = (max(ring_density) - backg_value)
    delta_backg = 0.2 * delta_total
    y_min, y_max = (backg_value - delta_backg) - (max(ring_density) -
    min(ring_density)) / 10, max(ring_density) + (max(ring_density) -
    min(ring_density)) / 10
    # Set plot limits
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    # Set axes labels
    plt.xlabel('radius (px)', fontsize=12)
    plt.ylabel("stars/px$^{2}$", fontsize=12)
    # Cluster's name.
    #text = str(clust_name)
    #plt.text(0.4, 0.9, text, transform=ax3.transAxes, fontsize=14)
    # Legend texts
    texts = ['Dens prof (%d px)' % bin_width,
            'backg = %.1E $st/px^{2}$' % backg_value,
            'r$_{cl}$ = %d $\pm$ %d px' % (clust_rad, round(bin_width))]
    # Plot density profile with the smallest bin size
    ax3.plot(radii, ring_density, 'ko-', zorder=3, label=texts[0])
    # Plot poisson error bars
    plt.errorbar(radii, ring_density, yerr=poisson_error, fmt='ko',
                 zorder=1)
    # Plot background level.
    ax3.hlines(y=backg_value, xmin=0, xmax=max(radii),
               label=texts[1], color='b', zorder=5)
    # Approx middle of the graph.
    arr_y_up = (y_max - y_min) / 2.3 + y_min
    # Length of arrow head.
    #head_w = abs((arr_y_up - backg_value) / 10.)
    head_w, head_l = x_max * 0.023, (y_max - y_min) * 0.045
    # Length of arrow.
    arr_y_dwn = -1. * abs(arr_y_up - backg_value) * 0.76
    # Plot radius.
    ax3.vlines(x=clust_rad, ymin=0, ymax=0., label=texts[2], color='r')
    ax3.arrow(clust_rad, arr_y_up, 0., arr_y_dwn, fc="r",
              ec="r", head_width=head_w, head_length=head_l, zorder=5)
    # get handles
    handles, labels = ax3.get_legend_handles_labels()
    # use them in the legend
    ax3.legend(handles, labels, loc='upper right', numpoints=2, fontsize=12)
    ax3.minorticks_on()

    plt.draw()
    print 'Plot displayed, waiting for it to be closed.'

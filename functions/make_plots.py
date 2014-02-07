"""
@author: gabriel
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle
from matplotlib.ticker import MultipleLocator
from itertools import cycle
from scipy.stats import norm
import matplotlib.mlab as mlab
        
from os.path import join


def make_plots(output_subdir, clust_name, x_data, y_data, center_cl, 
               cent_cl_err, x_center_bin, y_center_bin, h_filter, radii, 
               backg_value, inner_ring, outer_ring, delta_backg, ring_density,
               poisson_error, clust_rad, cont_index, width_bins,
               delta_percentage, mag_data, col1_data, bright_end, popt_mag,
               popt_umag, pol_mag, popt_col1, popt_ucol1, pol_col1,
               mag_val_left, mag_val_right, col1_val_left,
               col1_val_right, use_errors_fit, k_prof, k_pr_err,
               flag_king_no_conver, stars_in,
               stars_out, stars_in_rjct, stars_out_rjct, stars_in_mag,
               stars_in_all_mag, n_c, flag_area_stronger,
               cluster_region, field_region, pval_test_params, qq_params,
               clust_reg_prob_avrg, field_reg_box,
               kde_cl, kde_f, memb_prob_avrg_sort, iso_moved, zams_iso,
               cl_e_bv, cl_age, cl_feh, cl_dmod):
    '''
    Make all plots.
    '''
    
    def func(x, a, b, c):
        '''Exponential function.
        '''
        return a * np.exp(b * x) + c
        
    
    def line(x,slope, intercept):
        '''
        Linar function.
        '''
        y = slope*x + intercept
        return y

    def three_params(x):
        '''
        Three parameters King profile fit.
        '''
        a, b, c, d = k_prof[0], k_prof[1], k_prof[2], backg_value[0]
        return c*(1/np.sqrt(1+(x/a)**2) - 1/np.sqrt(1+(b/a)**2))**2 + d
        

    # Define plot limits for ALL CMD diagrams.
    col1_min, col1_max = max(-0.9, min(col1_data)-0.2),\
                         min(3.9, max(col1_data)+0.2)
    mag_min, mag_max = max(mag_data)+0.5, min(mag_data)-0.5



    # Plot all outputs
    # figsize(x1, y1), GridSpec(y2, x2) --> To have square plots: x1/x2 = 
    # y1/y2 = 2.5 
    fig = plt.figure(figsize=(20, 30)) # create the top-level container
    gs1 = gridspec.GridSpec(12, 8)  # create a GridSpec object
    #gs1.update(wspace=.09, hspace=.0)

    # 2D filtered histogram, d_b=100
    ax0 = plt.subplot(gs1[0:2, 0:2])
    #Set axis labels
    plt.xlabel('x (bins)', fontsize=12)
    plt.ylabel('y (bins)', fontsize=12)
    # Set minor ticks
    ax0.minorticks_on()
    # Set ticks font.
    #ax0.tick_params(axis='both', which='major', labelsize=7)
    # Set grid
    ax0.grid(b=True, which='major', color='k', linestyle='-', zorder=1)
    ax0.grid(b=True, which='minor', color='k', linestyle='-', zorder=1)
    # Add lines through the center of the cluster
    plt.axvline(x=x_center_bin[3], linestyle='-', color='white', linewidth=2, 
                zorder=3)
    plt.axhline(y=y_center_bin[3], linestyle='-', color='white', linewidth=2, 
                zorder=3)
    # Plot center coordinates in bin.
    #text = 'Center (%d, %d)' % (x_center_bin[3], y_center_bin[3])
    text = 'Bin: %d px' % (width_bins[3])
    plt.text(0.7, 0.92, text, transform = ax0.transAxes, 
             bbox=dict(facecolor='white', alpha=0.8), fontsize=12)
    plt.imshow(h_filter[3].transpose(), origin='lower')
        
        
    # 2D filtered histogram, d_b=75
    ax1 = plt.subplot(gs1[0:2, 2:4])
    plt.xlabel('x (bins)', fontsize=12)
    plt.ylabel('y (bins)', fontsize=12)
    ax1.minorticks_on()
    plt.axvline(x=x_center_bin[2], linestyle='--', color='white')
    plt.axhline(y=y_center_bin[2], linestyle='--', color='white')
    #text = 'Center (%d, %d)' % (x_center_bin[2], y_center_bin[2])
    text = 'Bin: %d px' % (width_bins[2])
    plt.text(0.7, 0.92, text, transform = ax1.transAxes, 
             bbox=dict(facecolor='white', alpha=0.8), fontsize=12)
    plt.imshow(h_filter[2].transpose(), origin='lower')
        
        
    # 2D filtered histogram, d_b=50
    ax2 = plt.subplot(gs1[0:2, 4:6])
    plt.xlabel('x (bins)', fontsize=12)
    plt.ylabel('y (bins)', fontsize=12)
    ax2.minorticks_on()
    plt.axvline(x=x_center_bin[1], linestyle='--', color='white')
    plt.axhline(y=y_center_bin[1], linestyle='--', color='white')
    #text = 'Center (%d, %d)' % (x_center_bin[1], y_center_bin[1])
    text = 'Bin: %d px' % (width_bins[1])
    plt.text(0.7, 0.92, text, transform = ax2.transAxes, 
             bbox=dict(facecolor='white', alpha=0.8), fontsize=12)
    plt.imshow(h_filter[1].transpose(), origin='lower')
        
        
    # 2D filtered histogram, d_b=25
    ax3 = plt.subplot(gs1[0:2, 6:8])
    plt.xlabel('x (bins)', fontsize=12)
    plt.ylabel('y (bins)', fontsize=12)
    ax3.minorticks_on()
    plt.axvline(x=x_center_bin[0], linestyle='--', color='white')
    plt.axhline(y=y_center_bin[0], linestyle='--', color='white')
    # Radius
    circle = plt.Circle((x_center_bin[0], y_center_bin[0]), clust_rad[0]/25., 
                        color='w', fill=False)
    fig.gca().add_artist(circle)
    #text = 'Center (%d, %d)' % (x_center_bin[0], y_center_bin[0])
    text = 'Bin: %d px' % (width_bins[0])
    plt.text(0.7, 0.92, text, transform = ax3.transAxes, 
             bbox=dict(facecolor='white', alpha=0.8), fontsize=12)
    plt.imshow(h_filter[0].transpose(), origin='lower')


    # x,y finding chart of full frame
    ax5 = plt.subplot(gs1[2:4, 0:2])
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
    ax5.minorticks_on()
    circle = plt.Circle((center_cl[0], center_cl[1]), clust_rad[0], color='r', 
                        fill=False)
    fig.gca().add_artist(circle)
    # Add text box
    text1 = r'$x_{cent} = %d \pm %d px$' % (center_cl[0], cent_cl_err[0])
    text2 = '\n'
    text3 = r'$y_{cent} = %d \pm %d px$' % (center_cl[1], cent_cl_err[1])
    text4 = text1+text2+text3
    plt.text(0.5, 0.85, text4, transform = ax5.transAxes, \
    bbox=dict(facecolor='white', alpha=0.85), fontsize=15)
    # Count the number of very bright stars in the field.
    range_10_perc = (max(mag_data)-min(mag_data))/10.+min(mag_data)
    bright_stars = len([i for i in mag_data if i < range_10_perc])
    # Set exponential factor for high and low density fields.
    exp_factor = -0.004 if backg_value[0] > 0.005 else -0.0035
    # Set a lower amplitude for fields with very bright stars.
    amplitude = 500 if bright_stars < 10 else 200
    plt.scatter(x_data, y_data, marker='o', c='black',
                s = amplitude*np.exp(exp_factor*mag_data**2.5))
    # Plot inner and outer rectangles used to calculate the background value. 
    # Inner ring
    boxes = plt.gca()
    boxes.add_patch(Rectangle(((center_cl[0]-inner_ring), 
                               (center_cl[1]-inner_ring)), inner_ring*2,\
                               inner_ring*2, facecolor='none', edgecolor='b',\
                               ls='dashed', lw=1.5))
    # Outer ring
    boxes.add_patch(Rectangle(((center_cl[0]-outer_ring), 
                               (center_cl[1]-outer_ring)), outer_ring*2,\
                               outer_ring*2, facecolor='none', edgecolor='b',\
                               ls='dashed', lw=1.5))

        
    # Density profile
#    ax4 = plt.subplot(gs1[2:4, 2:6])
#    # Get max and min values in x,y
#    x_max = max(radii[0])+10
#    y_min, y_max = (backg_value[0]-delta_backg)-(max(ring_density[0])-\
#    min(ring_density[0]))/10, max(ring_density[0])+(max(ring_density[0])-\
#    min(ring_density[0]))/10
#    # Set plot limits
#    plt.xlim(-10, x_max)
#    plt.ylim(y_min, y_max)
#    # Set axes labels
#    plt.xlabel('radius (px)', fontsize=12)
#    plt.ylabel(r"stars/px$^{2}$", fontsize=12)
#    # Cluster's name.
#    text = str(clust_name)
#    plt.text(0.4, 0.9, text, transform = ax4.transAxes, fontsize=14)
#    # Plot poisson error bars
#    plt.errorbar(radii[0], ring_density[0], yerr = poisson_error[0], fmt='ko', 
#                 zorder=1)
#    ## Plot the delta around the background value used to asses when the density
#    # has become stable
#    plt.hlines(y=(backg_value[0]-delta_backg), xmin=0, xmax=max(radii[0]), 
#               color='k', linestyles='dashed', zorder=2)
#    plt.hlines(y=(backg_value[0]+delta_backg), xmin=0, xmax=max(radii[0]), 
#               color='k', linestyles='dashed', zorder=2)
#    # Legend texts
#    texts = ['Density Profile (%d px)' % width_bins[0], 'Density Profile \
#(%d px)' % width_bins[1], \
#    'Density Profile (%d px)' % width_bins[2], 'Density Profile (%d px)' % \
#    width_bins[3]]
#    # Plot density profile with the smallest bin size
#    ax4.plot(radii[0], ring_density[0], 'ko-', zorder=3, label=texts[0])
#    # Plot density profiles for the last 3 (bigger) bin widths
#    ax4.plot(radii[1], ring_density[1], 'go--', zorder=6, label=texts[1])
#    ax4.plot(radii[2], ring_density[2], 'bo--', zorder=7, label=texts[2])
#    ax4.plot(radii[3], ring_density[3], 'ro--', zorder=8, label=texts[3])
#    # Plot background and radius lines
#    ax4.vlines(x=clust_rad[0], ymin=0, ymax=max(ring_density[0]),
#               label=r'Radius = %d $\pm$ %d px' % (clust_rad[0],\
#               round(width_bins[0]/2.)), color='r', zorder=4)
#    ax4.hlines(y=backg_value[0], xmin=0, xmax=max(radii[0]), 
#               label=r'Background ($\Delta$=%d%%)' % delta_percentage, \
#               color='k', zorder=5)
#    # get handles
#    handles, labels = ax4.get_legend_handles_labels()
#    # use them in the legend
#    ax4.legend(handles, labels, loc='upper right', numpoints=1, fontsize=10)


    ax4 = plt.subplot(gs1[2:4, 2:6])
    # Get max and min values in x,y
    x_max = max(radii[0])+10
    y_min, y_max = (backg_value[0]-delta_backg)-(max(ring_density[0])-\
    min(ring_density[0]))/10, max(ring_density[0])+(max(ring_density[0])-\
    min(ring_density[0]))/10
    # Set plot limits
    plt.xlim(-10, min(x_max, 500))
    plt.ylim(y_min, y_max)
    # Set axes labels
    plt.xlabel('radius (px)', fontsize=12)
    plt.ylabel(r"stars/px$^{2}$", fontsize=12)
    # Cluster's name.
    text = str(clust_name)
    plt.text(0.4, 0.9, text, transform = ax4.transAxes, fontsize=14)
    # Plot poisson error bars
    plt.errorbar(radii[0], ring_density[0], yerr = poisson_error[0], fmt='ko', 
                 zorder=1)
    ## Plot the delta around the background value used to asses when the density
    # has become stable
#    plt.hlines(y=(backg_value[0]-delta_backg), xmin=0, xmax=max(radii[0]), 
#               color='k', linestyles='dashed', zorder=2)
    plt.hlines(y=(backg_value[0]+delta_backg), xmin=0, xmax=max(radii[0]), 
               color='k', linestyles='dashed', zorder=2)
    # Legend texts
    # Calculate errors for fitted King parameters. Use sqrt since the error
    # given by 'curve_fit' is the variance and we want the standard deviation.
    if flag_king_no_conver == False:
        rc_err = round(np.sqrt(k_pr_err[0][0]))
        rt_err = round(np.sqrt(k_pr_err[1][1]))
    else:
        rc_err, rt_err = -1, -1
    texts = ['Density Profile (%d px)' % width_bins[0],
             r'$R_{cl}$ = %d $\pm$ %d px' % (clust_rad[0],
                                           round(width_bins[0]/2.)),\
            '3-P King profile',
             r'$R_c$ = %d $\pm$ %d px' % (k_prof[0], rc_err),\
             r'$R_t$ = %d $\pm$ %d px' % (k_prof[1], rt_err),\
             r'Background ($\Delta$=%d%%)' % delta_percentage]
    # Plot density profile with the smallest bin size
    ax4.plot(radii[0], ring_density[0], 'ko-', zorder=3, label=texts[0])
    # Plot King profile.
    arr_y_up = (y_max-y_min)/2. + y_min # Middle of the graph.
    head_l = abs((arr_y_up-backg_value[0])/7.) # Length of arrow head.
    arr_y_dwn = -1.*abs(arr_y_up-backg_value[0]-head_l*1.5) # Length of arrow.
    if flag_king_no_conver == False:
        ax4.plot(radii[0], three_params(radii[0]), 'b--', label=texts[2],
                 lw=1.5, zorder=3)
        # Plot R_c as a dashed line.
        ax4.vlines(x=k_prof[0], ymin=0, ymax=max(ring_density[0])/1.2,
                   label=texts[3], color='g', linestyles='dashed', zorder=4)
        # Plot R_t radius as an arrow. vline is there to show the label.
        ax4.vlines(x=k_prof[1], ymin=0., ymax=0., label=texts[4], color='b')
        ax4.arrow(k_prof[1], arr_y_up, 0., arr_y_dwn, fc="b", ec="b",\
                  head_width=10., head_length=head_l, zorder=5)
    # Plot radius.
    ax4.vlines(x=clust_rad[0], ymin=0, ymax=0., label=texts[1], color='r')
    ax4.arrow(clust_rad[0], arr_y_up, 0., arr_y_dwn, fc="r",
              ec="r", head_width=10., head_length=head_l, zorder=5)
    # Plot background level.
    ax4.hlines(y=backg_value[0], xmin=0, xmax=max(radii[0]), 
               label=texts[5], color='k', zorder=5)
    # get handles
    handles, labels = ax4.get_legend_handles_labels()
    # use them in the legend
    ax4.legend(handles, labels, loc='upper right', numpoints=1, fontsize=10)
    ax4.minorticks_on()


    # Zoom on x,y finding chart
    ax6 = plt.subplot(gs1[2:4, 6:8])
    #Set plot limits
    x_min, x_max = min(x_data), max(x_data)
    y_min, y_max = min(y_data), max(y_data)
    # If possible, zoom in.
    x_min, x_max = max(x_min, (center_cl[0]-1.5*clust_rad[0])), \
    min(x_max, (center_cl[0]+1.5*clust_rad[0]))
    y_min, y_max = max(y_min, (center_cl[1]-1.5*clust_rad[0])), \
    min(y_max, (center_cl[1]+1.5*clust_rad[0]))
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    #Set axis labels
    plt.xlabel('x (px)', fontsize=12)
    plt.ylabel('y (px)', fontsize=12)
    # Set minor ticks
    ax6.minorticks_on()
    # Add circle
    circle = plt.Circle((center_cl[0], center_cl[1]), clust_rad[0], color='r', 
                        fill=False)
    fig.gca().add_artist(circle)
    text1 = 'Cluster (zoom)\n'
    text2 = 'CI = %0.2f' % (cont_index)
    text = text1 + text2
    plt.text(0.6, 0.9, text, transform = ax6.transAxes, 
             bbox=dict(facecolor='white', alpha=0.85), fontsize=12)
    if backg_value[0] > 0.005:
        plt.scatter(x_data, y_data, marker='o', c='black', 
                    s=500*np.exp(-0.003*mag_data**2.5))
    else:
        plt.scatter(x_data, y_data, marker='o', c='black', 
                    s=500*np.exp(-0.0025*mag_data**2.5))
        
        
    # T1 magnitude error
    ax7 = plt.subplot(gs1[4, 0:2])
    #Set plot limits
    plt.xlim(min(mag_data)-0.5, max(mag_data)+0.5)
    plt.ylim(-0.005, 0.31)
    #Set axis labels
    plt.ylabel(r'$\sigma_{T_1}$', fontsize=18)
    plt.xlabel(r'$T_1$', fontsize=18)
    # Set minor ticks
    ax7.minorticks_on()
    mag_x = np.linspace(min(mag_data), max(mag_data), 50)
    # Condition to not plot the lines if the fit was rejected.
    if use_errors_fit:
        # Plot lower envelope.
        ax7.plot(mag_x, func(mag_x, *popt_mag), 'r-', zorder=3)
        # Plot left side of upper envelope (exponential).
        ax7.plot(mag_val_left, func(mag_val_left, *popt_umag), 'r--', lw=2.,
                 zorder=3)
        # Plot right side of upper envelope (polynomial).
        ax7.plot(mag_val_right, np.polyval(pol_mag, (mag_val_right)),
                 'r--', lw=2., zorder=3)
        # Plot rectangle.
        ax7.vlines(x=bright_end+0.05, ymin=-0.005, ymax=0.1, color='r', 
                   linestyles='dashed', zorder=2)
        ax7.vlines(x=min(mag_data)-0.05, ymin=-0.005, ymax=0.1, color='r',
                   linestyles='dashed', zorder=2)
        ax7.hlines(y=0.1, xmin=min(mag_data), xmax=bright_end, color='r',
                   linestyles='dashed', zorder=2)
    # Plot stars.
    stars_rjct_temp = [[], []]
    for star in stars_out_rjct:
        stars_rjct_temp[0].append(star[3])
        stars_rjct_temp[1].append(star[4])
    for star in stars_in_rjct:
        stars_rjct_temp[0].append(star[3])
        stars_rjct_temp[1].append(star[4])
    plt.scatter(stars_rjct_temp[0], stars_rjct_temp[1], marker='x', c='teal', 
                s=15, zorder=1)
    stars_acpt_temp = [[], []]
    for star in stars_out:
        stars_acpt_temp[0].append(star[3])
        stars_acpt_temp[1].append(star[4])
    for star in stars_in:
        stars_acpt_temp[0].append(star[3])
        stars_acpt_temp[1].append(star[4])
    plt.scatter(stars_acpt_temp[0], stars_acpt_temp[1], marker='o', c='k', 
                s=1, zorder=2)


    # C-T1 color error
    ax8 = plt.subplot(gs1[5, 0:2])
    #Set plot limits
    plt.xlim(min(mag_data)-0.5, max(mag_data)+0.5)
    plt.ylim(-0.005, 0.31)
    #Set axis labels
    plt.ylabel(r'$\sigma_{(C-T_1)}$', fontsize=18)
    plt.xlabel(r'$T_1$', fontsize=18)
    # Set minor ticks
    ax8.minorticks_on()
    # Condition to not plot the lines if the fit was rejected.
    if use_errors_fit:
        # Plot lower envelope.
        ax8.plot(mag_x, func(mag_x, *popt_col1), 'r-', zorder=3)
        # Plot left side of upper envelope (exponential).
        ax8.plot(col1_val_left, func(col1_val_left, *popt_ucol1), 'r--', lw=2.,
                 zorder=3)
        # Plot right side of upper envelope (polynomial).
        ax8.plot(col1_val_right, np.polyval(pol_col1, (col1_val_right)),
                 'r--', lw=2., zorder=3)
        # Plot rectangle.
        ax8.vlines(x=bright_end+0.05, ymin=-0.005, ymax=0.1, color='r', 
                   linestyles='dashed', zorder=2)
        ax8.vlines(x=min(mag_data)-0.05, ymin=-0.005, ymax=0.1, color='r',
                   linestyles='dashed', zorder=2)
        ax8.hlines(y=0.1, xmin=min(mag_data), xmax=bright_end, color='r',
                   linestyles='dashed', zorder=2)
    # Plot stars.
    stars_rjct_temp = [[], []]
    for star in stars_out_rjct:
        stars_rjct_temp[0].append(star[3])
        stars_rjct_temp[1].append(star[6])
    for star in stars_in_rjct:
        stars_rjct_temp[0].append(star[3])
        stars_rjct_temp[1].append(star[6])
    plt.scatter(stars_rjct_temp[0], stars_rjct_temp[1], marker='x', c='teal', 
                s=15, zorder=1)
    stars_acpt_temp = [[], []]
    for star in stars_out:
        stars_acpt_temp[0].append(star[3])
        stars_acpt_temp[1].append(star[6])
    for star in stars_in:
        stars_acpt_temp[0].append(star[3])
        stars_acpt_temp[1].append(star[6])
    plt.scatter(stars_acpt_temp[0], stars_acpt_temp[1], marker='o', c='k', 
                s=1, zorder=2)


    # Field stars CMD (stars outside cluster's radius)
    ax9 = plt.subplot(gs1[4:6, 2:4])
    #Set plot limits
    plt.xlim(col1_min, col1_max)
    plt.ylim(mag_min, mag_max)
    #Set axis labels
    plt.xlabel(r'$C-T_1$', fontsize=18)
    plt.ylabel(r'$T_1$', fontsize=18)
    # Set minor ticks
    ax9.minorticks_on()
    # Only draw units on axis (ie: 1, 2, 3)
    ax9.xaxis.set_major_locator(MultipleLocator(1.0))
    # Set grid
    ax9.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
    plt.text(0.8, 0.93, r'$r > R_{cl}$', transform = ax9.transAxes, \
    bbox=dict(facecolor='white', alpha=0.5), fontsize=16)
    # Plot stars.
    stars_rjct_temp = [[], []]
    for star in stars_out_rjct:
        stars_rjct_temp[0].append(star[5])
        stars_rjct_temp[1].append(star[3])
    plt.scatter(stars_rjct_temp[0], stars_rjct_temp[1], marker='x', c='teal', 
                s=10, zorder=1)
    stars_acpt_temp = [[], []]
    for star in stars_out:
        stars_acpt_temp[0].append(star[5])
        stars_acpt_temp[1].append(star[3])
    sz_pt = 0.2 if (len(stars_out_rjct)+len(stars_out)) > 5000 else 0.5
    plt.scatter(stars_acpt_temp[0], stars_acpt_temp[1], marker='o', c='k', 
                s=sz_pt, zorder=2)


    # Cluster's stars CMD (stars inside cluster's radius)
    ax10 = plt.subplot(gs1[4:6, 4:6])
    #Set plot limits
    plt.xlim(col1_min, col1_max)
    plt.ylim(mag_min, mag_max)
    #Set axis labels
    plt.xlabel(r'$C-T_1$', fontsize=18)
    plt.ylabel(r'$T_1$', fontsize=18)
    # Set minor ticks
    ax10.minorticks_on()
    # Only draw units on axis (ie: 1, 2, 3)
    ax10.xaxis.set_major_locator(MultipleLocator(1.0))
    # Set grid
    ax10.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
    # Calculate total number of stars whitin cluster's radius.
    tot_stars = len(stars_in_rjct) + len(stars_in)
    plt.text(0.55, 0.93, r'$r \leq R_{cl}\,|\,N=%d$' % tot_stars, 
             transform = ax10.transAxes, 
             bbox=dict(facecolor='white', alpha=0.5), fontsize=16)
    # Plot stars.
    stars_rjct_temp = [[], []]
    for star in stars_in_rjct:
        stars_rjct_temp[0].append(star[5])
        stars_rjct_temp[1].append(star[3])
    plt.scatter(stars_rjct_temp[0], stars_rjct_temp[1], marker='x', c='teal', 
                s=12, zorder=1)
    stars_acpt_temp = [[], []]
    for star in stars_in:
        stars_acpt_temp[0].append(star[5])
        stars_acpt_temp[1].append(star[3])
    sz_pt = 0.5 if (len(stars_in_rjct)+len(stars_in)) > 1000 else 1.
    plt.scatter(stars_acpt_temp[0], stars_acpt_temp[1], marker='o', c='k', 
                s=sz_pt, zorder=2)
        

    # T1 normalized histogram of stars outside cluster.
    ax11 = plt.subplot(gs1[4, 6:8])
    #Set plot limits
    x_min, x_max = min(mag_data)-0.5, max(mag_data)+0.5
    plt.xlim(x_min, x_max)
    plt.text(0.05, 0.8, r'$r > R_{cl}$', transform = ax11.transAxes, 
             bbox=dict(facecolor='white', alpha=0.5), fontsize=14)
    ax11.minorticks_on()
    #Set axis labels
    plt.xlabel(r'$T_1$', fontsize=18)
    # Plot stars.
    stars_out_temp = []
    for star in stars_out_rjct:
        stars_out_temp.append(star[3])
    for star in stars_out:
        stars_out_temp.append(star[3])
    binwidth = 0.25
    plt.hist(stars_out_temp, bins=np.arange(int(x_min), int(x_max+binwidth), 
                                        binwidth), normed=True)


    # T1 normalized histogram of stars inside cluster.
    ax12 = plt.subplot(gs1[5, 6:8])
    #Set plot limits
    x_min, x_max = min(mag_data)-0.5, max(mag_data)+0.5
    plt.xlim(x_min, x_max)
    plt.text(0.05, 0.8, r'$r \leq R_{cl}$', transform = ax12.transAxes, 
             bbox=dict(facecolor='white', alpha=0.5), fontsize=14)
    ax12.minorticks_on()
    #Set axis labels
    plt.xlabel(r'$T_1$', fontsize=18)
    # Plot stars.
    stars_in_temp = []
    for star in stars_in_rjct:
        stars_in_temp.append(star[3])
    for star in stars_in:
        stars_in_temp.append(star[3])
    binwidth = 0.25
    plt.hist(stars_in_temp, bins=np.arange(int(x_min), int(x_max+binwidth), 
                                        binwidth), normed=True)



    # Star regions as defined by the decontamination algorithm.
    # Check if decont algorithm was applied.
    if not(flag_area_stronger):
        ax13 = plt.subplot(gs1[6:8, 0:2])
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
        ax13.minorticks_on()
        ax13.grid(b=True, which='both', color='gray', linestyle='--', lw=0.5)
        # Radius
        circle = plt.Circle((center_cl[0], center_cl[1]), clust_rad[0], 
                            color='r', fill=False)
        fig.gca().add_artist(circle)
        plt.text(0.4, 0.92, 'Cluster + %d Field regions' % (len(field_region)), 
                 transform = ax13.transAxes, 
                 bbox=dict(facecolor='white', alpha=0.8), fontsize=12)
        # Plot cluster region.
        clust_reg_temp = [[], []]
        for star in cluster_region:
            dist = np.sqrt((center_cl[0]-star[1])**2 + \
            (center_cl[1]-star[2])**2)
            # Only plot stars inside the cluster's radius.
            if dist <= clust_rad[0]:
                clust_reg_temp[0].append(star[1])
                clust_reg_temp[1].append(star[2])
        plt.scatter(clust_reg_temp[0], clust_reg_temp[1], marker='o', c='black',
                    s=8, edgecolors='none')
        # Plot field stars regions.
        col = cycle(['red', 'darkgreen', 'blue', 'maroon'])
#        col = iter(plt.cm.rainbow(np.linspace(0, 1, len(field_region))))
        for i, reg in enumerate(field_region):
            stars_reg_temp = [[], []]
            for star in reg:
                # star[1] is the x coordinate and star[2] the y coordinate
                stars_reg_temp[0].append(star[1])
                stars_reg_temp[1].append(star[2])
            plt.scatter(stars_reg_temp[0], stars_reg_temp[1], marker='o', 
                        c=next(col), s=8, edgecolors='none')
   

    # Stars in the fields region with their boxes/KDE, depending on the
    # algorithm used.
    # Check if decont algorithm was applied.
    if not(flag_area_stronger):
        # Plot first field region.
        region = 0
        ax14 = plt.subplot(gs1[6:8, 2:4])
        plt.xlabel(r'$C-T_1$', fontsize=18)
        plt.ylabel(r'$T_1$', fontsize=18)
        col1_acpt_out, mag_acpt_out = [], []
        for star in stars_out:
            col1_acpt_out.append(star[5])
            mag_acpt_out.append(star[3])
        #Set plot limits
        plt.xlim(col1_min, col1_max)
        plt.ylim(mag_min, mag_max)
        # Set minor ticks
        ax14.minorticks_on()
        ax14.xaxis.set_major_locator(MultipleLocator(1.0))
        ax14.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
        plt.text(0.63, 0.93, 'Field region %d' % region, transform = \
        ax14.transAxes, bbox=dict(facecolor='white', alpha=0.8), fontsize=12)
        
        # Plot field star region for variable box algorithm.
#        boxes = plt.gca()
#        col = ['red', 'darkgreen', 'blue', 'maroon']
#        stars_reg_temp = [[], []]
#        for star in field_reg_box[region]:
#            # star[2] is the box size in T1 and star[3] in CT1.
#            boxes.add_patch(Rectangle(((star[1]-star[3]), \
#            (star[0]-star[2])), star[3]*2, star[2]*2, facecolor='none', \
#            edgecolor=col[0], lw=0.8))
#            # star[0] is the T1 coordinate and star[1] the CT1 coordinate.
#            stars_reg_temp[0].append(star[1])
#            stars_reg_temp[1].append(star[0])
#        plt.scatter(stars_reg_temp[0], stars_reg_temp[1], marker='o', 
#                    c='black', s=5, edgecolors='none', zorder=2)

        # Plot field stars and KDE for KDE algorithm.
        stars_reg_temp = [[], []]
        for star in field_region[region]:
            # star[3] is the T1 coordinate and star[5] the CT1 coordinate.
            stars_reg_temp[0].append(star[5])
            stars_reg_temp[1].append(star[3])
        plt.scatter(stars_reg_temp[0], stars_reg_temp[1], marker='o', 
                    c='black', s=4., edgecolors='none', zorder=2)
        # Plot field KDE.
        plt.imshow(np.rot90(kde_f), cmap=plt.cm.gist_earth_r, 
                   extent=[col1_min, col1_max, mag_min, mag_max],\
                   aspect='auto')
        # Plot ZAMS.
#        plt.plot(zams_iso[1], zams_iso[0], c='k', ls='--', lw=1.)



    # Cluster's stars CMD (stars inside cluster's radius) cleaned by the decont
    # algorithm.
    # Check if decont algorithm was applied.
    if not(flag_area_stronger):
        ax15 = plt.subplot(gs1[6:8, 4:6])
        #Set plot limits
        plt.xlim(col1_min, col1_max)
        plt.ylim(mag_min, mag_max)
        #Set axis labels
        plt.xlabel(r'$C-T_1$', fontsize=18)
        plt.ylabel(r'$T_1$', fontsize=18)
        # Add text box
        text1 = r'$E_{(B-V)} = %0.2f}$' '\n' % cl_e_bv
        text2 = r'$Age = %0.3f}$' '\n' % cl_age
        text3 = r'$[Fe/H] = %0.2f}$' '\n' % cl_feh
        text4 = r'$(m-M)_o = %0.2f}$' % cl_dmod
        text = text1+text2+text3+text4
        plt.text(0.05, 0.8, text, transform = ax15.transAxes, 
                 bbox=dict(facecolor='white', alpha=0.6), fontsize=12)
        tot_kde_clust = len(memb_prob_avrg_sort)
        plt.text(0.55, 0.93, r'$r \leq R_{cl}\,|\,N=%d$' % tot_kde_clust,
                 transform=ax15.transAxes, \
                 bbox=dict(facecolor='white', alpha=0.5), fontsize=16)
        # Set minor ticks
        ax15.minorticks_on()
        ax15.xaxis.set_major_locator(MultipleLocator(1.0))
        ax15.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
        stars_clust_temp = [[], [], []]
        # Plot stars inside cluster radius.
        for star in memb_prob_avrg_sort:
            stars_clust_temp[0].append(star[5])
            stars_clust_temp[1].append(star[3])
        plt.scatter(stars_clust_temp[0], stars_clust_temp[1], marker='o', 
                    c='black', s=4., edgecolors='none')
#        # Plot ZAMS.
#        plt.plot(zams_iso[1], zams_iso[0], c='k', ls='--', lw=1.)
#        # Plot isochrone.
#        plt.plot(iso_moved[1], iso_moved[0], c='crimson', ls='--', lw=1.5)
        # Plot KDE. Uncomment if KDE algorithm is used.
        plt.imshow(np.rot90(kde_cl), cmap=plt.cm.gist_earth_r, 
                   extent=[col1_min, col1_max, mag_min, mag_max], aspect='auto')



    # Finding chart of cluster region with decontamination algorithm applied.
    # Used for the finding chart with colors assigned according to the
    # probabilities obtained.
    # Check if decont algorithm was applied.
    if not(flag_area_stronger):
        ax16 = plt.subplot(gs1[6:8, 6:8])
        # Get max and min values in x,y
        x_min, x_max = 10000., -10000
        y_min, y_max = 10000., -10000
        for star in cluster_region:
            x_min, x_max = min(star[1], x_min), max(star[1], x_max)
            y_min, y_max = min(star[2], y_min), max(star[2], y_max)
        #Set plot limits
        plt.xlim(x_min, x_max)
        plt.ylim(y_min, y_max)
        #Set axis labels
        plt.xlabel('x (px)', fontsize=12)
        plt.ylabel('y (px)', fontsize=12)
        # Set minor ticks
        ax16.minorticks_on()
        # Radius
        circle = plt.Circle((center_cl[0], center_cl[1]), clust_rad[0], 
                            color='red', fill=False)
        fig.gca().add_artist(circle)
        plt.text(0.63, 0.93, 'Cluster region', transform = ax16.transAxes, \
        bbox=dict(facecolor='white', alpha=0.8), fontsize=12)
        # This color scheme makes higher prob stars look blacker.
        cm = plt.cm.get_cmap('gist_yarg')
        # Plot cluster region.
        stars_clust_temp = [[], [], []]
        for st_indx, star in enumerate(cluster_region):
            if backg_value[0] > 0.005:
                # Dense field
                star_size = 15
            else:
                star_size = 20
            stars_clust_temp[0].append(star[1])
            stars_clust_temp[1].append(star[2])
            stars_clust_temp[2].append(200*clust_reg_prob_avrg[st_indx]**8)       
        plt.scatter(stars_clust_temp[0], stars_clust_temp[1], marker='o', 
                    c=stars_clust_temp[2], s=star_size, edgecolors='black',\
                    cmap=cm)
            

    # Integrated magnitude distribution.
    ax17 = plt.subplot(gs1[8:10, 0:2])
    plt.xlim(min(stars_in_all_mag[0])-0.2, max(stars_in_all_mag[0])+0.2)
    plt.ylim(max(max(stars_in_all_mag[1]), max(stars_in_mag[1]))+0.2,
             min(stars_in_all_mag[1])-0.2)
    plt.xlabel(r'$T_1$', fontsize=18)
    plt.ylabel(r'$T_1^*$', fontsize=18)
    ax17.minorticks_on()
    ax17.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
    # Text.
    text1 = r'$stars_{acpt}\,(T_{1,max}^*\,=%.2f$)' '\n' % (min(stars_in_mag[1]))
    text2 =r'$stars_{all}\;(T_{1,max}^*\,=%.2f)$' % (min(stars_in_all_mag[1]))
    # Only accepted stars.
    plt.plot(stars_in_mag[0], stars_in_mag[1], 'r-', lw=1.5, label=text1)
    # All stars, including those rejected due to its large errors.
    plt.plot(stars_in_all_mag[0], stars_in_all_mag[1], 'b--', lw=1., label=text2)
    # get handles
    handles, labels = ax17.get_legend_handles_labels()
    # use them in the legend
    leg = ax17.legend(handles, labels, loc='lower right', numpoints=1,
                      fontsize=12)
    leg.get_frame().set_alpha(0.5)
    
    

    # Norm fit for KDE probability values.
    # Check if decont algorithm was applied.
    if not(flag_area_stronger):
        ax18 = plt.subplot(gs1[8:10, 2:4])
        plt.xlim(0., 1.)
        plt.xlabel('cluster membership prob', fontsize=12)
        ax18.minorticks_on()
        ax18.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
        prob_data = [star[7] for star in memb_prob_avrg_sort]
        # Best Gaussian fit of data.
        (mu, sigma) = norm.fit(prob_data)
        # Text.
        text = r'$\mu=%.3f,\ \sigma=%.3f$' %(mu, sigma)
        plt.text(0.05, 0.92, text, transform = ax18.transAxes, 
             bbox=dict(facecolor='white', alpha=0.85), fontsize=12)
        # Histogram of the data.
        n, bins, patches = plt.hist(prob_data, 60, normed=1, facecolor='green',
                                    alpha=0.75)
        # Best fit line.
        y = mlab.normpdf( bins, mu, sigma)
        plt.plot(bins, y, 'r--', linewidth=2)



    # Cluster's stars CMD. Check if decont algorithm was applied.
    if not(flag_area_stronger):
        ax19 = plt.subplot(gs1[8:10, 4:6])
        #Set plot limits
        plt.xlim(col1_min, col1_max)
        plt.ylim(mag_min, mag_max)
        #Set axis labels
        plt.xlabel(r'$C-T_1$', fontsize=18)
        plt.ylabel(r'$T_1$', fontsize=18)
        tot_kde_clust = len(memb_prob_avrg_sort)
        text = r'$N=%d\,|\,MI \geq 0.$' % tot_kde_clust
        plt.text(0.05, 0.93, text, transform = ax19.transAxes,
                 bbox=dict(facecolor='white', alpha=0.5), fontsize=14)
        # Set minor ticks
        ax19.minorticks_on()
        ax19.xaxis.set_major_locator(MultipleLocator(1.0))
        ax19.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
        # Plot ZAMS.
        plt.plot(zams_iso[1], zams_iso[0], c='k', ls='--', lw=1.)
        # Plot isochrone.
        plt.plot(iso_moved[1], iso_moved[0], 'g', lw=1.2)
        # This reversed colormap means higher prob stars will look redder.
        cm = plt.cm.get_cmap('RdYlBu_r')
        m_p_m_temp = [[], [], []]
        for star in memb_prob_avrg_sort:
            m_p_m_temp[0].append(star[5])
            m_p_m_temp[1].append(star[3])
            m_p_m_temp[2].append(star[7])
        # Create new list with inverted values so higher prob stars are on top.
        m_p_m_temp_inv = [i[::-1] for i in m_p_m_temp]
        plt.scatter(m_p_m_temp_inv[0], m_p_m_temp_inv[1], marker='o', 
                    c=m_p_m_temp_inv[2], s=40, cmap=cm, lw=0.5, vmin=0, vmax=1)
        # If list is not empty.
        if m_p_m_temp_inv[1]:
            # Plot error bars at several mag values.
            mag_y = np.arange(int(min(m_p_m_temp_inv[1])+0.5), 
                              int(max(m_p_m_temp_inv[1])+0.5) + 0.1)
            x_val = [min(3.9, max(col1_data)+0.2) - 0.4]*len(mag_y)
            plt.errorbar(x_val, mag_y, yerr=func(mag_y, *popt_mag), 
                         xerr=func(mag_y, *popt_col1), fmt='k.', lw=0.8, ms=0.,\
                         zorder=4)
            # Plot colorbar.
            cbaxes = fig.add_axes([0.678, 0.31, 0.04, 0.005])
            cbar = plt.colorbar(cax=cbaxes, ticks=[0,1], orientation='horizontal')
            cbar.ax.tick_params(labelsize=9)



    # Cluster's stars CMD. Check if decont algorithm was applied.
    if not(flag_area_stronger):
        ax20 = plt.subplot(gs1[8:10, 6:8])
        #Set plot limits
        plt.xlim(col1_min, col1_max)
        plt.ylim(mag_min, mag_max)
        #Set axis labels
        plt.xlabel(r'$C-T_1$', fontsize=18)
        plt.ylabel(r'$T_1$', fontsize=18)
        # Set minor ticks
        ax20.minorticks_on()
        ax20.xaxis.set_major_locator(MultipleLocator(1.0))
        ax20.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
        # Plot ZAMS.
        plt.plot(zams_iso[1], zams_iso[0], c='k', ls='--', lw=1.)
        # Plot isochrone.
        plt.plot(iso_moved[1], iso_moved[0], 'g', lw=1.2)
        # This reversed colormap means higher prob stars will look redder.
        cm = plt.cm.get_cmap('RdYlBu_r')
        m_p_m_temp = [[], [], []]
        for star in memb_prob_avrg_sort:
            # Only plot stars with MI>=mu
            if star[7] >= mu:
                m_p_m_temp[0].append(star[5])
                m_p_m_temp[1].append(star[3])
                m_p_m_temp[2].append(star[7])
        # Create new list with inverted values so higher prob stars are on top.
        m_p_m_temp_inv = [i[::-1] for i in m_p_m_temp]
        plt.text(0.05, 0.93, r'$N=%d\,|\,MI \geq \mu$' % len(m_p_m_temp[0]), 
                 transform = ax20.transAxes, 
                 bbox=dict(facecolor='white', alpha=0.5), fontsize=14)
        v_min = round(mu, 2)
        v_max = round(max(m_p_m_temp[2]),2) if m_p_m_temp[2] else 1.
        plt.scatter(m_p_m_temp_inv[0], m_p_m_temp_inv[1], marker='o', 
                    c=m_p_m_temp_inv[2], s=40, cmap=cm, lw=0.5, vmin=v_min, \
                    vmax=v_max)
        # If list is not empty.
        if m_p_m_temp_inv[1]:
            # Plot error bars at several mag values.
            mag_y = np.arange(int(min(m_p_m_temp_inv[1])+0.5), 
                              int(max(m_p_m_temp_inv[1])+0.5) + 0.1)
            x_val = [min(3.9, max(col1_data)+0.2) - 0.4]*len(mag_y)
            plt.errorbar(x_val, mag_y, yerr=func(mag_y, *popt_mag), 
                         xerr=func(mag_y, *popt_col1), fmt='k.', lw=0.8, ms=0.,\
                         zorder=4)
            # Plot colorbar.
            cbaxes20 = fig.add_axes([0.93, 0.31, 0.04, 0.005])
            cbar20 = plt.colorbar(cax=cbaxes20, ticks=[v_min,v_max],
                                  orientation='horizontal')
            cbar20.ax.tick_params(labelsize=9) 
        
        
        
    # Distribution of p_values.
    # pval_test_params[-1] is the flag that tells me if the block was processed.
    if not(flag_area_stronger) and pval_test_params[-1]:
        # Extract parameters from list.
        prob_cl_kde, p_vals_cl, p_vals_f, kde_cl_1d, kde_f_1d, x_kde, y_over = pval_test_params[:-1]
        ax21 = plt.subplot(gs1[10:12, 0:2])
        plt.xlim(-0.5, 1.5)
        plt.ylim(0, max(max(kde_f_1d), max(kde_cl_1d))+0.5)
        plt.xlabel('p-values', fontsize=12)
        ax21.minorticks_on()
        ax21.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
        # Grid to background.
        ax21.set_axisbelow(True)
        # Plot cluster vs field KDE.
        plt.plot(x_kde, kde_cl_1d, c='b', ls='-', lw =1., label=r'$KDE_{cl}$')
        # Plot field vs field KDE.
        plt.plot(x_kde, kde_f_1d, c='r', ls='-', lw=1., label=r'$KDE_{f}$')
        # Fill overlap.
        plt.fill_between(x_kde, y_over, 0, color='grey', alpha='0.5')
        text = r'$\;P_{cl}^{KDE} = %0.2f$' % round(prob_cl_kde,2)
        plt.text(0.05, 0.92, text, transform = ax21.transAxes, 
             bbox=dict(facecolor='white', alpha=0.6), fontsize=12)
        handles, labels = ax21.get_legend_handles_labels()
        leg = ax21.legend(handles, labels, loc='upper right', numpoints=1,
                          fontsize=12)
        leg.get_frame().set_alpha(0.6)
        
        
        # QQ-plot.
        # Extract parameters from list.
        ccc, quantiles, r_squared, slope, intercept = qq_params
        ax22 = plt.subplot(gs1[10:12, 2:4])
        plt.xlim(-0.05, 1.05)
        plt.ylim(-0.05, 1.05)
        plt.xlabel(r'$p-value_{cl}$', fontsize=16)
        plt.ylabel(r'$p-value_{f}$', fontsize=16)
        ax22.minorticks_on()
        ax22.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
        plt.scatter(quantiles[0], quantiles[1], marker='o', c='k', s=10.)
        text1 = r'$R^2\, = %0.2f$' '\n' % r_squared
        text2 = r'$CCC\, = %0.2f$' % ccc
        text = text1+text2
        plt.text(0.72, 0.87, text, transform = ax22.transAxes, 
             bbox=dict(facecolor='white', alpha=0.85), fontsize=12)
        plt.plot([0., 1.], [0., 1.], color='k', linestyle='--', linewidth=1.)
        # Plot qq-plot fitted line.
        x = np.arange(0.,1.1,0.1)
        y = line(x, slope, intercept)
        plt.plot(x,y, color='r', linestyle='--')
           
        
        





    fig.tight_layout()

    # Generate output file for each data file.
    plt.savefig(join(output_subdir, str(clust_name)+'.png'), dpi=150)
    
    # Close to release memory.
    plt.clf()
    plt.close()
    

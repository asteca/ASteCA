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
from scipy.ndimage.filters import gaussian_filter
from matplotlib.patches import Ellipse

from os.path import join


def make_plots(output_subdir, clust_name, x_data, y_data, center_cl, 
               cent_cl_err, x_center_bin, y_center_bin, h_filter, radii, 
               backg_value, inner_ring, outer_ring, radius_params,
               ring_density, poisson_error, cont_index, width_bins,
               mag_data, col1_data, popt_mag, popt_col1,
               err_plot, rjct_errors_fit, k_prof, k_pr_err,
               flag_king_no_conver, stars_in,
               stars_out, stars_in_rjct, stars_out_rjct, stars_in_mag,
               stars_in_all_mag, n_c, flag_area_stronger,
               cluster_region, field_region, pval_test_params, qq_params,
               clust_reg_prob_avrg, memb_prob_avrg_sort, bf_params, bf_return,
               ga_params, er_params, axes_params, ps_params):
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
        a, b, c, d = k_prof[0], k_prof[1], k_prof[2], backg_value
        return c*(1/np.sqrt(1+(x/a)**2) - 1/np.sqrt(1+(b/a)**2))**2 + d
        
    # Name for axes.
    x_ax, y_ax = axes_params[0], axes_params[1]
    
    # Define plot limits for *all* CMD diagrams.
    x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd = axes_params[2]
    col1_min, col1_max = max(x_min_cmd, min(col1_data)-0.2),\
    min(x_max_cmd, max(col1_data)+0.2)
    mag_min, mag_max = min(y_max_cmd, max(mag_data)+0.5),\
    max(y_min_cmd, min(mag_data)-0.5)

    # Parameters from get_radius function.
    clust_rad, delta_backg, delta_percentage = radius_params
   
    # Error parameters.
    be, be_e, e_max = er_params   
   
    # Parameters from error fitting.
    bright_end, popt_umag, pol_mag, popt_ucol1, pol_col1, mag_val_left,\
    mag_val_right, col1_val_left, col1_val_right = err_plot
    
    # Best isochrone fit params.
    bf_flag, best_fit_algor, N_b = bf_params
    
    # Genetic algorithm params.
    n_pop, n_gen, fdif, p_cross, cr_sel, p_mut, n_el, n_ei, n_es = ga_params
    
    # Best fitting process results.
    isoch_fit_params, isoch_fit_errors, shift_isoch, synth_clst = bf_return
    
    # Plot all outputs
    # figsize(x1, y1), GridSpec(y2, x2) --> To have square plots: x1/x2 = 
    # y1/y2 = 2.5 
    fig = plt.figure(figsize=(20, 35)) # create the top-level container
    gs1 = gridspec.GridSpec(14, 8)  # create a GridSpec object
    #gs1.update(wspace=.09, hspace=.0)

    # 2D filtered histogram, largest bin width.
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
        
    # 2D filtered histogram.
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
        
    # 2D filtered histogram.
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
        
    # 2D filtered histogram, smallest bin width.
    ax3 = plt.subplot(gs1[0:2, 6:8])
    plt.xlabel('x (bins)', fontsize=12)
    plt.ylabel('y (bins)', fontsize=12)
    ax3.minorticks_on()
    plt.axvline(x=x_center_bin[0], linestyle='--', color='white')
    plt.axhline(y=y_center_bin[0], linestyle='--', color='white')
    # Radius
    circle = plt.Circle((x_center_bin[0], y_center_bin[0]), clust_rad/25., 
                        color='w', fill=False)
    fig.gca().add_artist(circle)
    #text = 'Center (%d, %d)' % (x_center_bin[0], y_center_bin[0])
    text = 'Bin: %d px' % (width_bins[0])
    plt.text(0.7, 0.92, text, transform = ax3.transAxes, 
             bbox=dict(facecolor='white', alpha=0.8), fontsize=12)
    plt.imshow(h_filter[0].transpose(), origin='lower')

    # x,y finding chart of full frame
    ax4 = plt.subplot(gs1[2:4, 0:2])
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
    ax4.minorticks_on()
    circle = plt.Circle((center_cl[0], center_cl[1]), clust_rad, color='r', 
                        fill=False)
    fig.gca().add_artist(circle)
    # Add text box
    text1 = '$x_{cent} = %d \pm %d px$' % (center_cl[0], cent_cl_err[0])
    text2 = '\n'
    text3 = '$y_{cent} = %d \pm %d px$' % (center_cl[1], cent_cl_err[1])
    text4 = text1+text2+text3
    plt.text(0.5, 0.85, text4, transform = ax4.transAxes, \
    bbox=dict(facecolor='white', alpha=0.85), fontsize=15)
    # Count the number of very bright stars in the field.
    range_10_perc = (max(mag_data)-min(mag_data))/10.+min(mag_data)
    bright_stars = len([i for i in mag_data if i < range_10_perc])
    # Set exponential factor for high and low density fields.
    exp_factor = -0.004 if backg_value > 0.005 else -0.0035
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

    # Radial density plot.
    ax5 = plt.subplot(gs1[2:4, 2:6])
    # Get max and min values in x,y
    x_max = max(radii[0])+10
    y_min, y_max = (backg_value-delta_backg)-(max(ring_density[0])-\
    min(ring_density[0]))/10, max(ring_density[0])+(max(ring_density[0])-\
    min(ring_density[0]))/10
    # Set plot limits
    plt.xlim(-10, min(x_max, 500))
    plt.ylim(y_min, y_max)
    # Set axes labels
    plt.xlabel('radius (px)', fontsize=12)
    plt.ylabel("stars/px$^{2}$", fontsize=12)
    # Cluster's name.
    text = str(clust_name)
    plt.text(0.4, 0.9, text, transform = ax5.transAxes, fontsize=14)
    # Plot poisson error bars
    plt.errorbar(radii[0], ring_density[0], yerr = poisson_error[0], fmt='ko', 
                 zorder=1)
    ## Plot the delta around the background value used to asses when the density
    # has become stable
    plt.hlines(y=(backg_value+delta_backg), xmin=0, xmax=max(radii[0]), 
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
             '$R_{cl}$ = %d $\pm$ %d px' % (clust_rad,
                                           round(width_bins[0]/2.)),\
            '3-P King profile',
             '$R_c$ = %d $\pm$ %d px' % (k_prof[0], rc_err),\
             '$R_t$ = %d $\pm$ %d px' % (k_prof[1], rt_err),\
             'Background ($\Delta$=%d%%)' % delta_percentage]
    # Plot density profile with the smallest bin size
    ax5.plot(radii[0], ring_density[0], 'ko-', zorder=3, label=texts[0])
    # Plot King profile.
    arr_y_up = (y_max-y_min)/2. + y_min # Middle of the graph.
    head_l = abs((arr_y_up-backg_value)/7.) # Length of arrow head.
    arr_y_dwn = -1.*abs(arr_y_up-backg_value-head_l*1.5) # Length of arrow.
    if flag_king_no_conver == False:
        ax5.plot(radii[0], three_params(radii[0]), 'b--', label=texts[2],
                 lw=1.5, zorder=3)
        # Plot R_c as a dashed line.
        ax5.vlines(x=k_prof[0], ymin=0, ymax=max(ring_density[0])/1.2,
                   label=texts[3], color='g', linestyles='dashed', zorder=4)
        # Plot R_t radius as an arrow. vline is there to show the label.
        ax5.vlines(x=k_prof[1], ymin=0., ymax=0., label=texts[4], color='b')
        ax5.arrow(k_prof[1], arr_y_up, 0., arr_y_dwn, fc="b", ec="b",\
                  head_width=10., head_length=head_l, zorder=5)
    # Plot radius.
    ax5.vlines(x=clust_rad, ymin=0, ymax=0., label=texts[1], color='r')
    ax5.arrow(clust_rad, arr_y_up, 0., arr_y_dwn, fc="r",
              ec="r", head_width=10., head_length=head_l, zorder=5)
    # Plot background level.
    ax5.hlines(y=backg_value, xmin=0, xmax=max(radii[0]), 
               label=texts[5], color='k', zorder=5)
    # get handles
    handles, labels = ax5.get_legend_handles_labels()
    # use them in the legend
    ax5.legend(handles, labels, loc='upper right', numpoints=1, fontsize=10)
    ax5.minorticks_on()

    # Zoom on x,y finding chart
    ax6 = plt.subplot(gs1[2:4, 6:8])
    #Set plot limits
    x_min, x_max = min(x_data), max(x_data)
    y_min, y_max = min(y_data), max(y_data)
    # If possible, zoom in.
    x_min, x_max = max(x_min, (center_cl[0]-1.5*clust_rad)), \
    min(x_max, (center_cl[0]+1.5*clust_rad))
    y_min, y_max = max(y_min, (center_cl[1]-1.5*clust_rad)), \
    min(y_max, (center_cl[1]+1.5*clust_rad))
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    #Set axis labels
    plt.xlabel('x (px)', fontsize=12)
    plt.ylabel('y (px)', fontsize=12)
    # Set minor ticks
    ax6.minorticks_on()
    # Add circle
    circle = plt.Circle((center_cl[0], center_cl[1]), clust_rad, color='r', 
                        fill=False)
    fig.gca().add_artist(circle)
    text1 = 'Cluster (zoom)\n'
    text2 = 'CI = %0.2f' % (cont_index)
    text = text1 + text2
    plt.text(0.6, 0.9, text, transform = ax6.transAxes, 
             bbox=dict(facecolor='white', alpha=0.85), fontsize=12)
    if backg_value > 0.005:
        plt.scatter(x_data, y_data, marker='o', c='black', 
                    s=500*np.exp(-0.003*mag_data**2.5))
    else:
        plt.scatter(x_data, y_data, marker='o', c='black', 
                    s=500*np.exp(-0.0025*mag_data**2.5))
        
    # Cluster and field regions defined.
    ax7 = plt.subplot(gs1[4:6, 0:2])
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
    ax7.minorticks_on()
    ax7.grid(b=True, which='both', color='gray', linestyle='--', lw=0.5)
    # Radius
    circle = plt.Circle((center_cl[0], center_cl[1]), clust_rad, 
                        color='r', fill=False)
    fig.gca().add_artist(circle)
    plt.text(0.4, 0.92, 'Cluster + %d Field regions' % (len(field_region)), 
             transform = ax7.transAxes, 
             bbox=dict(facecolor='white', alpha=0.8), fontsize=12)
    # Plot cluster region.
    clust_reg_temp = [[], []]
    for star in cluster_region:
        dist = np.sqrt((center_cl[0]-star[1])**2 + \
        (center_cl[1]-star[2])**2)
        # Only plot stars inside the cluster's radius.
        if dist <= clust_rad:
            clust_reg_temp[0].append(star[1])
            clust_reg_temp[1].append(star[2])
    plt.scatter(clust_reg_temp[0], clust_reg_temp[1], marker='o', c='black',
                s=8, edgecolors='none')
    if not(flag_area_stronger):
        # Plot field stars regions.
        col = cycle(['red', 'darkgreen', 'blue', 'maroon'])
        for i, reg in enumerate(field_region):
            stars_reg_temp = [[], []]
            for star in reg:
                # star[1] is the x coordinate and star[2] the y coordinate
                stars_reg_temp[0].append(star[1])
                stars_reg_temp[1].append(star[2])
            plt.scatter(stars_reg_temp[0], stars_reg_temp[1], marker='o', 
                        c=next(col), s=8, edgecolors='none')        
        
    # Field stars CMD (stars outside cluster's radius)
    ax8 = plt.subplot(gs1[4:6, 2:4])
    #Set plot limits
    plt.xlim(col1_min, col1_max)
    plt.ylim(mag_min, mag_max)
    #Set axis labels
    plt.xlabel('$'+x_ax+'$', fontsize=18)
    plt.ylabel('$'+y_ax+'$', fontsize=18)
    # Set minor ticks
    ax8.minorticks_on()
    # Only draw units on axis (ie: 1, 2, 3)
    ax8.xaxis.set_major_locator(MultipleLocator(1.0))
    # Set grid
    ax8.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
    plt.text(0.8, 0.93, '$r > R_{cl}$', transform = ax8.transAxes, \
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
    ax9 = plt.subplot(gs1[4:6, 4:6])
    #Set plot limits
    plt.xlim(col1_min, col1_max)
    plt.ylim(mag_min, mag_max)
    #Set axis labels
    plt.xlabel('$'+x_ax+'$', fontsize=18)
    plt.ylabel('$'+y_ax+'$', fontsize=18)
    # Set minor ticks
    ax9.minorticks_on()
    # Only draw units on axis (ie: 1, 2, 3)
    ax9.xaxis.set_major_locator(MultipleLocator(1.0))
    # Set grid
    ax9.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
    # Calculate total number of stars whitin cluster's radius.
    tot_stars = len(stars_in_rjct) + len(stars_in)
    plt.text(0.55, 0.93, '$r \leq R_{cl}\,|\,N=%d$' % tot_stars, 
             transform = ax9.transAxes, 
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
        
    # T1 magnitude error
    ax10 = plt.subplot(gs1[4, 6:8])
    #Set plot limits
    plt.xlim(min(mag_data)-0.5, max(mag_data)+0.5)
    plt.ylim(-0.005, e_max+0.01)
    #Set axis labels
    plt.ylabel('$\sigma_{'+y_ax+'}$', fontsize=18)
    plt.xlabel('$'+y_ax+'$', fontsize=18)
    # Set minor ticks
    ax10.minorticks_on()
    mag_x = np.linspace(min(mag_data), max(mag_data), 50)
    # Condition to not plot the lines if the fit was rejected.
    # Plot lower envelope.
    ax10.plot(mag_x, func(mag_x, *popt_mag), 'r-', zorder=3)
    if not rjct_errors_fit:
        # Plot left side of upper envelope (exponential).
        ax10.plot(mag_val_left, func(mag_val_left, *popt_umag), 'r--', lw=2.,
                 zorder=3)
        # Plot right side of upper envelope (polynomial).
        ax10.plot(mag_val_right, np.polyval(pol_mag, (mag_val_right)),
                 'r--', lw=2., zorder=3)
    # Plot rectangle.
    ax10.vlines(x=bright_end+0.05, ymin=-0.005, ymax=be_e, color='r', 
               linestyles='dashed', zorder=2)
    ax10.vlines(x=min(mag_data)-0.05, ymin=-0.005, ymax=be_e, color='r',
               linestyles='dashed', zorder=2)
    ax10.hlines(y=be_e, xmin=min(mag_data), xmax=bright_end, color='r',
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
    ax11 = plt.subplot(gs1[5, 6:8])
    #Set plot limits
    plt.xlim(min(mag_data)-0.5, max(mag_data)+0.5)
    plt.ylim(-0.005, e_max+0.01)
    #Set axis labels
    plt.ylabel('$\sigma_{'+x_ax+'}$', fontsize=18)
    plt.xlabel('$'+y_ax+'$', fontsize=18)
    # Set minor ticks
    ax11.minorticks_on()
    # Condition to not plot the lines if the fit was rejected.
    # Plot lower envelope.
    ax11.plot(mag_x, func(mag_x, *popt_col1), 'r-', zorder=3)
    if not rjct_errors_fit:
        # Plot left side of upper envelope (exponential).
        ax11.plot(col1_val_left, func(col1_val_left, *popt_ucol1), 'r--', lw=2.,
                 zorder=3)
        # Plot right side of upper envelope (polynomial).
        ax11.plot(col1_val_right, np.polyval(pol_col1, (col1_val_right)),
                 'r--', lw=2., zorder=3)
    # Plot rectangle.
    ax11.vlines(x=bright_end+0.05, ymin=-0.005, ymax=be_e, color='r', 
               linestyles='dashed', zorder=2)
    ax11.vlines(x=min(mag_data)-0.05, ymin=-0.005, ymax=be_e, color='r',
               linestyles='dashed', zorder=2)
    ax11.hlines(y=be_e, xmin=min(mag_data), xmax=bright_end, color='r',
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

    # LF of stars in cluster region and outside.
    ax12 = plt.subplot(gs1[6:8, 0:2])
    #Set plot limits
    x_min, x_max = min(mag_data)-0.5, max(mag_data)+0.5
    plt.xlim(x_max, x_min)
    ax12.minorticks_on()
    # Only draw units on axis (ie: 1, 2, 3)
    ax12.xaxis.set_major_locator(MultipleLocator(2.0))
    # Set grid
    ax12.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
    #Set axis labels
    plt.xlabel('$'+y_ax+'$', fontsize=18)
    plt.ylabel('$log(N^{\star})$', fontsize=18)
    # Plot create lists with mag values.
    stars_out_temp = []
    for star in stars_out_rjct:
        stars_out_temp.append(star[3])
    for star in stars_out:
        stars_out_temp.append(star[3])
    stars_in_temp = []
    for star in stars_in_rjct:
        stars_in_temp.append(star[3])
    for star in stars_in:
        stars_in_temp.append(star[3])        
    # Plot histograms.
    binwidth = 0.25
    plt.hist(stars_out_temp,
             bins=np.arange(int(x_min), int(x_max+binwidth), binwidth),\
             log=True, histtype='step', label='$r > R_{cl}$', color='b')
    plt.hist(stars_in_temp,
             bins=np.arange(int(x_min), int(x_max+binwidth), binwidth),\
             log=True, histtype='step', label='$r \leq R_{cl}$', color='r')
    # Force y axis min to 1.
    plt.ylim(1., plt.ylim()[1])
    # Legends.
    leg11 = plt.legend(fancybox=True, loc='upper right', numpoints=1,
                       fontsize=16)
    # Set the alpha value of the legend.
    leg11.get_frame().set_alpha(0.7)

    # Integrated magnitude distribution.
    ax13 = plt.subplot(gs1[6:8, 2:4])
    plt.xlim(min(stars_in_all_mag[0])-0.2, max(stars_in_all_mag[0])+0.2)
    plt.ylim(max(max(stars_in_all_mag[1]), max(stars_in_mag[1]))+0.2,
             min(stars_in_all_mag[1])-0.2)
    plt.xlabel('$'+y_ax+'$', fontsize=18)
    plt.ylabel('$'+y_ax+'^*$', fontsize=18)             
    ax13.minorticks_on()
    ax13.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
    # Text.
    text1 = '$'+y_ax+'^*\,(max_{acpt})\,=%.2f$' '\n' % (min(stars_in_mag[1]))
    text2 = '$'+y_ax+'^*\,(max_{all})\;=%.2f$' % (min(stars_in_all_mag[1]))
    text = text1+text2
    plt.text(0.5, 0.25, text, transform = ax13.transAxes, 
         bbox=dict(facecolor='white', alpha=0.85), fontsize=14)
    # Only accepted stars.
    plt.plot(stars_in_mag[0], stars_in_mag[1], 'r-', lw=1.5,
             label='$stars_{acpt}$')
    # All stars, including those rejected due to its large errors.
    plt.plot(stars_in_all_mag[0], stars_in_all_mag[1], 'b--', lw=1.,
             label='$stars_{all}$')
    # get handles
    handles, labels = ax13.get_legend_handles_labels()
    # use them in the legend
    leg = ax13.legend(handles, labels, loc='lower right', numpoints=1,
                      fontsize=14)
    leg.get_frame().set_alpha(0.5)

    # Distribution of p_values.
    # pval_test_params[-1] is the flag that indicates if the block was processed.
    if pval_test_params[-1]:
        # Extract parameters from list.
        prob_cl_kde, p_vals_cl, p_vals_f, kde_cl_1d, kde_f_1d, x_kde, y_over\
        = pval_test_params[:-1]
        ax14 = plt.subplot(gs1[6:8, 4:6])
        plt.xlim(-0.5, 1.5)
        plt.ylim(0, max(max(kde_f_1d), max(kde_cl_1d))+0.5)
        plt.xlabel('p-values', fontsize=12)
        ax14.minorticks_on()
        ax14.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
        # Grid to background.
        ax14.set_axisbelow(True)
        # Plot cluster vs field KDE.
        plt.plot(x_kde, kde_cl_1d, c='b', ls='-', lw =1., label='$KDE_{cl}$')
        # Plot field vs field KDE.
        plt.plot(x_kde, kde_f_1d, c='r', ls='-', lw=1., label='$KDE_{f}$')
        # Fill overlap.
        plt.fill_between(x_kde, y_over, 0, color='grey', alpha='0.5')
        text = '$\;P_{cl}^{KDE} = %0.2f$' % round(prob_cl_kde,2)
        plt.text(0.05, 0.92, text, transform = ax14.transAxes, 
             bbox=dict(facecolor='white', alpha=0.6), fontsize=12)
        # Legend.
        handles, labels = ax14.get_legend_handles_labels()
        leg = ax14.legend(handles, labels, loc='upper right', numpoints=1,
                          fontsize=12)
        leg.get_frame().set_alpha(0.6)
        
        # QQ-plot.
        # Extract parameters from list.
        ccc, quantiles = qq_params
        ax15 = plt.subplot(gs1[6:8, 6:8])
        plt.xlim(-0.05, 1.05)
        plt.ylim(-0.05, 1.05)
        plt.xlabel('$p-value_{cl}$', fontsize=16)
        plt.ylabel('$p-value_{f}$', fontsize=16)
        ax15.minorticks_on()
        ax15.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
        text = '$CCC\, = %0.2f$' % ccc
        plt.text(0.05, 0.92, text, transform = ax15.transAxes, 
             bbox=dict(facecolor='white', alpha=0.85), fontsize=12)
        # Plot quantiles.
        plt.scatter(quantiles[0], quantiles[1], marker='o', c='k', s=10.)
        # Identity line.
        plt.plot([0., 1.], [0., 1.], color='k', linestyle='--', linewidth=1.)   

    # Stars in the first field region with their KDE.
    # Check if decont algorithm was applied.
#    if not(flag_area_stronger):
#        # Plot first field region.
#        region = 0
#        ax14 = plt.subplot(gs1[6:8, 2:4])
#        plt.xlabel('$C-T_1$', fontsize=18)
#        plt.ylabel('$T_1$', fontsize=18)
#        col1_acpt_out, mag_acpt_out = [], []
#        for star in stars_out:
#            col1_acpt_out.append(star[5])
#            mag_acpt_out.append(star[3])
#        #Set plot limits
#        plt.xlim(col1_min, col1_max)
#        plt.ylim(mag_min, mag_max)
#        # Set minor ticks
#        ax14.minorticks_on()
#        ax14.xaxis.set_major_locator(MultipleLocator(1.0))
#        ax14.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
#        plt.text(0.63, 0.93, 'Field region %d' % region, transform = \
#        ax14.transAxes, bbox=dict(facecolor='white', alpha=0.8), fontsize=12)
#        # Plot field stars and KDE for KDE algorithm.
#        stars_reg_temp = [[], []]
#        for star in field_region[region]:
#            # star[3] is the T1 coordinate and star[5] the CT1 coordinate.
#            stars_reg_temp[0].append(star[5])
#            stars_reg_temp[1].append(star[3])
#        plt.scatter(stars_reg_temp[0], stars_reg_temp[1], marker='o', 
#                    c='black', s=4., edgecolors='none', zorder=2)
#        # Plot field KDE.
#        plt.imshow(np.rot90(kde_f), cmap=plt.cm.gist_earth_r, 
#                   extent=[col1_min, col1_max, mag_min, mag_max],\
#                   aspect='auto')
        # Plot ZAMS.
#        plt.plot(zams_iso[1], zams_iso[0], c='k', ls='--', lw=1.)

    # Norm fit for decontamination algorithm probability values.
    # Check if decont algorithm was applied.
    if not(flag_area_stronger):
        ax16 = plt.subplot(gs1[8:10, 0:2])
        plt.xlim(0., 1.)
        plt.xlabel('cluster membership prob', fontsize=12)
        ax16.minorticks_on()
        ax16.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
        prob_data = [star[7] for star in memb_prob_avrg_sort]
        # Best Gaussian fit of data.
        (mu, sigma) = norm.fit(prob_data)
        # Text.
        text = '$\mu=%.3f,\ \sigma=%.3f$' %(mu, sigma)
        plt.text(0.05, 0.92, text, transform = ax16.transAxes, 
             bbox=dict(facecolor='white', alpha=0.85), fontsize=12)
        # Histogram of the data.
        n, bins, patches = plt.hist(prob_data, 50, normed=1, facecolor='green',
                                    alpha=0.75)
        # Best fit line.
        y = mlab.normpdf(bins, mu, sigma)
        plt.plot(bins, y, 'r--', linewidth=2)

        # Finding chart of cluster region with decontamination algorithm applied.
        # Used for the finding chart with colors assigned according to the
        # probabilities obtained.
        # Check if decont algorithm was applied.
        ax17 = plt.subplot(gs1[8:10, 2:4])
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
        ax17.minorticks_on()
        # Radius
        circle = plt.Circle((center_cl[0], center_cl[1]), clust_rad, 
                            color='red', fill=False)
        fig.gca().add_artist(circle)
        plt.text(0.63, 0.93, 'Cluster region', transform = ax17.transAxes, \
        bbox=dict(facecolor='white', alpha=0.8), fontsize=12)
        # Color map, higher prob stars look redder.
        cm = plt.cm.get_cmap('RdYlBu_r')
        # Star sizes for dense and not dense regions.
        star_size = 20 if backg_value > 0.005 else 35
        m_p_m_temp = [[], [], []]
        for star in memb_prob_avrg_sort:
            m_p_m_temp[0].append(star[1])
            m_p_m_temp[1].append(star[2])
            m_p_m_temp[2].append(star[7])
        # Create new list with inverted values so higher prob stars are on top.
        m_p_m_temp_inv = [i[::-1] for i in m_p_m_temp]
        plt.scatter(m_p_m_temp_inv[0], m_p_m_temp_inv[1], marker='o', 
                    c=m_p_m_temp_inv[2], s=star_size, edgecolors='black',\
                    cmap=cm, lw=0.5)
        out_clust_rad = [[], []]
        for star in cluster_region:
            dist = np.sqrt((center_cl[0]-star[1])**2 + (center_cl[1]-star[2])**2)
            # Only plot stars outside the cluster's radius.
            if dist >= clust_rad:
                out_clust_rad[0].append(star[1])
                out_clust_rad[1].append(star[2])
        plt.scatter(out_clust_rad[0], out_clust_rad[1], marker='o', 
                    s=star_size, edgecolors='black', facecolors='none', lw=0.5)

        # Star's membership probabilities on cluster's CMD.
        ax18 = plt.subplot(gs1[8:10, 4:6])
        #Set plot limits
        plt.xlim(col1_min, col1_max)
        plt.ylim(mag_min, mag_max)
        #Set axis labels
        plt.xlabel('$'+x_ax+'$', fontsize=18)
        plt.ylabel('$'+y_ax+'$', fontsize=18)
        tot_kde_clust = len(memb_prob_avrg_sort)
        text = '$r \leq R_{cl}\,|\,N=%d$' % tot_kde_clust
        plt.text(0.05, 0.93, text, transform = ax18.transAxes,
                 bbox=dict(facecolor='white', alpha=0.5), fontsize=14)
        # Set minor ticks
        ax18.minorticks_on()
        ax18.xaxis.set_major_locator(MultipleLocator(1.0))
        ax18.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
        if bf_flag:
            # Plot isochrone if best fit process was used.
            plt.plot(shift_isoch[0], shift_isoch[1], 'g', lw=1.2)
        # This reversed colormap means higher prob stars will look redder.
        cm = plt.cm.get_cmap('RdYlBu_r')
        m_p_m_temp = [[], [], []]
        for star in memb_prob_avrg_sort:
            m_p_m_temp[0].append(star[5])
            m_p_m_temp[1].append(star[3])
            m_p_m_temp[2].append(star[7])
        # Create new list with inverted values so higher prob stars are on top.
        m_p_m_temp_inv = [i[::-1] for i in m_p_m_temp]
        v_min, v_max = round(min(m_p_m_temp[2]), 2), round(max(m_p_m_temp[2]),2)
        plt.scatter(m_p_m_temp_inv[0], m_p_m_temp_inv[1], marker='o', 
                    c=m_p_m_temp_inv[2], s=40, cmap=cm, lw=0.5, vmin=v_min,
                    vmax=v_max)
        # If list is not empty.
        if m_p_m_temp_inv[1]:
            # Plot error bars at several mag values.
            mag_y = np.arange(int(min(m_p_m_temp_inv[1])+0.5), 
                              int(max(m_p_m_temp_inv[1])+0.5) + 0.1)
            x_val = [min(x_max_cmd, max(col1_data)+0.2) - 0.4]*len(mag_y)
            plt.errorbar(x_val, mag_y, yerr=func(mag_y, *popt_mag), 
                         xerr=func(mag_y, *popt_col1), fmt='k.', lw=0.8, ms=0.,\
                         zorder=4)
            # Plot colorbar.
            cbar_posit = [0.65, 0.417, 0.04, 0.005] if bf_flag and \
            best_fit_algor == 'genet' else [0.65, 0.408, 0.04, 0.005]
            cbaxes = fig.add_axes(cbar_posit)
            cbar = plt.colorbar(cax=cbaxes, ticks=[v_min,v_max],
                                orientation='horizontal')
            cbar.ax.tick_params(labelsize=9)
            
            # Synthetic cluster.
            ax19 = plt.subplot(gs1[8:10, 6:8])
            #Set plot limits
            plt.xlim(col1_min, col1_max)
            plt.ylim(mag_min, mag_max)
            #Set axis labels
            plt.xlabel('$'+x_ax+'$', fontsize=18)
            plt.ylabel('$'+y_ax+'$', fontsize=18)
            # Set minor ticks
            ax19.minorticks_on()
            ax19.xaxis.set_major_locator(MultipleLocator(1.0))
            ax19.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
            # Add text box
            m, a, e, d = isoch_fit_params[0]
            e_m, e_a, e_e, e_d = isoch_fit_errors
            text1 = '$z = %0.4f \pm %0.4f$' '\n' % (m, e_m)
            text2 = '$log(age) = %0.2f \pm %0.2f$' '\n' % (a, e_a)
            text3 = '$E_{(B-V)} = %0.2f \pm %0.2f$' '\n' % (e, e_e)
            text4 = '$(m-M)_o = %0.2f \pm %0.2f$' % (d, e_d)
            text = text1+text2+text3+text4
            plt.text(0.1, 0.8, text, transform = ax19.transAxes, 
                     bbox=dict(facecolor='white', alpha=0.6), fontsize=12)
            # Plot isochrone.
            plt.plot(shift_isoch[0], shift_isoch[1], 'r', lw=1.2)
            # Plot synth clust.
            plt.scatter(synth_clst[0], synth_clst[2], marker='o', s=30, c='#4682b4',
                        lw=0.5)

    # Best fitting process plots for GA.
    if bf_flag and best_fit_algor == 'genet':

        # Set ranges used by plots below.
        m_min, m_max = ps_params[3][0], ps_params[3][1]
        a_min, a_max = ps_params[4][0], ps_params[4][1]
        e_min, e_max = ps_params[5][0], ps_params[5][1]
        d_min, d_max = ps_params[6][0], ps_params[6][1]
        if m_min == m_max:
            m_min, m_max = m_min-0.1*m_min, m_max+0.1*m_min
        if a_min == a_max:
            a_min, a_max = a_min-0.1*a_min, a_max+0.1*a_min
        if e_min == e_max:
            e_min, e_max = e_min-0.1*e_min, e_max+0.1*e_min
        if d_min == d_max:
            d_min, d_max = d_min-0.1*d_min, d_max+0.1*d_min
        
        # Age vs metallicity GA diagram.
        isoch_done = isoch_fit_params[3]
        plt.subplot(gs1[10:12, 0:2])
        # Axis limits.
        plt.xlim(m_min, m_max)
        plt.ylim(a_min, a_max)
        plt.xlabel('$z$', fontsize=16)
        plt.ylabel('$log(age)$', fontsize=16)
        # Plot best fit point.
        plt.scatter(m, a, marker='o', c='r', s=30)
        # Plot ellipse error.
        ax20 = plt.gca()
        ellipse = Ellipse(xy=(m, a), width=2*e_m, height=2*e_a, 
                                edgecolor='r', fc='None', lw=1.)
        ax20.add_patch(ellipse)
        # Plot density map.
        hist, xedges, yedges = np.histogram2d(zip(*isoch_done[0])[0],
                                              zip(*isoch_done[0])[1], bins=100)
        # H_g is the 2D histogram with a gaussian filter applied
        h_g = gaussian_filter(hist, 2, mode='constant')        
        plt.imshow(h_g.transpose(), origin='lower',
                   extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],\
                   cmap=plt.get_cmap('Blues'), aspect='auto')
        
        # GA diagram.
        lkl_old, ext_imm_indx = isoch_fit_params[1], isoch_fit_params[2]
        ax21 = plt.subplot(gs1[10:12, 2:6])
        plt.xlim(-0.5, n_gen+int(0.01*n_gen))
        plt.ylim(max(0, min(lkl_old[0])-0.3*min(lkl_old[0])),
                 max(lkl_old[1])+min(lkl_old[0])/2.)
        ax21.tick_params(axis='y', which='major', labelsize=9)
        ax21.grid(b=True, which='major', color='gray', linestyle='--', lw=0.6)
        plt.xlabel('Generation', fontsize=12)
        plt.ylabel('Likelihood', fontsize=12)
        text1 = '$N = %d\,;\,L_{min}=%0.2f$' '\n' % (len(lkl_old[0]),\
        min(lkl_old[0]))
        text2 = '$n_{gen}=%d\,;\,n_{pop}=%d$' '\n' % (n_gen, n_pop)
        text3 = '$f_{dif}=%0.2f\,;\,cr_{sel}=%s$' '\n' % (fdif, cr_sel)
        text4 = '$p_{cross}=%0.2f\,;\,p_{mut}=%0.2f$' '\n' % (p_cross, p_mut)
        text5 = '$n_{el}=%d\,;\,n_{ei}=%d\,;\,n_{es}=%d$' % (n_el, n_ei, n_es)
        text = text1+text2+text3+text4+text5
        plt.text(0.05, 0.75, text, transform = ax21.transAxes, \
        bbox=dict(facecolor='white', alpha=0.75), fontsize=12)
        # Plot likelihood minimum and mean lines.
        ax21.plot(range(len(lkl_old[0])), lkl_old[0], lw=1., c='black',
                  label='$L_{min}$')
        ax21.plot(range(len(lkl_old[0])), lkl_old[1], lw=1., c='blue',
                  label='$L_{mean}$')
        # Plot extinction/immigration lines.
        for lin in ext_imm_indx:
            plt.axvline(x=lin, linestyle='--', lw=2., color='green')
        # Legend.
        handles, labels = ax21.get_legend_handles_labels()
        leg = ax21.legend(handles, labels, loc='upper right', numpoints=1,
                          fontsize=12)
        leg.get_frame().set_alpha(0.6)
        
        # Extinction vs distance modulus GA diagram.
        isoch_done = isoch_fit_params[3]
        plt.subplot(gs1[10:12, 6:8])
        plt.xlim(e_min, e_max)
        plt.ylim(d_min, d_max)
        plt.xlabel('$E_{(B-V)}$', fontsize=16)
        plt.ylabel('$(m-M)_o$', fontsize=16)
        # Plot best fit point.
        plt.scatter(e, d, marker='o', c='b', s=30)
        # Plot ellipse error.
        ax21 = plt.gca()
        ellipse = Ellipse(xy=(e, d), width=2*e_e, height=2*e_d, 
                                edgecolor='b', fc='None', lw=1.)
        ax21.add_patch(ellipse)
        # Plot density map.
        hist, xedges, yedges = np.histogram2d(zip(*isoch_done[0])[2],
                                              zip(*isoch_done[0])[3], bins=100)
        # H_g is the 2D histogram with a gaussian filter applied
        h_g = gaussian_filter(hist, 2, mode='constant')        
        plt.imshow(h_g.transpose(), origin='lower',
                   extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],\
                   cmap=plt.get_cmap('Reds'), aspect='auto')


        ax22 = plt.subplot(gs1[12:14, 0:2])
        plt.ylim(max(0, min(lkl_old[0])-0.3*min(lkl_old[0])), max(lkl_old[1]))
        plt.xlim(m_min, m_max)
        ax22.tick_params(axis='y', which='major', labelsize=9)
        plt.ylabel('Likelihood', fontsize=12)
        plt.xlabel('$z$', fontsize=16)
        text = '$z = %0.4f \pm %0.4f$' % (m, e_m)
        plt.text(0.1, 0.93, text, transform = ax22.transAxes, \
        bbox=dict(facecolor='white', alpha=0.5), fontsize=12)
        hist, xedges, yedges = np.histogram2d(zip(*isoch_done[0])[0],
                                              isoch_done[1], bins=100)
        # H_g is the 2D histogram with a gaussian filter applied
        h_g = gaussian_filter(hist, 2, mode='constant')
        y_min_edge = max(0, min(lkl_old[0])-0.3*min(lkl_old[0]))
        plt.imshow(h_g.transpose(), origin='lower',
                   extent=[xedges[0], xedges[-1], y_min_edge, yedges[-1]],\
                   cmap=plt.get_cmap('gist_yarg'), aspect='auto')
        plt.axvline(x=m, linestyle='--', color='blue')
        plt.axvline(x=m+e_m, linestyle='--', color='red')
        plt.axvline(x=m-e_m, linestyle='--', color='red')
    
        
        ax23 = plt.subplot(gs1[12:14, 2:4])
        plt.ylim(max(0, min(lkl_old[0])-0.3*min(lkl_old[0])), max(lkl_old[1]))
        plt.xlim(a_min, a_max)
        ax23.tick_params(axis='y', which='major', labelsize=9)
        plt.ylabel('Likelihood', fontsize=12)
        plt.xlabel('$log(age)$', fontsize=16)
        text = '$log(age) = %0.2f \pm %0.2f$' % (a, e_a)
        plt.text(0.1, 0.93, text, transform = ax23.transAxes, \
        bbox=dict(facecolor='white', alpha=0.5), fontsize=12)
        hist, xedges, yedges = np.histogram2d(zip(*isoch_done[0])[1],
                                              isoch_done[1], bins=100)
        # H_g is the 2D histogram with a gaussian filter applied
        h_g = gaussian_filter(hist, 2, mode='constant')        
        plt.imshow(h_g.transpose(), origin='lower',
                   extent=[xedges[0], xedges[-1], y_min_edge, yedges[-1]],\
                   cmap=plt.get_cmap('gist_yarg'), aspect='auto')
        plt.axvline(x=a, linestyle='--', color='blue')
        plt.axvline(x=a+e_a, linestyle='--', color='red')
        plt.axvline(x=a-e_a, linestyle='--', color='red')
    
    
        ax24 = plt.subplot(gs1[12:14, 4:6])
        plt.ylim(max(0, min(lkl_old[0])-0.3*min(lkl_old[0])), max(lkl_old[1]))
        plt.xlim(e_min, e_max)
        ax24.tick_params(axis='y', which='major', labelsize=9)
        plt.ylabel('Likelihood', fontsize=12)
        plt.xlabel('$E_{(B-V)}$', fontsize=16)
        text = '$E_{(B-V)} = %0.2f \pm %0.2f$' % (e, e_e)
        plt.text(0.1, 0.93, text, transform = ax24.transAxes, \
        bbox=dict(facecolor='white', alpha=0.5), fontsize=12)
        hist, xedges, yedges = np.histogram2d(zip(*isoch_done[0])[2],
                                              isoch_done[1], bins=100)
        # H_g is the 2D histogram with a gaussian filter applied
        h_g = gaussian_filter(hist, 2, mode='constant')        
        plt.imshow(h_g.transpose(), origin='lower',
                   extent=[xedges[0], xedges[-1], y_min_edge, yedges[-1]],\
                   cmap=plt.get_cmap('gist_yarg'), aspect='auto')
        plt.axvline(x=e, linestyle='--', color='blue')
        plt.axvline(x=e+e_e, linestyle='--', color='red')
        plt.axvline(x=e-e_e, linestyle='--', color='red')
        
        
        ax25 = plt.subplot(gs1[12:14, 6:8])
        plt.ylim(max(0, min(lkl_old[0])-0.3*min(lkl_old[0])), max(lkl_old[1]))
        plt.xlim(d_min, d_max)
        ax25.tick_params(axis='y', which='major', labelsize=9)
        plt.ylabel('Likelihood', fontsize=12)
        plt.xlabel('$(m-M)_o$', fontsize=16)
        text = '$(m-M)_o = %0.2f \pm %0.2f$' % (d, e_d)
        plt.text(0.1, 0.93, text, transform = ax25.transAxes, \
        bbox=dict(facecolor='white', alpha=0.5), fontsize=12)
        hist, xedges, yedges = np.histogram2d(zip(*isoch_done[0])[3],
                                              isoch_done[1], bins=100)
        # H_g is the 2D histogram with a gaussian filter applied
        h_g = gaussian_filter(hist, 2, mode='constant')        
        plt.imshow(h_g.transpose(), origin='lower',
                   extent=[xedges[0], xedges[-1], y_min_edge, yedges[-1]],\
                   cmap=plt.get_cmap('gist_yarg'), aspect='auto')
        plt.axvline(x=d, linestyle='--', color='blue')
        plt.axvline(x=d+e_d, linestyle='--', color='red')
        plt.axvline(x=d-e_d, linestyle='--', color='red')


    fig.tight_layout()

    # Generate output file for each data file.
    plt.savefig(join(output_subdir, str(clust_name)+'.png'), dpi=150)
    
    # Close to release memory.
    plt.clf()
    plt.close()
    

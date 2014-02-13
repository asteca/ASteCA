"""
@author: gabriel
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def func(x, a, b, c):
    '''
    Exponential function.
    '''
    return a * np.exp(b * x) + c   


def disp_errors(mag_data, popt_mag,  popt_col1, acpt_stars, rjct_stars,
                err_plot, er_params):
    '''
    Plot errors diagrams.
    '''

    bright_end, popt_umag, pol_mag, popt_ucol1, pol_col1, mag_val_left, \
    mag_val_right, col1_val_left, col1_val_right = err_plot

    be, be_e, e_max = er_params
    
    # Plot all outputs
    plt.figure(figsize=(10, 8)) # create the top-level container
    gs = gridspec.GridSpec(2, 2)  # create a GridSpec object
        
    # T1 magnitude error

    # Store magnitude, colors and their errors.
    stars_rjct_temp = [[], []]
    for star in rjct_stars:
        stars_rjct_temp[0].append(star[3])
        stars_rjct_temp[1].append(star[4])
    stars_acpt_temp = [[], []]
    for star in acpt_stars:
        stars_acpt_temp[0].append(star[3])
        stars_acpt_temp[1].append(star[4])      
    
    ax7 = plt.subplot(gs[0, 0:2])
    #Set plot limits
    plt.xlim(min(mag_data)-0.5, max(mag_data)+0.5)
    plt.ylim(-0.005, min(0.5, stars_rjct_temp[1]))
    #Set axis labels
    plt.ylabel(r'$\sigma_{T_1}$', fontsize=18)
    plt.xlabel(r'$T_1$', fontsize=18)
    # Set minor ticks
    ax7.minorticks_on()
    # Plot exponential fitted function.
    mag_x = np.linspace(min(mag_data), max(mag_data), 50)
    # Plot lower envelope.
    ax7.plot(mag_x, func(mag_x, *popt_mag), 'r-', zorder=3)
    # Plot left side of upper envelope (exponential).
    ax7.plot(mag_val_left, func(mag_val_left, *popt_umag), 'r--', lw=2.,
             zorder=3)
    # Plot right side of upper envelope (polynomial).
    ax7.plot(mag_val_right, np.polyval(pol_mag, (mag_val_right)), 'r--', lw=2.,
             zorder=3)
    # Plot rectangle.
    ax7.vlines(x=bright_end+0.05, ymin=-0.005, ymax=be_e, color='r', 
               linestyles='dashed', zorder=2)
    ax7.vlines(x=min(mag_data)-0.05, ymin=-0.005, ymax=be_e, color='r',
               linestyles='dashed', zorder=2)
    ax7.hlines(y=be_e, xmin=min(mag_data), xmax=bright_end, color='r',
               linestyles='dashed', zorder=2)
    # Plot stars.
    plt.scatter(stars_rjct_temp[0], stars_rjct_temp[1], marker='x', c='teal', 
                s=15, zorder=1)
    plt.scatter(stars_acpt_temp[0], stars_acpt_temp[1], marker='o', c='k', 
                s=1, zorder=2)


    # C-T1 color error

    stars_rjct_temp = [[], []]
    for star in rjct_stars:
        stars_rjct_temp[0].append(star[3])
        stars_rjct_temp[1].append(star[6])
    stars_acpt_temp = [[], []]
    for star in acpt_stars:
        stars_acpt_temp[0].append(star[3])
        stars_acpt_temp[1].append(star[6])    
    
    ax8 = plt.subplot(gs[1, 0:2])
    #Set plot limits
    plt.xlim(min(mag_data)-0.5, max(mag_data)+0.5)
    plt.ylim(-0.005, min(0.5, stars_rjct_temp[1]))
    #Set axis labels
    plt.ylabel(r'$\sigma_{(C-T_1)}$', fontsize=18)
    plt.xlabel(r'$T_1$', fontsize=18)
    # Set minor ticks
    ax8.minorticks_on()
    # Plot lower envelope.
    ax8.plot(mag_x, func(mag_x, *popt_col1), 'r-', zorder=3)
    # Plot left side of upper envelope (exponential).
    ax8.plot(col1_val_left, func(col1_val_left, *popt_ucol1), 'r--', lw=2.,
             zorder=3)
    # Plot right side of upper envelope (polynomial).
    ax8.plot(col1_val_right, np.polyval(pol_col1, (col1_val_right)), 'r--',
             lw=2., zorder=3)
    # Plot rectangle.
    ax8.vlines(x=bright_end+0.05, ymin=-0.005, ymax=be_e, color='r', 
               linestyles='dashed', zorder=2)
    ax8.vlines(x=min(mag_data)-0.05, ymin=-0.005, ymax=be_e, color='r',
               linestyles='dashed', zorder=2)
    ax8.hlines(y=be_e, xmin=min(mag_data), xmax=bright_end, color='r',
               linestyles='dashed', zorder=2)
    # Plot stars.
    plt.scatter(stars_rjct_temp[0], stars_rjct_temp[1], marker='x', c='teal', 
                s=15, zorder=1)
    plt.scatter(stars_acpt_temp[0], stars_acpt_temp[1], marker='o', c='k', 
                s=1, zorder=2)

    plt.draw()
    print 'Plot displayed, waiting for it to be closed.'


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from ..inp import input_params as g
from ..math_f import exp_function


def main(er_mode, mag, err_plot, acpt_stars, rjct_stars, err_pck):
    '''
    Plot errors diagrams.
    '''

    # Error parameters.
    be_m, interv_mag, mag_value = err_pck
    e_max, be, be_e = g.er_params[1:4]
    # Define names for CMD axes.
    y_ax, x_ax0, m_ord = g.axes_params
    if m_ord == 21:
        x_ax = '(' + x_ax0 + '-' + y_ax + ')'
    elif m_ord == 12:
        x_ax = '(' + y_ax + '-' + x_ax0 + ')'

    # Plot all outputs
    plt.figure(figsize=(10, 8))  # create the top-level container
    gs = gridspec.GridSpec(2, 2)  # create a GridSpec object

    # Magnitude error
    stars_rjct_temp = [[], []]
    for star in rjct_stars:
        stars_rjct_temp[0].append(star[3])
        stars_rjct_temp[1].append(star[4])
    stars_acpt_temp = [[], []]
    for star in acpt_stars:
        stars_acpt_temp[0].append(star[3])
        stars_acpt_temp[1].append(star[4])

    axm = plt.subplot(gs[0, 0:2])
    # Set plot limits
    x_min, x_max = min(mag) - 0.5, max(mag) + 0.5
    plt.xlim(x_min, x_max)
    plt.ylim(-0.005, e_max + (e_max / 5.))
    # Set axis labels
    plt.ylabel('$\sigma_{' + y_ax + '}$', fontsize=18)
    plt.xlabel('$' + y_ax + '$', fontsize=18)
    # Set minor ticks
    axm.minorticks_on()
    # Plot e_max line.
    axm.hlines(y=e_max, xmin=x_min, xmax=x_max, color='r',
               linestyles='dashed', zorder=2)
    # Plot rectangle.
    bright_end = min(mag) + be
    axm.vlines(x=bright_end + 0.05, ymin=-0.005, ymax=be_e, color='r',
               linestyles='dashed', zorder=2)
    axm.vlines(x=min(mag) - 0.05, ymin=-0.005, ymax=be_e, color='r',
               linestyles='dashed', zorder=2)
    axm.hlines(y=be_e, xmin=min(mag), xmax=bright_end, color='r',
               linestyles='dashed', zorder=2)
    # Plot curve(s) according to the method used.
    if er_mode == 'eyefit':
        # Unpack params.
        val_mag, pol_mag, val_col, pol_col, mag_val_left, \
            mag_val_right, col_val_left, col_val_right = err_plot
        # Combine left + right values.
        m_v, e_v = mag_val_left + mag_val_right, \
            [val_mag for _ in mag_val_left] + list(np.polyval(pol_mag,
                                                   (mag_val_right)))
        axm.plot(m_v, e_v, 'k-', zorder=3)

    elif er_mode == 'lowexp':
        # Unpack params.
        popt_mag, popt_col1 = err_plot
        # Plot exponential curve.
        mag_x = np.linspace(bright_end, max(mag), 50)
        axm.plot(mag_x, exp_function.exp_2p(mag_x, *popt_mag), 'r-', zorder=3)
    # Plot stars.
    plt.scatter(stars_rjct_temp[0], stars_rjct_temp[1], marker='x', c='teal',
                s=15, zorder=1)
    plt.scatter(stars_acpt_temp[0], stars_acpt_temp[1], marker='o', c='k',
                s=1, zorder=2)

    # Color error
    stars_rjct_temp = [[], []]
    for star in rjct_stars:
        stars_rjct_temp[0].append(star[3])
        stars_rjct_temp[1].append(star[6])
    stars_acpt_temp = [[], []]
    for star in acpt_stars:
        stars_acpt_temp[0].append(star[3])
        stars_acpt_temp[1].append(star[6])

    axc1 = plt.subplot(gs[1, 0:2])
    # Set plot limits
    x_min, x_max = min(mag) - 0.5, max(mag) + 0.5
    plt.xlim(x_min, x_max)
    plt.ylim(-0.005, e_max + (e_max / 5.))
    # Set axis labels
    plt.ylabel('$\sigma_{' + x_ax + '}$', fontsize=18)
    plt.xlabel('$' + y_ax + '$', fontsize=18)
    # Set minor ticks
    axc1.minorticks_on()
    # Plot e_max line.
    axc1.hlines(y=e_max, xmin=x_min, xmax=x_max, color='r',
                linestyles='dashed', zorder=2)
    # Plot rectangle.
    bright_end = min(mag) + be
    axc1.vlines(x=bright_end + 0.05, ymin=-0.005, ymax=be_e, color='r',
                linestyles='dashed', zorder=2)
    axc1.vlines(x=min(mag) - 0.05, ymin=-0.005, ymax=be_e, color='r',
                linestyles='dashed', zorder=2)
    axc1.hlines(y=be_e, xmin=min(mag), xmax=bright_end, color='r',
                linestyles='dashed', zorder=2)
    # Plot curve(s) according to the method used.
    if er_mode == 'eyefit':
        # Unpack params.
        val_mag, pol_mag, val_col, pol_col, mag_val_left, \
            mag_val_right, col_val_left, col_val_right = err_plot
        # Combine left + right values.
        m_v, e_v = col_val_left + col_val_right, \
            [val_col for _ in col_val_left] + list(np.polyval(pol_col,
                                                   (col_val_right)))
        axc1.plot(m_v, e_v, 'k-', zorder=3)

    elif er_mode == 'lowexp':
        # Unpack params.
        popt_mag, popt_col1 = err_plot
        # Plot exponential curve.
        axc1.plot(mag_x, exp_function.exp_2p(mag_x, *popt_col1), 'r-',
                  zorder=3)
    # Plot stars.
    plt.scatter(stars_rjct_temp[0], stars_rjct_temp[1], marker='x', c='teal',
                s=15, zorder=1)
    plt.scatter(stars_acpt_temp[0], stars_acpt_temp[1], marker='o', c='k',
                s=1, zorder=2)

    plt.draw()
    print 'Plot displayed, waiting for it to be closed.'

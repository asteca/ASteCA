
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.offsetbox as offsetbox
from ..out import prep_plots


def main(x_data, y_data, mag_data, center_cl, cent_bin, h_filter):
    '''
    Show plot of cluster with value of center obtained.
    '''

    coord, x_name, y_name = prep_plots.coord_syst()
    st_sizes_arr = prep_plots.star_size(mag_data)

    # Plot all outputs
    plt.figure(figsize=(18, 8))  # create the top-level container
    gs = gridspec.GridSpec(1, 2)  # create a GridSpec object

    # 1 subplot: 2D filtered histogram
    ax1 = plt.subplot(gs[0, 0])
    # Set axis labels.
    plt.xlabel('{} (bins)'.format(x_name), fontsize=12)
    plt.ylabel('{} (bins)'.format(y_name), fontsize=12)
    # Set minor ticks
    ax1.minorticks_on()
    # Set grid
    ax1.grid(b=True, which='major', color='k', linestyle='--', zorder=3)
    ax1.grid(b=True, which='minor', color='k', linestyle='--', zorder=3)
    plt.axvline(x=cent_bin[0], linestyle='--', color='white')
    plt.axhline(y=cent_bin[1], linestyle='--', color='white')
    plt.imshow(h_filter.transpose(), origin='lower')
    if coord == 'deg':
        # If RA is used, invert axis.
        plt.gca().invert_xaxis()

    # 2 subplot: x,y finding chart of full frame
    ax2 = plt.subplot(gs[0, 1])
    ax2.set_aspect('equal')
    # Get max and min values in x,y
    x_min, x_max = min(x_data), max(x_data)
    y_min, y_max = min(y_data), max(y_data)
    # Set plot limits
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    if coord == 'deg':
        # If RA is used, invert axis.
        plt.gca().invert_xaxis()
    # Set axis labels.
    plt.xlabel('{} ({})'.format(x_name, coord), fontsize=12)
    plt.ylabel('{} ({})'.format(y_name, coord), fontsize=12)
    # Set minor ticks
    ax2.minorticks_on()
    # Set grid
    ax2.grid(b=True, which='major', color='k', linestyle='-', zorder=1)
    ax2.grid(b=True, which='minor', color='k', linestyle='-', zorder=1)
    # Add lines through the center of the cluster
    plt.axvline(x=center_cl[0], linestyle='--', color='red', lw=2.)
    plt.axhline(y=center_cl[1], linestyle='--', color='red', lw=2.)
    # Add text box.
    text1 = '${0}_{{cent}} = {1:g}\,{2}$'.format(x_name, center_cl[0], coord)
    text2 = '${0}_{{cent}} = {1:g}\,{2}$'.format(y_name, center_cl[1], coord)
    text = text1 + '\n' + text2
    ob = offsetbox.AnchoredText(text, loc=2, prop=dict(size=11))
    ob.patch.set(boxstyle='square,pad=-0.2', alpha=0.85)
    plt.gca().add_artist(ob)
    # Plot stars.
    plt.scatter(x_data, y_data, marker='o', c='black', s=st_sizes_arr)

    plt.draw()
    print('Plot displayed, waiting for it to be closed.')

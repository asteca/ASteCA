
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.offsetbox as offsetbox
from ..out import prep_plots


def main(phot_data, bin_width, center_params, clust_rad, e_rad, field_dens,
         rdp_params):
    '''
    Plot cluster and its radius.
    '''

    # Unpack.
    x_data, y_data, mag_data = phot_data[1], phot_data[2], phot_data[3]
    cent_bin, kde_cent, e_cent = center_params[:3]
    hist_2d_g = center_params[5]
    center_cl = [kde_cent[0], kde_cent[1]]
    x_center_bin, y_center_bin = cent_bin
    radii, ring_density, poisson_error = rdp_params[:3]

    coord, x_name, y_name = prep_plots.coord_syst()
    st_sizes_arr = prep_plots.star_size(mag_data)

    # Plot all outputs
    fig = plt.figure(figsize=(12, 12))
    gs = gridspec.GridSpec(2, 2)  # create a GridSpec object

    # 2D filtered histogram, smallest bin width.
    ax1 = plt.subplot(gs[0, 0])
    # Set axis labels.
    plt.xlabel('{} (bins)'.format(x_name), fontsize=12)
    plt.ylabel('{} (bins)'.format(y_name), fontsize=12)
    ax1.minorticks_on()
    plt.axvline(x=x_center_bin, linestyle='--', color='white')
    plt.axhline(y=y_center_bin, linestyle='--', color='white')
    # Radius
    circle = plt.Circle((x_center_bin, y_center_bin),
                        clust_rad / bin_width, color='w', fill=False)
    fig.gca().add_artist(circle)
    # Add text boxs.
    text = 'Bin $\simeq$ {:.0f} {}'.format(bin_width, coord)
    ob = offsetbox.AnchoredText(text, loc=1, prop=dict(size=10))
    ob.patch.set(boxstyle='square,pad=-0.2', alpha=0.85)
    fig.gca().add_artist(ob)
    plt.imshow(hist_2d_g.transpose(), origin='lower')
    if coord == 'deg':
        # If RA is used, invert axis.
        plt.gca().invert_xaxis()

    # Finding chart.
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
    circle = plt.Circle((center_cl[0], center_cl[1]), clust_rad, color='r',
                        fill=False)
    fig.gca().add_artist(circle)
    # Add text box
    text1 = '${0}_{{cent}} = {1:g}\,{2}$'.format(x_name, center_cl[0], coord)
    text2 = '${0}_{{cent}} = {1:g}\,{2}$'.format(y_name, center_cl[1], coord)
    text = text1 + '\n' + text2
    ob = offsetbox.AnchoredText(text, loc=2, prop=dict(size=11))
    ob.patch.set(boxstyle='square,pad=-0.2', alpha=0.85)
    plt.gca().add_artist(ob)
    # Plot stars.
    plt.scatter(x_data, y_data, marker='o', c='black', s=st_sizes_arr)

    # Radial density plot.
    ax3 = plt.subplot(gs[1, 0:2])
    # Get max and min values in x,y
    x_min, x_max = min(radii) - (max(radii) / 10.), \
        max(radii) + (max(radii) / 10.)
    max(radii) + (max(radii) / 20.)
    delta_total = (max(ring_density) - field_dens)
    delta_backg = 0.2 * delta_total
    y_min = (field_dens - delta_backg) - (max(ring_density) -
                                          min(ring_density)) / 10
    y_max = max(ring_density) + (max(ring_density) - min(ring_density)) / 10
    # Set plot limits
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    # Set axes labels.
    plt.xlabel('radius ({})'.format(coord), fontsize=12)
    plt.ylabel("stars/{}$^{{2}}$".format(coord), fontsize=12)
    # Legend texts
    texts = ['RDP ($\sim${:.0f} {})'.format(bin_width, coord),
             '$d_{{field}}$ = {:.1E} $st/{}^{{2}}$'.format(field_dens, coord),
             'r$_{{cl}}$ = {0:g} $\pm$ {1:g} {2}'.format(clust_rad, e_rad,
                                                         coord)]
    # Plot density profile with the smallest bin size
    ax3.plot(radii, ring_density, 'ko-', zorder=3, label=texts[0])
    # Plot poisson error bars
    plt.errorbar(radii, ring_density, yerr=poisson_error, fmt='ko',
                 zorder=1)
    # Plot background level.
    ax3.hlines(y=field_dens, xmin=0, xmax=max(radii),
               label=texts[1], color='b', zorder=5)
    # Approx middle of the graph.
    arr_y_up = (y_max - y_min) / 2.3 + y_min
    # Length of arrow head.
    head_w, head_l = x_max * 0.023, (y_max - y_min) * 0.045
    # Length of arrow.
    arr_y_dwn = -1. * abs(arr_y_up - field_dens) * 0.76
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
    print('Plot displayed, waiting for it to be closed.')

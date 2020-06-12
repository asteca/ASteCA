
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.offsetbox as offsetbox


def pl_full_frame(
    N, gs, fig, x_name, y_name, coord, x_min, x_max, y_min, y_max, asp_ratio,
        x, y, st_sizes_arr, mag_range, main_mag):
    """
    x,y finding chart of stars in frame within 'mag_range'
    """
    gs_map = {
        0: gs[0:2, 0:2], 1: gs[2:4, 0:2], 2: gs[4:6, 0:2],
        3: gs[6:8, 0:2], 4: gs[8:10, 0:2]}
    ax = plt.subplot(gs_map.get(N))

    ax.set_aspect(aspect=asp_ratio)
    plt.title("{}: {}".format(main_mag, mag_range))
    # Set plot limits
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    # If RA is used, invert axis.
    if coord == 'deg':
        ax.invert_xaxis()
    # Set axis labels
    plt.xlabel('{} ({})'.format(x_name, coord))
    plt.ylabel('{} ({})'.format(y_name, coord))
    # Set minor ticks
    ax.minorticks_on()
    # Plot stars.
    plt.scatter(x, y, marker='o', c='black', s=st_sizes_arr)


def pl_densxy(
    N, gs, fig, asp_ratio, x_name, y_name, coord, bw_list, kde_pl,
        cent_xy):
    """
    2D Gaussian convolved histogram.
    """
    gs_map = {
        5: gs[0:2, 2:4], 6: gs[0:2, 4:6], 7: gs[0:2, 6:8],
        8: gs[2:4, 2:4], 9: gs[2:4, 4:6], 10: gs[2:4, 6:8],
        11: gs[4:6, 2:4], 12: gs[4:6, 4:6], 13: gs[4:6, 6:8],
        14: gs[6:8, 2:4], 15: gs[6:8, 4:6], 16: gs[6:8, 6:8],
        17: gs[8:10, 2:4], 18: gs[8:10, 4:6], 19: gs[8:10, 6:8]
    }
    ax = plt.subplot(gs_map.get(N))

    if N in [5, 6, 7]:
        r_frmt = '{:.0f}' if coord == 'px' else '{:.4f}'
        plt.title(("Bandwidth: " + r_frmt + ' [{}]').format(
            bw_list[N - 5], coord), fontsize=11)
    plt.xlabel('{} ({})'.format(x_name, coord), fontsize=11)
    plt.ylabel('{} ({})'.format(y_name, coord), fontsize=11)

    plt.axvline(x=cent_xy[0], linestyle='--', lw=.85, color='green')
    plt.axhline(y=cent_xy[1], linestyle='--', lw=.85, color='green')
    plt.scatter(*cent_xy, marker='x', color='w', s=20)

    if kde_pl:
        ext_range, x, y, k_pos = kde_pl
        kde = np.reshape(k_pos.T, x.shape)
        plt.imshow(
            np.rot90(kde), cmap=plt.get_cmap('RdYlBu_r'), extent=ext_range)
        plt.contour(x, y, kde, colors='#551a8b', linewidths=0.5)

        # Add text box
        r_frmt = '{:.0f}' if coord == 'px' else '{:.3f}'
        t1 = ('${}_{{c}} =$' + r_frmt + '$\,{}$').format(
            x_name, cent_xy[0], coord)
        t2 = ('${}_{{c}} =$' + r_frmt + '$\,{}$').format(
            y_name, cent_xy[1], coord)
        text = t1 + '\n' + t2
        ob = offsetbox.AnchoredText(text, pad=0.2, loc=2, prop=dict(size=10))
        ob.patch.set(alpha=0.85)
        ax.add_artist(ob)

    # If RA is used, invert axis.
    if coord == 'deg':
        ax.invert_xaxis()
    ax.set_aspect(aspect=asp_ratio)


def plot(N, *args):
    """
    Handle each plot separately.
    """
    if N < 5:
        plt_map = dict.fromkeys(
            [0, 1, 2, 3, 4], [pl_full_frame, 'full frame'])
    else:
        plt_map = dict.fromkeys(
            [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19],
            [pl_densxy, 'xy density map'])

    fxn = plt_map.get(N, None)[0]
    if fxn is None:
        raise ValueError("  ERROR: there is no plot {}.".format(N))

    try:
        fxn(N, *args)
    except Exception:
        print("  WARNING: error when plotting {}".format(plt_map.get(N)[1]))
        import traceback
        print(traceback.format_exc())


import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


def main(filters, colors, mmag, acpt_stars, rjct_stars, err_max):
    '''
    Plot errors diagrams.
    '''
    # Define names for CMD axes.
    x_ax = colors[0][1].replace(',', '-')
    y_ax = filters[0][1]

    # Plot all outputs
    plt.figure(figsize=(10, 8))  # create the top-level container
    gs = gridspec.GridSpec(2, 2)  # create a GridSpec object

    # Magnitude error
    axm = plt.subplot(gs[0, 0:2])
    # Set plot limits
    x_min, x_max = min(mmag) - 0.5, max(mmag) + 0.5
    plt.xlim(x_min, x_max)
    plt.ylim(-0.005, err_max + (err_max / 5.))
    # Set axis labels
    plt.ylabel('$\sigma_{' + y_ax + '}$', fontsize=18)
    plt.xlabel('$' + y_ax + '$', fontsize=18)
    # Set minor ticks
    axm.minorticks_on()
    # Plot err_max line.
    axm.hlines(y=err_max, xmin=x_min, xmax=x_max, color='r',
               linestyles='dashed', zorder=2)
    # Plot stars.
    # Index 3,4 point to magnitude and its error. Index 0 takes the first
    # one defined and stored.
    plt.scatter(zip(*zip(*rjct_stars)[3])[0], zip(*zip(*rjct_stars)[4])[0],
                marker='x', c='teal', s=15, zorder=1)
    plt.scatter(zip(*zip(*acpt_stars)[3])[0], zip(*zip(*acpt_stars)[4])[0],
                marker='o', c='k', s=1, zorder=2)

    # Color error
    axc1 = plt.subplot(gs[1, 0:2])
    # Set plot limits
    x_min, x_max = min(mmag) - 0.5, max(mmag) + 0.5
    plt.xlim(x_min, x_max)
    plt.ylim(-0.005, err_max + (err_max / 5.))
    # Set axis labels
    plt.ylabel('$\sigma_{' + x_ax + '}$', fontsize=18)
    plt.xlabel('$' + y_ax + '$', fontsize=18)
    # Set minor ticks
    axc1.minorticks_on()
    # Plot err_max line.
    axc1.hlines(y=err_max, xmin=x_min, xmax=x_max, color='r',
                linestyles='dashed', zorder=2)
    # Plot stars.
    plt.scatter(zip(*zip(*rjct_stars)[3])[0], zip(*zip(*rjct_stars)[6])[0],
                marker='x', c='teal', s=15, zorder=1)
    plt.scatter(zip(*zip(*acpt_stars)[3])[0], zip(*zip(*acpt_stars)[6])[0],
                marker='o', c='k', s=1, zorder=2)

    plt.draw()
    print("<<Plot displayed. Will continue after it is closed.>>")

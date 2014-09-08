
import numpy as np
from scipy.ndimage.filters import gaussian_filter


def center_fun(x_data, y_data, d_b):
    '''
    Function that returns 2D histograms for the observed field.
    '''

    xmin, xmax = min(x_data), max(x_data)
    ymin, ymax = min(y_data), max(y_data)
    rang = [[xmin, xmax], [ymin, ymax]]

    # Number of bins in x,y given the bin width 'd_b'
    binsxy = [int((xmax - xmin) / d_b), int((ymax - ymin) / d_b)]
    print 'binsxy', binsxy

    # hist is the 2D histogran, *edges store the edges of the bins.
    hist, xedges, yedges = np.histogram2d(x_data, y_data,
        range=rang, bins=binsxy)

    # H_g is the 2D histogram with a gaussian filter applied.
    h_g = gaussian_filter(hist, 2, mode='constant')

    return hist, xedges, yedges, h_g


def iqr(x):
    '''
    '''
    return np.subtract(*np.percentile(x, [75, 25]))


def fd_rule(x_data, y_data):
    '''
    '''
    h_xy = []
    for data_lst in [x_data, y_data]:
        h_xy.append(2 * iqr(data_lst) * (len(data_lst) ** (-1. / 3.)))

    xmin, xmax = min(x_data), max(x_data)
    ymin, ymax = min(y_data), max(y_data)
    Mxy = [int((xmax - xmin) / h_xy[0]), int((ymax - ymin) / h_xy[1])]
    max_val = Mxy.index(max(Mxy))
    #print h_xy
    #print Mxy, max_val
    #print h_xy[max_val]

    return h_xy[max_val]


def get_2d_hists(x_data, y_data):
    '''
    '''

    xmin, xmax = min(x_data), max(x_data)
    ymin, ymax = min(y_data), max(y_data)
    rang = [[xmin, xmax], [ymin, ymax]]

    # Calculate the number of bins used.
    min_rang = min((rang[0][1] - rang[0][0]), (rang[1][1] - rang[1][0]))
    # Number of bins given by 1%, 2%, 3% and 4% of the minimum axis range.
    bin_list = [(i * min_rang / 100.) for i in range(1, 5)]
    #bw1 = np.sqrt(len(x_data))  # Square root
    ##bw1 = 3.5 * len(x_data) ** (-1./4.)  # Scotts' rule
    #bw2 = fd_rule(x_data, y_data)  # Freedman-Diaconis' rule
    #bw3 = min_rang / 100.  # 1% range
    #bw4 = np.log2(len(x_data)) + 1  # Sturges' rule
    #bin_list = [bw1, bw2, bw3, bw4]

    print bin_list

    hists_2d = []
    # Iterate for the defined bin widths.
    for indx, d_b in enumerate(bin_list):
        print 'indx', indx

        hist, xedges, yedges, h_g = center_fun(x_data, y_data, d_b)
        hists_2d.append(hist)

    import matplotlib.pyplot as plt
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)  # , sharex='col', sharey='row')
    # 2D filtered histograms.
    st_dev = 1.5
    h_g = gaussian_filter(hists_2d[0], st_dev, mode='constant')
    ax1.imshow(h_g.transpose(), origin='lower')
    h_g = gaussian_filter(hists_2d[1], st_dev, mode='constant')
    ax2.imshow(h_g.transpose(), origin='lower')
    h_g = gaussian_filter(hists_2d[2], st_dev, mode='constant')
    ax3.imshow(h_g.transpose(), origin='lower')
    h_g = gaussian_filter(hists_2d[3], st_dev, mode='constant')
    ax4.imshow(h_g.transpose(), origin='lower')
    #ax1.imshow(hists_2d[0].transpose(), origin='lower')
    #ax2.imshow(hists_2d[1].transpose(), origin='lower')
    #ax3.imshow(hists_2d[2].transpose(), origin='lower')
    #ax4.imshow(hists_2d[3].transpose(), origin='lower')

    plt.show()

    raw_input()
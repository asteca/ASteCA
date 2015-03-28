# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 2015

@author: gabriel
"""

import numpy as np
from astroML.plotting import hist
import operator
from decont_algor_bys import sort_members


def rm_stars(decont_algor_return, field_region, local_bin):
    '''
    Takes the photometric diagram (CMD, CCD, etc.) of the cluster region with
    assigned MPs, divides it into su--regions according to the density within
    it, and removes in each sub-regiona number of stars equal to the average
    excess due to field star contamination.
    '''

    memb_prob_avrg_sort = decont_algor_return[0]

    # Remove ID's and zip.
    P = np.array(zip(*memb_prob_avrg_sort)[1:], dtype='float')
    # Create list with all magnitudes and colors defined.
    mag_col_cl = [P[4], P[2]]

    # Obtain bin edges for each photometric dimension using the cluster region
    # diagram.
    bin_edges = []
    if local_bin in ['sturges', 'sqrt']:
        if local_bin == 'sturges':
            b_num = 1 + np.log2(len(P[0]))
        else:
            b_num = np.sqrt(len(P[0]))

        for mag_col in mag_col_cl:
            bin_edges.append(np.histogram(mag_col, bins=b_num)[1])
    else:
        for mag_col in mag_col_cl:
            bin_edges.append(hist(mag_col, bins=local_bin)[1])

    # Cluster region N-dimensional histogram.
    cl_hist = np.histogramdd(np.array(zip(*mag_col_cl)), bins=bin_edges)[0]
    # np.shape(cl_hist) gives the tuple containing one element per dimension,
    # indicating how many cells that dimension was divided into.

    # Add a very small amount to each outer-most edge so the 'np.digitize'
    # function will position the stars on the edges correctly.
    for i, b_e in enumerate(bin_edges):
        bin_edges[i][0] = b_e[0] - (abs(b_e[0]) / 100.)
        bin_edges[i][-1] = b_e[-1] + (b_e[-1] / 100.)

    # Position each cluster region star in its corresponding N-dimensional
    # cell/bin.
    cl_st_indx = []
    # Store indexes for each dimension.
    for i, mag_col in enumerate(mag_col_cl):
        # Set correct indexes for array substracting 1, since 'np.digitize'
        # counts one more bin to the right by default.
        cl_st_indx.append(np.digitize(mag_col, bin_edges[i]) - 1)

    # Create empty list with the same dimension as the photometric N-histogram.
    cl_hist_p = np.empty(shape=cl_hist.shape + (0,)).tolist()

    # Position stars in their corresponding N-histogram cells. Since the stars
    # are already sorted by their MPs, they will be correctly sorted in the
    # final list here too.
    for i, h_indx in enumerate(zip(*cl_st_indx)):
        # Store stars.
        reduce(operator.getitem, list(h_indx), cl_hist_p).append(
            memb_prob_avrg_sort[i])

    # Flatten array to access all of its elements.
    cl_hist_p_flat = np.asarray(cl_hist_p).flatten()

    # Empty field region array shaped like the cluster region array.
    f_hists = np.zeros(shape=np.shape(cl_hist))
    # Add stars in all the defined field regions.
    for freg in field_region:
        Q = np.array(zip(*freg)[1:])
        f_hists = f_hists + np.histogramdd(np.array(zip(*[Q[4], Q[2]])),
            bins=bin_edges)[0]

    # Average number of stars in each cell/bin and round to integer.
    f_hists = np.around(f_hists / len(field_region), 0)

    # Iterate through each N-dimensional cell of both cluster and field region
    # arrays and remove the excess of field stars in each one, selecting
    # those stars with the lowest assigned MPs.
    red_memb_fit, red_memb_no_fit = [], []
    # For each cell defined.
    for i in range(len(cl_hist.flatten())):

        # Get average number of field regions in this cell.
        N_fl_reg = f_hists.flatten()[i]

        if N_fl_reg > 0.:
            # Discard the excess of N_reg_fl stars from this cluster region.
            # Discard those with the smallest MPs, and keep the rest.
            red_memb_fit.append(cl_hist_p_flat[i][:-int(N_fl_reg)])
            red_memb_no_fit.append(cl_hist_p_flat[i][-int(N_fl_reg):])
        else:
            # No field region stars in this cell, keep all stars.
            red_memb_fit.append(cl_hist_p_flat[i])

    # Flatten lists of stars and re-sort according to highest MPs.
    red_memb_fit = sort_members([i for sublst in red_memb_fit for i in sublst])
    red_memb_no_fit = sort_members([i for sublst in red_memb_no_fit for i in
        sublst])

    #import matplotlib.pyplot as plt
    #fig, ax = plt.subplots()
    #plt.scatter(P[4], P[2])
    #ax.set_xticks(bin_edges[0], minor=False)
    #ax.set_yticks(bin_edges[1], minor=False)
    #ax.xaxis.grid(True, which='major')
    #ax.yaxis.grid(True, which='major')
    #ax.invert_yaxis()
    #plt.show()

    # Minimum probability in the list that contains the stars selected.
    red_plot_pars = [red_memb_fit[-1][-1]]

    #return red_memb_fit, red_memb_no_fit, bin_edges
    return red_memb_fit, red_memb_no_fit, red_plot_pars
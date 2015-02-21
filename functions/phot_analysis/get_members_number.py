# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 09:58:13 2013

@author: gabriel
"""
import numpy as np


def get_memb_num(n_clust, field_dens, clust_rad, rdp_params, bw):
    """
    Calculate the approximate number of cluster's members inside the cluster's
    radius. The process is as follows: field_dens is the density of field
    stars per unit area so multiplying that by the cluster's area we get the
    approx number of field stars inside the cluster. We can count
    the total number of stars inside the radius, n_t and then
    obtain n_c (which is the approx number of members) with:

    n_c = n_clust - [field_dens * cluster_area]
    """

    square_rings = rdp_params[-1]

    # Index of point in RDP closer to the calculated cluster radius.
    sq_indx = int(round(((clust_rad - (bw / 2.)) / bw) + 1))

    sum_bins_in_rad = 0.
    # For each square ring until reachong the limit imposed by the cluster's
    # radius.
    for sq_ring in square_rings[:sq_indx]:
        # For each bin within this square ring.
        for bin_coords in sq_ring:
            # Coordinates of bin with center in (0, 0)
            x_b, y_b = bin_coords
            # Sign of each coordinate.
            x_s = x_b / abs(float(x_b)) if x_b != 0 else 1.
            y_s = y_b / abs(float(y_b)) if y_b != 0 else 1.
            # Coordinates of farthest corner of bin.
            x_cor, y_cor = x_b + x_s / 2., y_b + y_s / 2.
            # Coordinates of corner of bin with center in (0., 0.)
            x_c, y_c = (x_cor * bw) + (bw / 2.), (y_cor * bw) + (bw / 2.)
            # Distance of center to corner of bin.
            bin_dist = np.sqrt(x_c ** 2 + y_c ** 2)
            # Length of bin diagonal.
            bin_diag = np.sqrt(2 * bw ** 2)
            if bin_dist - clust_rad > bin_diag:
                # The entire bin is outside of the cluster radius range.
                pass
            elif 0. < bin_dist - clust_rad <= bin_diag:
                # Add a portion of the bin to the total sum of bins.
                sum_bins_in_rad += min(1., (bin_dist - clust_rad) / bw)
            else:
                # Add entire bin.
                sum_bins_in_rad += 1.

    # Cluster's area.
    a_clust = (sum_bins_in_rad * bw ** 2)
    # Fraction of cluster's area present in frame.
    frac_cl_area = a_clust / (np.pi * clust_rad ** 2)

    print sum_bins_in_rad, a_clust, frac_cl_area
    raw_input()
    #import matplotlib.pyplot as plt
    #from itertools import cycle
    #cols = cycle(['red', 'blue', 'green', 'black', 'cyan'])
    #for lst in square_rings[:sq_indx]:
        #plt.scatter(zip(*lst)[0], zip(*lst)[1], c=next(cols), s=100)
    #plt.show()

    # Approx number of members.
    n_c = max(int(round(n_clust - (field_dens * a_clust))), 0)

    # Raise a flag if the number of members is < 10
    flag_num_memb_low = False if n_c >= 10 else True

    return n_c, frac_cl_area, flag_num_memb_low, a_clust, n_clust

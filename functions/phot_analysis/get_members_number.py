# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 09:58:13 2013

@author: gabriel
"""
import numpy as np


def get_memb_num(field_dens, clust_rad, rdp_params, bin_width):
    """
    Calculate the approximate number of cluster's members inside the cluster's
    radius. The process is as follows: field_dens is the density of field
    stars per unit area so multiplying that by the cluster's area we get the
    approx number of field stars inside the cluster. We can count
    the total number of stars inside the radius, n_t and then
    obtain n_c (which is the approx number of members) with:

    n_c = n_clust - [field_dens * cluster_area]
    """

    bins_in_rings = rdp_params[3]
    # Index of point in RDP closer to the calculated cluster radius.
    sq_indx = int(round(((clust_rad - (bin_width / 2.)) / bin_width) + 1))

    # Count the total number of stars inside the cluster's radius by
    # summing stars inside the bins around the center.
    n_clust = np.sum(bins_in_rings[:sq_indx], axis=0)[1]

    # Approximate cluster's area multiplying the number of bins around the
    # center up to the cluster's radius by the area of each bin.
    n_b = np.sum(bins_in_rings[:sq_indx], axis=0)[0]
    a_clust = n_b * bin_width ** 2

    # Approx number of members.
    n_c = max(int(round(n_clust - (field_dens * a_clust))), 0)

    # Raise a flag if the number of members is < 10
    flag_num_memb_low = False if n_c >= 10 else True

    return n_c, flag_num_memb_low, a_clust, n_clust

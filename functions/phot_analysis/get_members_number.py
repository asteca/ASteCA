# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 09:58:13 2013

@author: gabriel
"""


def get_memb_num(n_clust, cl_area, field_dens):
    """
    Calculate the approximate number of cluster's members inside the cluster's
    radius. The process is as follows: field_dens is the density of field
    stars per unit area so multiplying that by the cluster's area we get the
    approx number of field stars inside the cluster. We can count
    the total number of stars inside the radius, n_t and then
    obtain n_c (which is the approx number of members) with:

    n_c = n_clust - [field_dens * cluster_area]
    """

    # Approx number of members.
    n_memb = max(int(round(n_clust - (field_dens * cl_area))), 0)

    # Raise a flag if the number of members is < 10
    flag_num_memb_low = False if n_memb >= 10 else True

    if flag_num_memb_low:
        print ("  WARNING: only {:0.f} true members estimated in cluster"
            " region.".format(n_memb))
    else:
        print ("Approximate number of members in cluster obtained "
            "({:.0f}).".format(n_memb))

    return n_memb, flag_num_memb_low

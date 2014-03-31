# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 09:58:13 2013

@author: gabriel
"""
import numpy as np


def get_memb_num(backg_value, clust_rad, stars_in, stars_in_rjct):
    """
    Calculate the approximate number of cluster's members inside the cluster's
    radius. The process is as follows: backg_value is the density of field
    stars per px^2 so multiplying that by the cluster's radius we get the
    approx number of field stars inside the cluster's radius N_f. We can count
    the total number of stars inside the radius, n_t and then
    obtain n_c (which is the approx number of members) with:

    n_t = n_c + N_f
    then:
    n_c = n_t - N_f = n_t - [backg_value*(Pi*clust_rad**2)]
    """

    # Count the total number of stars inside the cluster's radius: n_t.
    n_t = len(stars_in) + len(stars_in_rjct)

    # Cluster's area.
    a_c = np.pi * (clust_rad ** 2)

    # Approx number of members.
    n_c = round(n_t - (backg_value * a_c))

    # Raise a flag if the number of members is <10
    flag_num_memb_low = False if n_c >= 10 else True

    return int(n_c), flag_num_memb_low

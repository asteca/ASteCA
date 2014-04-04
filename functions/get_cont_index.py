# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 15:51:34 2013

@author: gabriel
"""

import numpy as np


def cont_indx(backg_val, rdp_params, clust_rad, stars_in, stars_in_rjct):
    '''
    Calculate and return the contamination index value cont_index.
    The contamination index is defined as the ratio of field stars that should
    be present in the cluster region (calculated by means of the background
    value 'backg_val' and the cluster's radius 'r') to the total number of
    stars in the cluster region 'n_tot'.

    CI = [backg_val*PI*r**2]/n_tot

    If this number equals 1, it means that all of the stars in the cluster
    region are expected to be field stars. A small number (close to zero)
    means the field contamination in the cluster region is very small. A
    number larger than 1 means there are more stars in average in the
    background than there are inside the cluster region (which isn't a good
    sign).
    '''

    # Cluster's area.
    a_c = np.pi * (clust_rad ** 2)
    # Total number of stars in cluster region.
    n_tot = len(stars_in) + len(stars_in_rjct)
    # Calculate contamination index.
    cont_index = backg_val * a_c / n_tot
    print cont_index

    radii, ring_density = rdp_params[:2]
    num_stars = 0
    for indx, dens in enumerate(ring_density):
        # Count stars inside cluster region (ie: up to clust_rad)
        if radii[indx] <= clust_rad:
            # Area of previous square.
            prev_area = (2 * radii[indx - 1]) ** 2 if indx != 0 else 0.
            # Area of square ring.
            sq_rng_area = (2 * radii[indx]) ** 2 - prev_area
            # Number of stars in the square ring.
            num_stars = num_stars + (dens * sq_rng_area)

    # Star density in the cluster region.
    cl_dens = num_stars / (2 * clust_rad) ** 2

    # Final contamination index.
    cont_index = backg_val / cl_dens
    print cont_index

    return cont_index
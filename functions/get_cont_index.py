# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 15:51:34 2013

@author: gabriel
"""


def cont_indx(backg_val, rdp_params, clust_rad):
    '''
    Calculate the contamination index value. This parameter is defined as the
    ratio of field stars density over density of stars in the cluster region.

    If this number equals 1, it means that all of the stars in the cluster
    region are expected to be field stars. A small number (close to zero)
    means the field contamination in the cluster region is very small. A
    number larger than 1 means there are more stars in average in the
    background than there are inside the cluster region (which isn't a good
    sign).
    '''

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

    return cont_index
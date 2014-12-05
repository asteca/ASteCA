"""
@author: gabriel
"""

import numpy as np


def get_in_out(center_cl, clust_rad, acpt_stars, rjct_stars):
    """
    Separate stars between those inside the cluster's radius and those outside.
    """

    # Create new empty lists to store items inside and outside of the cluster's
    # radius limit.
    # cl_region will contain those stars within the radius value and
    # with accepted photometric errors.
    cl_region, stars_out, cl_region_rjct, stars_out_rjct = [], [], [], []

    # Iterate through all stars with accepted photom errors.
    for star in acpt_stars:

        # Separate in and out of cluster's boundaries.
        dist = np.sqrt((center_cl[0] - star[1]) ** 2 +
            (center_cl[1] - star[2]) ** 2)

        if dist > clust_rad:
            # Star is out of the cluster's radius limit.
            stars_out.append(star)
        else:
            # Star is inside the cluster's radius limit.
            cl_region.append(star)

    # Iterate through all stars with rejected photom errors.
    for star in rjct_stars:

        # Separate in and out of cluster's boundaries.
        dist = np.sqrt((center_cl[0] - star[1]) ** 2 +
            (center_cl[1] - star[2]) ** 2)

        if dist > clust_rad:
            # Star is out of the cluster's radius limit.
            stars_out_rjct.append(star)
        else:
            # Star is inside the cluster's radius limit.
            cl_region_rjct.append(star)

    return cl_region, stars_out, cl_region_rjct, stars_out_rjct


import numpy as np


def main(i_c, clp):
    """
    Separate stars between those inside the cluster's radius and those outside.
    """

    # cl_region will contain those stars within the radius value and
    # with accepted photometric errors.
    cl_region, stars_out, cl_region_rjct, stars_out_rjct = [], [], [], []

    # Iterate through all stars with accepted photometric errors.
    for star in clp['acpt_stars_' + i_c[0]]:

        # Separate in and out of cluster's boundaries.
        dist = np.sqrt((clp['kde_cent'][0] - star[1]) ** 2 +
                       (clp['kde_cent'][1] - star[2]) ** 2)

        if dist > clp['clust_rad']:
            # Star is outside of the cluster's radius limit.
            stars_out.append(star)
        else:
            # Star is inside the cluster's radius limit.
            cl_region.append(star)

    # Catch empty cluster region.
    if len(cl_region) <= 1:
        print("\n  ERROR: <=1 stars left in cluster region after the\n"
              "  removal of stars with large photometric errors.")
        raise ValueError('Empty cluster region.')
    elif 1 < len(cl_region) < 10:
        print("    WARNING: less than 10 stars present in cluster region.")
    else:
        print("  Stars separated in/out of cluster's boundaries.")

    # Iterate through all stars with rejected photometric errors.
    for star in clp['rjct_stars_' + i_c[0]]:

        # Separate in and out of cluster's boundaries.
        dist = np.sqrt((clp['kde_cent'][0] - star[1]) ** 2 +
                       (clp['kde_cent'][1] - star[2]) ** 2)

        if dist > clp['clust_rad']:
            # Star is out of the cluster's radius limit.
            stars_out_rjct.append(star)
        else:
            # Star is inside the cluster's radius limit.
            cl_region_rjct.append(star)

    # Add parameters to dictionary.
    clp['cl_region_' + i_c[0]], clp['stars_out_' + i_c[0]],\
        clp['cl_region_rjct_' + i_c[0]], clp['stars_out_rjct_' + i_c[0]] =\
        cl_region, stars_out, cl_region_rjct, stars_out_rjct
    return clp

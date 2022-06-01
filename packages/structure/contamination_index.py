
import numpy as np
from scipy.spatial.distance import cdist


def main(clp, x, y, mags, **kwargs):
    """
    Calculate the contamination index value. This parameter is defined as the
    ratio of field stars density over the density of stars in the cluster
    region.

    CI ~ 0.0 -> the field contamination in the cluster region is very small.
    CI ~ 0.5 -> an equal number of field stars and cluster members are
    expected inside the cluster region.
    CI ~ 1.0 -> there are no expected cluster members inside the cluster region
    (not a good sign).
    """

    # If the cluster's area is less than 10% of the full frame's area, don't
    # estimate the CI.
    frame_area = np.ptp(x) * np.ptp(y)
    if (frame_area - clp['cl_area']) / frame_area > .1:

        # Count the total number of stars within the defined cluster region
        # (including stars with rejected photometric errors)
        dist = cdist([clp['kde_cent']], np.array([x, y]).T)[0]
        n_in_cl_reg = (dist < clp['clust_rad']).sum()

        # Final contamination index and estimated number of members.
        cont_index = CIfunc(n_in_cl_reg, clp['field_dens'], clp['cl_area'])

        if cont_index >= 1.:
            print("  WARNING: contamination index value is very large: "
                  "{:.2f}".format(cont_index))
        else:
            print("Contamination index obtained ({:.2f})".format(cont_index))
    else:
        print("  WARNING: cluster radius is too large to obtain\n"
              "  a reliable contamination index value")
        cont_index = 0.

    clp['cont_index'] = cont_index
    return clp


def CIfunc(n_in_cl_reg, field_dens, area):
    """
    Contamination index
    """
    # Star density in the cluster region.
    cl_dens = n_in_cl_reg / area
    CI = field_dens / cl_dens

    return CI


import numpy as np
from scipy.spatial.distance import cdist


def main(clp, x, y, **kwargs):
    """
    Calculate the contamination index value. This parameter is defined as the
    ratio of field stars density over the density of stars in the cluster
    region. Uses the 'incomplete' data.

    A small number (close to zero) means the field contamination in the
    cluster region is very small.
    If this number equals 0.5, it means that an equal number of field stars
    and cluster members are expected inside the cluster region. A value of
    1 means there are no expected cluster members inside the cluster region
    (which isn't a good sign).
    """

    # If the cluster radius exceeds the length of the area where the field
    # density value was obtained (ie: the extension of the RDP), then do not
    # obtain the 'cont_index' parameter since the field density does not
    # represent the density of the field but rather the density of the
    # outermost regions of the cluster.
    if clp['clust_rad'] < clp['rdp_radii'][-1]:

        # Count the total number of stars within the defined cluster region
        # (including stars with rejected photometric errors)
        dist = cdist([clp['kde_cent']], np.array([x, y]).T)[0]
        # cdist(np.array([x, y]).T, np.atleast_2d(clp['kde_cent']))
        n_in_cl_reg = (dist < clp['clust_rad']).sum()

        # Final contamination index.
        cont_index, n_memb_i, _ = CIfunc(
            n_in_cl_reg, clp['field_dens'], clp['cl_area'])
        n_memb_i = int(n_memb_i)

        if cont_index >= 1.:
            print("  WARNING: contamination index value is very large: "
                  "{:.2f}".format(cont_index))
        else:
            print("Contamination index obtained ({:.2f})".format(cont_index))
    else:
        print("  WARNING: cluster radius is too large to obtain\n"
              "  a reliable contamination index value")
        cont_index, n_memb_i = np.nan, np.nan

    clp['cont_index'], clp['n_memb_i'] = cont_index, n_memb_i
    return clp


def CIfunc(n_in_cl_reg, field_dens, area):
    """
    """
    # Estimated number of field stars in the area.
    n_fl = field_dens * area
    # Estimated number of members in the area.
    n_memb = n_in_cl_reg - n_fl

    # Star density in the cluster region.
    cl_dens = n_in_cl_reg / area
    # Contamination index.
    CI = field_dens / cl_dens

    return CI, n_memb, n_fl

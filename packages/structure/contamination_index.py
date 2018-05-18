
import numpy as np


def main(clp, x, y, **kwargs):
    '''
    Calculate the contamination index value. This parameter is defined as the
    ratio of field stars density over the density of stars in the cluster
    region.

    A small number (close to zero) means the field contamination in the
    cluster region is very small.
    If this number equals 0.5, it means that an equal number of field stars
    and cluster members are expected inside the cluster region. A value of
    1 means there are no expected cluster members inside the cluster region
    (which isn't a good sign).
    '''

    # If the cluster radius exceeds the length of the area where the field
    # density value was obtained (ie: the extension of the RDP), then do not
    # obtain the 'cont_index' parameter since the field density does not
    # represent the density of the field but rather the density of the
    # outermost regions of the cluster.
    if clp['clust_rad'] < clp['rdp_length'] / 2.:

        # Count the total number of stars within the defined cluster region
        # (including stars with rejected photometric errors)
        n_clust = 0
        for xy in zip(*[x, y]):
            # Separate in and out of cluster's boundaries.
            dist = np.sqrt((clp['kde_cent'][0] - xy[0]) ** 2 +
                           (clp['kde_cent'][1] - xy[1]) ** 2)
            if dist <= clp['clust_rad']:
                n_clust += 1

        # Star density in the cluster region.
        cl_dens = n_clust / clp['cl_area']

        # Final contamination index.
        cont_index = clp['field_dens'] / cl_dens

        if cont_index >= 1.:
            print("  WARNING: contamination index value is very large: "
                  "{:.2f}".format(cont_index))
        else:
            print('Contamination index obtained ({:.2f}).'.format(cont_index))
    else:
        print("  WARNING: cluster radius is too large to obtain\n"
              "  a reliable contamination index value.")
        cont_index = -1.

    clp['cont_index'] = cont_index
    return clp

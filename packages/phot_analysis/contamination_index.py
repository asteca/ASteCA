

def main(n_clust, cl_area, field_dens, clust_rad, rdp_length):
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
    # obtain the n_memb parameter since the field density does not represent
    # the density of the field but rather the density of the outermost regions
    # of the cluster.
    if clust_rad < rdp_length / 2.:

        # Star density in the cluster region.
        cl_dens = n_clust / cl_area

        # Final contamination index.
        cont_index = field_dens / cl_dens

        if cont_index >= 1.:
            print("  WARNING: CI value obtained is too high: "
                  "{:.2f}".format(cont_index))
        else:
            print('Contamination index obtained ({:.2f}).'.format(cont_index))
    else:
        print("  WARNING: cluster radius is too large to obtain\n"
              "  a reliable contamination index value.")
        cont_index = -1.

    return cont_index

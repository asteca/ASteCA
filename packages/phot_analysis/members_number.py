

def main(clp):
    """
    Calculate the approximate number of cluster's members inside the cluster's
    radius as follows:

    n_c = n_clust - [field_dens * cluster_area]
    """

    # If the cluster radius exceeds the length of the area where the field
    # density value was obtained (ie: the extension of the RDP), then do not
    # obtain the n_memb parameter since the field density does not represent
    # the density of the field but rather the density of the outermost regions
    # of the cluster.
    if clp['clust_rad'] < clp['rdp_length'] / 2.:

        # Approx number of members.
        n_memb = max(int(round(
            clp['n_clust'] - (clp['field_dens'] * clp['cl_area']))), 0)

        # Raise a flag if the number of members is < 10
        flag_num_memb_low = False if n_memb >= 10 else True

        if flag_num_memb_low:
            print("  WARNING: only {:.0f} true members estimated in cluster"
                  " region.".format(n_memb))
        else:
            print("Approximate number of members in cluster obtained "
                  "({:.0f}).".format(n_memb))
    else:
        print("  WARNING: cluster radius is too large to obtain\n"
              "  the approximate number of members.")
        n_memb, flag_num_memb_low = -1., True

    clp['n_memb'], clp['flag_num_memb_low'] = n_memb, flag_num_memb_low
    return clp



def main(n_clust, cl_area, field_dens, clust_rad, rdp_length):
    """
    Calculate the approximate number of cluster's members inside the cluster's
    radius. The process is as follows: field_dens is the density of field
    stars per unit area so multiplying that by the cluster's area we get the
    approx number of field stars inside the cluster. We can count
    the total number of stars inside the radius, n_t and then
    obtain n_c (which is the approx number of members) with:

    n_c = n_clust - [field_dens * cluster_area]
    """

    # If the cluster radius exceeds the length of the area where the field
    # density value was obtained (ie: the extension of the RDP), then do not
    # obtain the n_memb parameter since the field density does not represent
    # the density of the field but rather the density of the outermost regions
    # of the cluster.
    if clust_rad < rdp_length / 2.:

        # Approx number of members.
        n_memb = max(int(round(n_clust - (field_dens * cl_area))), 0)

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

    return n_memb, flag_num_memb_low


import numpy as np

# DEPRECATED 05/22

def main(clp):
    """
    Calculate the approximate number of cluster's members inside the cluster's
    radius as follows:

    n_c = N_cluster_region - mean(N_field_regions)
    """

    # If no field regions were defined, this parameter can not be obtained.
    if clp['flag_no_fl_regs']:

        print("  WARNING: no field regions defined, can not estimate\n"
              "  the number of cluster members for the complete set")
        n_memb = np.nan

    else:
        # Approx number of members.
        n_memb = max(int(round(
            len(clp['cl_region'])
            - np.mean([len(_) for _ in clp['field_regions']]))), 0)

        if n_memb < 10:
            print("  WARNING: only {:.0f} true members estimated in cluster"
                  "  region".format(n_memb))

    clp['n_memb'] = n_memb
    return clp

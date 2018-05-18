
import numpy as np


def main(clp):
    """
    Calculate the approximate number of cluster's members inside the cluster's
    radius as follows:

    n_c = N_cluster_region - mean(N_field_regions)
    """

    # If no field regions were defined, this parameter can not be obtained.
    if clp['flag_no_fl_regs_c']:

        print("  WARNING: no field regions defined, can not estimate\n"
              "  the approximate number of cluster members.")
        n_memb, flag_num_memb_low = -1., True

    else:
        # Approx number of members.
        n_memb = max(int(round(
            len(clp['cl_region_c']) -
            np.mean([len(_) for _ in clp['field_regions_c']]))), 0)

        # Raise a flag if the number of members is < 10
        flag_num_memb_low = False if n_memb >= 10 else True

        if flag_num_memb_low:
            print("  WARNING: only {:.0f} true members estimated in cluster"
                  " region.".format(n_memb))
        else:
            print("Approximate number of members in cluster obtained "
                  "({:.0f}).".format(n_memb))

    clp['n_memb'], clp['flag_num_memb_low'] = n_memb, flag_num_memb_low
    return clp

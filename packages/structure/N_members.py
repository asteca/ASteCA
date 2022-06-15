
import numpy as np


def main(clp):
    """
    Distribution for the number of members estimated. Uses the uncertainties
    in the field density.

    The estimated number of members 'n_memb' is used by the King's profile
    function.
    """

    N_in_rad = (clp['xy_cent_dist'] < clp['clust_rad']).sum()

    n_memb = 0
    # The field density was manually fixed
    if np.isnan(clp['field_dens_std']):
        n_memb = int(N_in_rad - clp['cl_area'] * clp['field_dens'])
        N_memb_all = np.array([])
    else:
        field_dens_s = np.random.normal(
            clp['field_dens'], clp['field_dens_std'], 1000)
        N_memb_all = N_in_rad - clp['cl_area'] * field_dens_s
        N_memb_all = N_memb_all[N_memb_all > 0]
        if N_memb_all.size > 0:
            n_memb = int(np.median(N_memb_all))

    if n_memb < 10:
        print("  WARNING: The estimated number of members is < 10")
        n_memb = 10

    clp['n_memb'], clp['members_dist'] = n_memb, N_memb_all

    return clp

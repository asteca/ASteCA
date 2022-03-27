
import numpy as np


def main(clp):
    """
    Distribution for the number of members estimated. Uses the uncertainties
    in the field density.

    The estimated number of members 'n_memb' is used by the King's profile
    function.
    """

    # Field density and radius were st manually
    if np.isnan(clp['field_dens_std']):
        clp['n_memb'], clp['members_dist'] = 1, np.array([])
        return clp

    N_in_rad = (clp['xy_cent_dist'] < clp['clust_rad']).sum()
    field_dens_s = np.random.normal(
        clp['field_dens'], clp['field_dens_std'], 1000)
    N_memb_all = N_in_rad - clp['cl_area'] * field_dens_s
    N_memb_all = N_memb_all[N_memb_all > 0]

    # Used in the King profile fitting
    if len(N_memb_all) > 10:
        n_memb_i = int(np.median(N_memb_all))
    else:
        print("  WARNING: The estimated number of members is < 10")
        n_memb_i = 10
    clp['n_memb'], clp['members_dist'] = n_memb_i, N_memb_all

    return clp

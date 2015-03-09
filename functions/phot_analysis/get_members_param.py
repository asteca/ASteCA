# -*- coding: utf-8 -*-
"""
Created on Mon Mar 09 12:00:00 2015

@author: gabriel
"""


def mp_members(n_memb, decont_algor_return):
    '''
    Compare the estimated number of true members obtained via the stars density
    analysis done in `get_members_number` with the number of stars in the
    cluster region that are assigned a MP of 0.5 or more. These stars are the
    ones with a greater probability of being cluster members than field region
    stars.
    '''

    memb_prob_avrg_sort, flag_decont_skip = decont_algor_return

    memb_par, n_memb_da = float("inf"), -1
    # Obtain parameter if the DA was applied.
    if not flag_decont_skip:

        n_memb_da = 0
        # Number of stars assigned a MP>=0.5.
        for star in memb_prob_avrg_sort:
            if star[7] >= 0.5:
                n_memb_da += 1

        if n_memb != 0 or n_memb_da != 0:
            # Obtain parameter.
            memb_par = (float(n_memb) - float(n_memb_da)) / \
                (float(n_memb) + float(n_memb_da))

    return memb_par, n_memb_da

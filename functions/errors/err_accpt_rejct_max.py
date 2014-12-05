# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 17:41:28 2013

@author: gabriel
"""


def err_a_r_m(e_mag, e_col1, err_pck):
    """
    Accept stars with photom errors < e_max both in mag and in color.
    """

    # Unpack params.
    e_max = err_pck[0][1]

    # Initialize empty list to hold accepted/rejected stars' indexes.
    acpt_indx, rjct_indx = [], []

    # Iterate through all stars
    for st_ind in range(len(e_mag)):

        # Reject stars with at least one error >= e_max.
        if e_mag[st_ind] >= e_max or e_col1[st_ind] >= e_max:
            rjct_indx.append(st_ind)
        else:
            # Accept star.
            acpt_indx.append(st_ind)

    # Pass empty list.
    err_plot = []
    return acpt_indx, rjct_indx, err_plot
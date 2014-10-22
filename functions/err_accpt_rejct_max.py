# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 17:41:28 2013

@author: gabriel
"""
import get_in_params as g


def err_a_r_m(e_mag, e_col):
    """
    Accept stars with photom errors < e_max both in mag and in color.
    """

    # Unpack e_max.
    e_max = g.er_params[1]

    # Initialize empty list to hold accepted/rejected stars' indexes.
    acpt_indx, rjct_indx = [], []

    # Create list of combined errors for all stars.
    em_z, ec_z = zip(*e_mag), zip(*e_col)
    e_mc = [em_z[i] + ec_z[i] for i in range(len(em_z))]

    # Iterate through all stars
    for st_ind in range(len(e_mag[0])):

        # Reject stars with at least one error >= e_max.
        if any(e >= e_max for e in e_mc[st_ind]):
            rjct_indx.append(st_ind)
        else:
            # Accept star.
            acpt_indx.append(st_ind)

    # Pass empty list.
    err_plot = []
    return acpt_indx, rjct_indx, err_plot
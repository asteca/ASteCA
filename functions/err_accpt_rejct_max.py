# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 17:41:28 2013

@author: gabriel
"""

def err_a_r_m(phot_data, er_params):
    """    
    Accept stars with photom errors < e_max both in mag and in color.
    """

    id_star, x_data, y_data, mag, e_mag, col1, e_col1 = phot_data
    e_max = er_params[2]

    # Initialize empty list to hold accepted/rejected stars.
    acpt_stars, rjct_stars = [], []    
    
    # Iterate through all stars
    for st_ind, star_id in enumerate(id_star):

        # Reject stars with at least one error >= e_max.
        if e_mag[st_ind] >= e_max or e_col1[st_ind] >= e_max:

            rjct_stars.append([star_id, x_data[st_ind], y_data[st_ind],
                               mag[st_ind], e_mag[st_ind], col1[st_ind], \
                               e_col1[st_ind]])
        
        else:
            # Accept star.
            acpt_stars.append([star_id, x_data[st_ind], y_data[st_ind],
                               mag[st_ind], e_mag[st_ind], col1[st_ind], \
                               e_col1[st_ind]])
                               
    return acpt_stars, rjct_stars
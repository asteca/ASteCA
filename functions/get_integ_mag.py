# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 20:55:57 2013

@author: gabriel
"""

import numpy as np

def calc_integ_mag(mag_list, mag_limit):
    '''
    Calculate integrated magnitude up to a certain maximum magnitude value.
    '''
    sort_lis = sorted(mag_list)
    trim_lis = []
    for i in sort_lis:
        if i <= mag_limit:
            trim_lis.append(i)
    n_st = float(len(trim_lis))
    
    if n_st == 0.:
        int_mag_val = mag_limit
    else:
        sum_mag = sum(trim_lis)
        sum_log_energ = 0.
        for mag_j in trim_lis:
            energ_sum = 0.
            for mag_i in trim_lis:
                energ_sum = energ_sum + 10**((mag_i-mag_j)/-2.5)
            sum_log_energ = sum_log_energ + np.log10(energ_sum)
        int_mag_val = (1/n_st)*(-2.5*sum_log_energ + sum_mag)
    
    return int_mag_val
    

def integ_mag(stars_in, stars_in_rjct):
    '''
    Obtain integrated magnitude using all stars inside the cluster's radius for
    several limits in magnitude.
    '''
    # Generate lists holding only mag values.
    mags_in = [i[3] for i in stars_in]
    mags_in_all = mags_in + [i[3] for i in stars_in_rjct]

    # Define range in magnitude.
    min_mag, max_mag = min(min(mags_in), min(mags_in_all)), max(max(mags_in),\
                           max(mags_in_all))
    mag_range = np.arange(min_mag, max_mag, 0.1)
    
    # Define final output lists.
    stars_in_mag, stars_in_all_mag = [[],[]], [[],[]]
    # Combine into one list.
    all_lists = [stars_in_mag, stars_in_all_mag]
                 
    # Loop through each list.
    for indx, out_list in enumerate(all_lists):

        # Define the list to be used.        
        if indx == 0:
            mag_list = mags_in
        elif indx == 1:
            mag_list = mags_in_all

        for mag_limit in mag_range:
            # Calculate integrated magnitude up to this mag value.
            int_mag_val = calc_integ_mag(mag_list, mag_limit)
            out_list[0].append(mag_limit)
            out_list[1].append(int_mag_val)
    
    return stars_in_mag, stars_in_all_mag

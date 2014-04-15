# -*- coding: utf-8 -*-
"""
Created on Thu Mar 20 2014

@author: gabriel
"""


def red_memb(flag_area_stronger, memb_prob_avrg_sort, rm_params):
    '''
    Reduce number of stars according to a given membership probability
    lower limit.
    '''

    flag_red_memb, min_prob = rm_params

    memb_prob_avrg_sort2 = []

    if flag_area_stronger:
        memb_prob_avrg_sort2 = memb_prob_avrg_sort
        print "WARNING: no field regions found. Skipping membership reduction."

    else:

        if flag_red_memb == 'auto':
            # Reject stars in the lower half of the membership probabilities
            # list.
            middle_indx = int(len(memb_prob_avrg_sort) / 2)
            memb_prob_avrg_sort2 = memb_prob_avrg_sort[:middle_indx]

        elif flag_red_memb == 'manual':

            for star in memb_prob_avrg_sort:
                if star[7] >= min_prob:
                    memb_prob_avrg_sort2.append(star)

        else:
            memb_prob_avrg_sort2 = memb_prob_avrg_sort


    return memb_prob_avrg_sort2
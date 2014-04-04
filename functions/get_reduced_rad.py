# -*- coding: utf-8 -*-
"""
Created on Thu Mar 20 2014

@author: gabriel
"""

from get_cont_index import cont_indx as ci
import numpy as np


def red_rad(flag_red_rad, backg_val, clust_rad, cont_index, center_cl,
    rdp_params, memb_prob_avrg_sort):
    '''
    Attempts to reduce the radius until a value of CI<=0.5 is
    achieved.
    '''

    radii, ring_density = rdp_params[:2]

    new_rad, flag_new_rad = clust_rad, False
    i = 1
    if flag_red_rad is True:

        stars_in = [0. for _ in memb_prob_avrg_sort]
        while cont_index > 0.5 and len(stars_in) > \
        0.05 * len(memb_prob_avrg_sort):
            # Reduce radius by 1 px.
            new_rad = clust_rad - i

            stars_in, stars_in_rjct = [], []
            for star in memb_prob_avrg_sort:
                dist = np.sqrt((center_cl[0] - star[1]) ** 2 +
                    (center_cl[1] - star[2]) ** 2)
                if dist <= new_rad:
                    stars_in.append(0.)

            print new_rad, backg_val, len(stars_in), np.pi * (new_rad ** 2), \
            (float(len(stars_in)) / np.pi * (new_rad ** 2))

            # Get new cont_indx with new radius value.
            cont_index = ci(backg_val, new_rad, stars_in, stars_in_rjct)
            if cont_index <= 0.5:
                flag_new_rad = True
            i += 1
            print cont_index

        if flag_new_rad is True:
            memb_prob_avrg_sort2 = []
            for star in memb_prob_avrg_sort:
                dist = np.sqrt((center_cl[0] - star[1]) ** 2 +
                    (center_cl[1] - star[2]) ** 2)
                if dist <= new_rad:
                    memb_prob_avrg_sort2.append(star)

    else:
        memb_prob_avrg_sort2 = memb_prob_avrg_sort

    print len(memb_prob_avrg_sort), len(memb_prob_avrg_sort2)
    raw_input()
    return memb_prob_avrg_sort2
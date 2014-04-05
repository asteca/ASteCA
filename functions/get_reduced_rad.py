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

    flag_new_rad = False
    indx_rad = min(range(len(radii)), key=lambda i: abs(radii[i] - clust_rad))
    new_rad = radii[indx_rad]

    i = 1
    if flag_red_rad is True and cont_index > 0.5:

        # The condition (indx_rad - i) > 0 means that the minimum radius
        # value is the one associated with the second point in the raddi
        # list.
        while flag_new_rad is False and (indx_rad - i) > 0:

            # Reduce radius.
            new_rad = radii[indx_rad - i]
            i += 1
            # Get new cont_indx with new radius value.
            cont_index = ci(backg_val, rdp_params, new_rad)

            if cont_index <= 0.5:
                # Obtain new members.
                memb_prob_avrg_sort2 = []
                for star in memb_prob_avrg_sort:
                    dist = np.sqrt((center_cl[0] - star[1]) ** 2 +
                        (center_cl[1] - star[2]) ** 2)
                    if dist <= new_rad:
                        memb_prob_avrg_sort2.append(star)
                # Break out of while.
                flag_new_rad = True
                print 'Reduced radius found (%0.2f).' % new_rad

    if flag_red_rad is not True or flag_new_rad is not True:
        print 'No reduced radius value found.'
        # Return unchanged list of members.
        memb_prob_avrg_sort2 = memb_prob_avrg_sort

    return memb_prob_avrg_sort2
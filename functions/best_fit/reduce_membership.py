# -*- coding: utf-8 -*-
"""
Created on Thu Mar 20 2014

@author: gabriel
"""

from .._in import get_in_params as g


def auto_select(n_memb, memb_prob_avrg_sort):
    '''
    Auto algorithm to select which stars to use by the best fit funtcion.
    Will set the minimum probability value such that an equal number of
    stars are used in the best fit process, as the approximate number of
    members found when comparing the density of the cluster region with that
    of the field regions defined.
    '''
    red_memb_fit, red_memb_no_fit, min_prob = memb_prob_avrg_sort, [], -1.

    # Check approximate number of true members obtained by the structural
    # analysis.
    if n_memb < 10:
        print ("  WARNING: less than 10 stars identified as true\n"
        "  cluster members. Using full list.")
    else:

        # Total number of stars in the cluster region.
        n_tot = len(memb_prob_avrg_sort)

        # If there are less stars in the cluster region than than n_memb stars,
        # use all stars in the cluster region.
        if n_memb >= n_tot:
            # Use all stars in the cluster region.
            indx, min_prob = n_tot, 0.
        else:
            # Use the first n_memb stars, ie: those stars with the highest
            # membership probability values.
            indx, min_prob = n_memb, zip(*memb_prob_avrg_sort)[-1][n_memb]

    # DEPRECATED
    #elif n_memb_da >= n_memb:
        ## There are more stars with MP>=0.5 than true members estimated by
        ## the structural analysis. Use all stars with MP>=0.5.
        #indx, min_prob = 0, 0.5
        #for prob in zip(*memb_prob_avrg_sort)[-1]:
            #if prob >= 0.5:
                #indx += 1
            #else:
                #break

        red_memb_fit, red_memb_no_fit = memb_prob_avrg_sort[:indx], \
        memb_prob_avrg_sort[indx:]

    return red_memb_fit, red_memb_no_fit, min_prob


def red_memb(n_memb, decont_algor_return):
    '''
    Reduce number of stars according to a given membership probability
    lower limit or minimum magnitude limit.
    '''

    bf_flag = g.bf_params[0]
    memb_prob_avrg_sort, flag_decont_skip = decont_algor_return

    # Skip reduction process.
    red_memb_fit, red_memb_no_fit, min_prob = memb_prob_avrg_sort, [], -1.
    # Only run if best fit process is set to run.
    if bf_flag:

        mode_red_memb, min_prob_man = g.rm_params

        if flag_decont_skip is True:
            print ("  WARNING: decontamination algorithm was skipped.\n"
            "  Can't apply membership reduction. Using full list.")

        elif mode_red_memb == 'auto':
            # Select stars to use automatically.
            red_memb_fit, red_memb_no_fit, min_prob = auto_select(n_memb,
                memb_prob_avrg_sort)

        elif mode_red_memb == 'top-h':
            # Reject stars in the lower half of the membership
            # probabilities list.
            middle_indx = int(len(memb_prob_avrg_sort) / 2)
            red_fit = memb_prob_avrg_sort[:middle_indx]
            # Check number of stars left.
            if len(red_fit) > 10:
                red_memb_fit, red_memb_no_fit, min_prob = red_fit, \
                memb_prob_avrg_sort[middle_indx:], \
                memb_prob_avrg_sort[middle_indx][-1]
            else:
                print ("  WARNING: less than 10 stars left after reducing\n"
                "  by top half membership probability. Using full list.")

        elif mode_red_memb == 'man':
            # Find index of star with membership probability < min_prob_man
            indx = 0
            for star in memb_prob_avrg_sort:
                if star[-1] < min_prob_man:
                    break
                else:
                    indx += 1

            if len(memb_prob_avrg_sort[:indx]) > 10:
                red_memb_fit, red_memb_no_fit, min_prob = \
                memb_prob_avrg_sort[:indx], memb_prob_avrg_sort[indx:],\
                min_prob_man
            else:
                print ("  WARNING: less than 10 stars left after reducing\n"
                "  by manual membership probability. Using full list.")

        elif mode_red_memb == 'mag':
            # Reject stars beyond the given magnitude limit.
            red_fit, red_not_fit = [], []
            for star in memb_prob_avrg_sort:
                if star[3] <= min_prob_man:
                    red_fit.append(star)
                else:
                    red_not_fit.append(star)

            # Check number of stars left.
            if len(red_fit) > 10:
                red_memb_fit, red_memb_no_fit, min_prob = red_fit, red_not_fit,\
                min_prob_man
            else:
                print ("  WARNING: less than 10 stars left after reducing\n"
                "  by magnitude limit. Using full list.")

    return red_memb_fit, red_memb_no_fit, min_prob
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 20 2014

@author: gabriel
"""

from scipy import stats


def memb_mu(memb_prob_avrg_sort, flag_memb_mu):
    '''
    Use only those stars in the upper half of the gaussian fit to the
    distribution of membership probabilities.
    '''

    prob_data = [star[7] for star in memb_prob_avrg_sort]
    # Best Gaussian fit of data.
    (mu, sigma) = stats.norm.fit(prob_data)

    memb_prob_avrg_sort2 = []
    if flag_memb_mu is True:

        for star in memb_prob_avrg_sort:
            if star[7] >= mu:
                memb_prob_avrg_sort2.append(star)

    else:
        memb_prob_avrg_sort2 = memb_prob_avrg_sort

    return memb_prob_avrg_sort2
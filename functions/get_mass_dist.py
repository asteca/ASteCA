# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 15:26:39 2014

@author: gabriel
"""

# http://www.astro.ru.nl/~slarsen/teaching/Galaxies/cmd.pdf
# http://python4mpia.github.io/fitting_data/MC-sampling-from-Salpeter.html

import numpy as np


def weighted_fast_samp(refs, weights, lim):

    # Sample the masses with replacement according to the probabilities
    # given to each by the PDF.
    blocksize = 10000
    samp = np.random.choice(refs, size=blocksize, replace=True, p=weights)
    samp_sum = np.cumsum(samp)

    # Is the sum of our current block of samples >= lim?
    while samp_sum[-1] < lim:

        # Sample another block.
        newsamp = np.random.choice(refs, size=blocksize, replace=True,
                                   p=weights)
        # Stack arrays in sequence horizontally (column wise).
        samp = np.hstack((samp, newsamp))
        # Update sum of masses.
        samp_sum = np.hstack((samp_sum, np.cumsum(newsamp) + samp_sum[-1]))

    # Find index where the lim value is found in samp_sum.
    last = np.searchsorted(samp_sum, lim, side='right')

    return samp[:last + 1]


def mass_dist(mass_params):
    '''
    Returns a mass distribution according to a given IMF and a total cluster
    mass.
    '''
    imf_cdf, M_total = mass_params
    # Sample masses from IMF according to the PDF obtained from it until
    # the sum of all stars reaches M_total.
    dist_mass = weighted_fast_samp(imf_cdf[0], imf_cdf[1], M_total)

    return dist_mass
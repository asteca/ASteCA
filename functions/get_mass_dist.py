# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 15:26:39 2014

@author: gabriel
"""

# http://www.astro.ru.nl/~slarsen/teaching/Galaxies/cmd.pdf
# http://python4mpia.github.io/fitting_data/MC-sampling-from-Salpeter.html

import numpy as np


def weighted_fast_samp(mass_init, probs, M_total):

    # Sample the masses with replacement according to the probabilities
    # given to each by the PDF.
    blocksize = 10000
    masses = np.random.choice(mass_init, size=blocksize, replace=True, p=probs)

    #from scipy.stats import rv_discrete
    #indexes = np.arange(len(mass_init))
    #random_masses = rv_discrete(values=(indexes, probs))
    #masses_indx = random_masses.rvs(size=blocksize)
    #masses = mass_init[masses_indx]

    # Sum all the random masses selected (CDF)
    masses_sum = np.cumsum(masses)

    # Is the sum of our current block of masses >= M_total?
    while masses_sum[-1] < M_total:

        # Sample another block.
        newsamp = np.random.choice(mass_init, size=blocksize, replace=True,
                                   p=probs)
        # Stack arrays in sequence horizontally (column wise).
        masses = np.hstack((masses, newsamp))
        # Update sum of masses and go check again.
        masses_sum = np.hstack((masses_sum, np.cumsum(newsamp) +
            masses_sum[-1]))

    # Find index where the M_total value is found in samp_sum.
    last = np.searchsorted(masses_sum, M_total, side='right')
    masses_pass = masses[:last + 1]

    return masses_pass


def mass_dist(imf_pdf, M_total):
    '''
    Returns a mass distribution according to a given IMF and a total cluster
    mass.
    '''
    # Sample masses from IMF according to the PDF obtained from it until
    # the sum of all stars reaches M_total.
    dist_mass = weighted_fast_samp(imf_pdf[0], imf_pdf[1], M_total)

    return dist_mass
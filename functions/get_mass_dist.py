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
    blocksize=10000
    samp = np.random.choice(refs, size=blocksize, replace=True, p=weights)
    samp_sum = np.cumsum(samp)

    # Is the sum of our current block of samples >= lim?
    while samp_sum[-1] < lim:

        # if not, we'll sample another block and try again until it is
        newsamp = np.random.choice(refs, size=blocksize, replace=True, 
                                   p=weights)
        samp = np.hstack((samp, newsamp))
        samp_sum = np.hstack((samp_sum, np.cumsum(newsamp) +  samp_sum[-1]))

    last = np.searchsorted(samp_sum, lim, side='right')
    
    return samp[:last + 1]
    

def mass_dist(mass_params):
    '''
    Main function.
    
    Returns a mass distribution according to a given IMF and a total cluster
    mass.
    '''
    imf_cdf, M_total = mass_params
    
    dist_mass = weighted_fast_samp(imf_cdf[0], imf_cdf[1], M_total)
    
    return dist_mass
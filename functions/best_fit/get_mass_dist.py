# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 15:26:39 2014

@author: gabriel
"""

# http://www.astro.ru.nl/~slarsen/teaching/Galaxies/cmd.pdf
# http://python4mpia.github.io/fitting_data/MC-sampling-from-Salpeter.html
import numpy as np


def mass_dist(st_dist, M_total):
    '''
    Returns a mass distribution according to a given IMF and a total cluster
    mass.
    '''
    # Generate a number N of stars for each interval (m, m+dm) with masses
    # randomly distributed within the interval.

    # Normalize number of stars per interval  of mass according to total mass.
    mass_up, N_stars = st_dist[0], st_dist[1] * M_total

    # Generate stars in each interval.
    dist_mass = []
    m_low, N_st_add = 0.01, 0.
    for m_up, N_st in zip(*[mass_up, N_stars]):
        N_st += N_st_add
        # If the number of stars is less than 1, add as many intervals as
        # necessary until it is at least one and then generate that star(s)
        # with a random mass in the m_low, m_up interval.
        if N_st < 1.:
            # Store this fraction to be added with the next one.
            N_st_add = N_st
        else:
            # Generate N_st stars with masses randomly distributed between
            # m_low and m_up and store them in the dist_mass list.
            dist_mass.extend(np.random.uniform(m_low, m_up, int(round(N_st))))
            # Reset parameters.
            N_st_add, m_low = 0., m_up

    return dist_mass
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 15:26:39 2014

@author: gabriel
"""

# http://www.astro.ru.nl/~slarsen/teaching/Galaxies/cmd.pdf
# http://python4mpia.github.io/fitting_data/MC-sampling-from-Salpeter.html

import numpy as np
#import time


def get_mass_cdf(imf_cdf):
    '''
    Draw a random mass value from the CDF.
    '''
    # Draw random [0,1) value.
    u = np.random.random()
    # Find closest value in CDF array.
    idx = np.argmin(np.abs(u - np.asarray(imf_cdf[1])))
    # Identify mass value corresponding to this CDF value.
    m_upper = imf_cdf[0][idx]

    return m_upper


def mass_dist(mass_params):
    '''
    Main function.
    
    Returns a mass distribution according to a given IMF and
    a total number of stars or total cluster mass.
    '''
    
    imf_cdf, limit_sel, MN_total = mass_params
    
    # Prepare array for output masses.
#    tik=time.time()
    dist_mass = []
    # Fill in array according to selected limits.
    if limit_sel == 'total_number':
        while (len(dist_mass) < MN_total):
            # Draw random mass from the IMF's CDF.
            m_upper = get_mass_cdf(imf_cdf)
            # Store mass in array.
            dist_mass.append(round(m_upper,3))
        
    elif limit_sel == 'total_mass':
        while (sum(dist_mass) < MN_total):
            # Draw random mass from the IMF's CDF.
            m_upper = get_mass_cdf(imf_cdf)
            # Store mass in array.
            dist_mass.append(round(m_upper,3))
            
#    print '0', time.time()-tik
#    raw_input()    
    
    return dist_mass
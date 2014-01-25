# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 15:26:39 2014

@author: gabriel
"""

# http://www.astro.ru.nl/~slarsen/teaching/Galaxies/cmd.pdf
# http://python4mpia.github.io/fitting_data/MC-sampling-from-Salpeter.html

import numpy as np
from scipy.integrate import quad


def imfs(imf_name, m_star, norm_const):
    '''
    Define any number of IMFs.
    '''
    if imf_name == 'kroupa_1993':
        # Kroupa, Tout & Gilmore. (1993) piecewise IMF.
        power = [-1.3, -2.2, -2.7]
        if 0.08<m_star<=0.5:
            i = 0
        elif 0.5<m_star<=1.:
            i = 1
        elif 1.<m_star:
            i = 2
        imf_val = norm_const*m_star**power[i]
        
    elif imf_name == 'chabrier_2001':
        # Chabrier (2001) exponential form of the IMF.
        imf_val = norm_const*3.*m_star**(-3.3)*np.exp(-(716.4/m_star)**0.25)

    return imf_val
    

def integral_IMF_M(m_star, imf_sel, norm_const):
    '''
    Return the properly normalized function to perform the integration of the
    selected IMF. Returns mass values.
    '''
    imf_val = m_star*imfs(imf_sel, m_star, norm_const)
    return imf_val
    
    
def integral_IMF_N(m_star, imf_sel, norm_const):
    '''
    Return the properly normalized function to perform the integration of the
    selected IMF. Returns number of stars.
    '''
    imf_val = imfs(imf_sel, m_star, norm_const)
    return imf_val
    
    
def get_mass_cdf(cdf_arr):
    '''
    Draw a random mass value from the CDF.
    '''
    # Draw random [0,1) value.
    u = np.random.random()
    # Find closest value in CDF array.
    idx = np.argmin(np.abs(u - np.asarray(cdf_arr[1])))
    # Identify mass value corresponding to this CDF value.
    m_upper = cdf_arr[0][idx]

    return m_upper


def mass_dist(imf_sel, limit_sel, MN_total):
    '''
    Main function.
    
    Returns a mass distribution according to a given IMF and
    a total number of stars or total cluster mass.
    '''
    
    # Set IMFs limits.
    if imf_sel == 'kroupa_1993':
        m_low, m_high = 0.081, 100.
    elif imf_sel == 'chabrier_2001':
        m_low, m_high = 0.001, 100.
    
    # Normalize to a total unit mass.
    M = 1.
    norm_const = 1./quad(integral_IMF_M, m_low, m_high, args=(imf_sel, 1./M))[0]
    
    # Generate CDF for the given normalized IMF.
    cdf_arr = [[], []]
    for m_upper in np.arange(m_low, m_high, 0.01):
        cdf_arr[0].append(m_upper)
        cdf_arr[1].append(quad(integral_IMF_M, m_low, m_upper,\
        args=(imf_sel, norm_const))[0])
    
    # Prepare array for output masses.
    dist_mass = []
    # Fill in array according to selected limits.
    if limit_sel == 'total_number':
        while (len(dist_mass) < MN_total):
            # Draw random mass from the IMF's CDF.
            m_upper = get_mass_cdf(cdf_arr)
            # Store mass in array.
            dist_mass.append(round(m_upper,3))
        
    elif limit_sel == 'total_mass':
        while (sum(dist_mass) < MN_total):
            # Draw random mass from the IMF's CDF.
            m_upper = get_mass_cdf(cdf_arr)
            # Store mass in array.
            dist_mass.append(round(m_upper,3))

    return dist_mass
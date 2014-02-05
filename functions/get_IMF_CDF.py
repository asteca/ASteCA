# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 17:11:11 2014

@author: gabriel
"""

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


def IMF_CDF(imf_sel):
    '''
    Main function.
    
    Returns the selected IMF's cumulative distrubution function (CDF).
    '''
    
    # Set IMFs limits.
    if imf_sel == 'kroupa_1993':
        m_low, m_high, m_step = 0.081, 100., 0.01
    elif imf_sel == 'chabrier_2001':
        m_low, m_high, m_step = 0.001, 100., 0.01
    
    # Normalize to a total unit mass.
    M = 1.
    norm_const = 1./quad(integral_IMF_M, m_low, m_high, args=(imf_sel, 1./M))[0]
    
    # Generate CDF for the given normalized IMF.
    IMF_vals = []
    for m_upper in np.arange(m_low, m_high, m_step):
        IMF_vals.append(integral_IMF_M(m_upper, imf_sel, norm_const))
    # First sublist contains the masses, second the CDF values.
    cdf_arr = [np.arange(m_low, m_high, m_step), np.cumsum(IMF_vals)/m_high]


    return cdf_arr
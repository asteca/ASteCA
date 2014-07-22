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
        # 1993MNRAS.262..545K
        power = [-1.3, -2.2, -2.7]
        factor = [0.035, 0.019, 0.019]
        if 0.08 < m_star <= 0.5:
            i = 0
        elif 0.5 < m_star <= 1.:
            i = 1
        elif 1. < m_star:
            i = 2
        imf_val = norm_const * factor[i] * m_star ** power[i]

    elif imf_name == 'kroupa_2002':
        # Kroupa (2002) piecewise IMF (taken from MASSCLEAN article).
        power = [-0.3, -1.3, -2.3]
        m0, m1, m2 = 0.01, 0.08, 0.5
        factor = [0.4687, 0.0375, 0.01875]
        if 0.01 < m_star <= 0.08:
            i = 0
        elif 0.08 < m_star <= 0.5:
            i = 1
        elif 0.5 < m_star:
            i = 2
        imf_val = norm_const * factor[i] * m_star ** power[i]

    elif imf_name == 'chabrier_2001':
        # Chabrier (2001) exponential form of the IMF.
        imf_val = norm_const * 3. * m_star ** (-3.3) * \
            np.exp(-(716.4 / m_star) ** 0.25)

    return imf_val


def integral_IMF_M(m_star, imf_sel, tot_mass):
    '''
    Return the properly normalized function to perform the integration of the
    selected IMF. Returns mass values.
    '''
    imf_val = m_star * imfs(imf_sel, m_star, tot_mass)
    return imf_val


#def integral_IMF_N(m_star, imf_sel, norm_const):
#    '''
#    Return the properly normalized function to perform the integration of the
#    selected IMF. Returns number of stars.
#    '''
#    imf_val = imfs(imf_sel, m_star, norm_const)
#    return imf_val


def IMF_PDF(imf_sel):
    '''
    Main function.

    Returns the selected IMF's probability distribution function (PDF)
    normalized to 1 solar mass.
    '''

    # Set IMFs limits.
    if imf_sel == 'kroupa_1993':
        m_low = 0.081
    elif imf_sel == 'kroupa_2002':
        m_low = 0.011
    elif imf_sel == 'chabrier_2001':
        m_low = 0.001

    m_high, m_step = 100., 0.01
    # Normalize IMF to a total unit mass.
    M = 1.
    norm_const = 1. / quad(integral_IMF_M, m_low, m_high,
        args=(imf_sel, 1. / M))[0]

    # Generate PDF for the given normalized IMF. First sublist contains
    # the masses, second the PDF values.
    pdf_arr = [[], []]
    pdf_sum = 0.
    m_upper = m_low
    while pdf_sum < M:
        pdf_val = integral_IMF_M(m_upper, imf_sel, norm_const) * m_step
        pdf_arr[0].append(m_upper)
        pdf_arr[1].append(pdf_val)
        pdf_sum += pdf_val
        m_upper = m_upper + m_step

    return pdf_arr
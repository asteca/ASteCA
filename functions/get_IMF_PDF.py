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
        # http://adsabs.harvard.edu/abs/1993MNRAS.262..545K
        # Eq. (13), p. 572 (28)
        alpha = [-1.3, -2.2, -2.7]
        m_i = [0.08, 0.5, 1.]
        m0, m1, m2 = m_i
        factor = [0.035, 0.019, 0.019]
        if m0 < m_star <= m1:
            i = 0
        elif m1 < m_star <= m2:
            i = 1
        elif m2 < m_star:
            i = 2
        imf_val = norm_const * factor[i] * (m_star ** alpha[i])

    elif imf_name == 'kroupa_2002':
        # Kroupa (2002) piecewise IMF (taken from MASSCLEAN article).
        # http://adsabs.harvard.edu/abs/2002Sci...295...82K
        # Eq. (2) & (3), p. 1725
        alpha = [-0.3, -1.3, -2.3]
        m0, m1, m2 = [0.01, 0.08, 0.5]
        factor = [(1. / m1) ** alpha[0], (1. / m1) ** alpha[1],
            ((m2 / m1) ** alpha[1]) * ((1. / m2) ** alpha[2])]
        if m0 < m_star <= m1:
            i = 0
        elif m1 < m_star <= m2:
            i = 1
        elif m2 < m_star:
            i = 2
        imf_val = norm_const * factor[i] * (m_star ** alpha[i])

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
    Returns the selected IMF's probability distribution function (PDF)
    normalized to 1 solar mass.
    '''

    # Low mass limits are defined for each IMF to avoid numerical
    # issues when integrating.
    imfs_dict = {'kroupa_1993': (0.081), 'chabrier_2001': (0.001),
        'kroupa_2002': (0.011)}

    if imf_sel not in imfs_dict:
        print ("  WARNING: Name of IMF ({}) is incorrect.\n"
        "  Defaulting to Chabrier (2001).".format(imf_sel))
        imf_sel = 'chabrier_2001'

    # Set IMF low mass limit.
    m_low = imfs_dict[imf_sel]
    # Set IMF max mass limit and interpolation step.
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

    # Normalize probabilities to avoid 'np.random.choice' error if PDF doesn't
    # add up *exactly* to 1.  with numpy > 1.8.0.
    pdf_arr[1] /= np.asarray(pdf_arr[1]).sum()

    return np.asarray(pdf_arr)
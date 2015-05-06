# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 17:11:11 2014

@author: gabriel
"""

import numpy as np
from scipy.integrate import quad
from .._in import get_in_params as g


def imfs(imf_name, m_star):
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
        imf_val = factor[i] * (m_star ** alpha[i])

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
        imf_val = factor[i] * (m_star ** alpha[i])

    elif imf_name == 'chabrier_2001_log':
        # Chabrier (2001) lognormal form of the IMF.
        # http://adsabs.harvard.edu/abs/2001ApJ...554.1274C
        # Eq (7)
        imf_val = (1. / (np.log(10) * m_star)) * 0.141 * \
            np.exp(-((np.log10(m_star) - np.log10(0.1)) ** 2) /
                (2 * 0.627 ** 2))

    elif imf_name == 'chabrier_2001_exp':
        # Chabrier (2001) exponential form of the IMF.
        # http://adsabs.harvard.edu/abs/2001ApJ...554.1274C
        # Eq (8)
        imf_val = 3. * m_star ** (-3.3) * np.exp(-(716.4 / m_star) ** 0.25)

    return imf_val


def integral_IMF_M(m_star, imf_sel):
    '''
    Return the properly normalized function to perform the integration of the
    selected IMF. Returns mass values.
    '''
    imf_val = m_star * imfs(imf_sel, m_star)
    return imf_val


def integral_IMF_N(m_star, imf_sel):
    '''
    Return the properly normalized function to perform the integration of the
    selected IMF. Returns number of stars.
    '''
    imf_val = imfs(imf_sel, m_star)
    return imf_val


def N_IMF():
    '''
    Returns the number of stars per interval of mass for the selected IMF.
    '''

    imf_sel, m_high = g.sc_params[0], g.sc_params[1]

    # Low mass limits are defined for each IMF to avoid numerical
    # issues when integrating.
    imfs_dict = {'chabrier_2001_exp': (0.01), 'chabrier_2001_log': (0.01),
        'kroupa_1993': (0.081), 'kroupa_2002': (0.011)}

    # Set IMF low mass limit.
    m_low = imfs_dict[imf_sel]
    # Set IMF max mass limit and interpolation step.
    # For m_high > 100 Mo the differences in the resulting normalization
    # constant are negligible. This is because th IMF drops very rapidly for
    # high masses.
    # The step (m_step) should not be too small since it will have an impact
    # on the performance of the get_mass_dist function.
    m_step = 0.1

    # Obtain normalization constant. This is equivalent to 'k' in Eq. (7)
    # of Popescu & Hanson 2009 (138:1724-1740; PH09)
    norm_const = 1. / quad(integral_IMF_M, m_low, m_high, args=(imf_sel))[0]

    # Obtain number of stars in each mass interval. Equivalent to the upper
    # fraction of Eq. (8) in PH09, without the total mass.
    st_dist = [[], []]
    m_upper = m_low
    while m_upper < m_high:
        m_lower = m_upper
        m_upper = m_upper + m_step
        # Number of stars in the (m_lower, m_upper) interval.
        N_stars = quad(integral_IMF_N, m_lower, m_upper, args=(imf_sel))[0]
        st_dist[0].append(m_upper)
        st_dist[1].append(N_stars)

    # Normalize number of stars by constant.
    st_dist[1] = np.asarray(st_dist[1]) * norm_const

    return st_dist

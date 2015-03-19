"""
@author: gabriel
"""

import numpy as np


def exp_3p(x, a, b, c):
    '''
    Three-params exponential function.
    '''
    return a * np.exp(b * np.asarray(x)) + c


def exp_2p(x, a, b):
    '''
    Two-params exponential function.
    '''
    return a * np.exp(b * np.asarray(x))

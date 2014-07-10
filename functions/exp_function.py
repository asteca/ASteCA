"""
@author: gabriel
"""

import numpy as np


# Define exponential function.
def exp_func(x, a, b, c):
    '''
    Exponential function.
    '''
    return a * np.exp(b * x) + c

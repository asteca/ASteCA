
import numpy as np


def exp_3p(x, a, b, c):
    '''
    Three-params exponential function.
    '''
    return a * np.exp(b * np.asarray(x)) + c


def exp_2p(x, a, b):
    '''
    Two-params exponential function.
    This function is tied to the 'lowexp' error rejection function and the
    'synt_cl_err' function.
    '''
    return a * np.exp(b * np.asarray(x))

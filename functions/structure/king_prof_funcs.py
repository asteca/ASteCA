# -*- coding: utf-8 -*-
"""
Created on Sat Feb 28 12:00:00 2015

@author: gabriel
"""

import numpy as np


def two_params(x, cd, rc, fd):
    '''
    Two parameters King profile fit.
    '''
    rc = 0.0001 if rc < 0.0001 else rc
    return fd + cd / (1 + (np.asarray(x) / rc) ** 2)


def three_params(x, rt, cd, rc, fd):
    '''
    Three parameters King profile fit.
    '''
    rc = 0.0001 if rc < 0.0001 else rc
    return cd * (1 / np.sqrt(1 + (np.asarray(x) / rc) ** 2) -
                 1 / np.sqrt(1 + (rt / rc) ** 2)) ** 2 + fd

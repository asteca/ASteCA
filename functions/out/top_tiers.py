# -*- coding: utf-8 -*-
"""
Created on Wed Dic 10 12:00:00 2014

@author: gabriel
"""

import numpy as np
from .._in import get_in_params as g


def get_top_tiers(bf_return):
    '''
    Obtain top tier models, produce data file and output image.
    '''

    isoch_fit_params, isoch_fit_errors, shift_isoch, synth_clst = bf_return

    return
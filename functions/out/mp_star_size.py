# -*- coding: utf-8 -*-
"""
Created on Tue Dic 16 12:00:00 2014

@author: gabriel
"""

import numpy as np


def star_size(mag_data):
    '''
    Convert magnitudes into intensities and define sizes of stars in
    finding chart.
    '''
    # Scale factor.
    factor = 500. * (1 - 1 / (1 + 150 / len(mag_data) ** 0.85))
    return 0.1 + factor * 10 ** ((np.array(mag_data) - min(mag_data)) / -2.5)
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 14:11:22 2013

@author: gabriel
"""

'''
Calculate the QQ-plot for the distribution of p-values obtained comparing
the cluster's KDE with the field region's KDEs.
'''

import numpy as np
from scipy.stats.mstats import mquantiles


def ppoints(vector):
    '''
    Analogue to R's `ppoints` function
    see details at 'http://stat.ethz.ch/R-manual/R-patched/library/stats/html/
    ppoints.html'
    '''
    try:
        n = np.float(len(vector))
    except TypeError:
        n = np.float(vector)
    a = 3. / 8. if n <= 10 else 1. / 2

    return (np.arange(n) + 1 - a) / (n + 1 - 2 * a)


def qqplot(p_vals_cl, p_vals_f):

    # Interpolate the larger list.
    if len(p_vals_f) >= len(p_vals_cl):
        A, B = p_vals_f, p_vals_cl
    else:
        B, A = p_vals_f, p_vals_cl
    # Calculate the quantiles, using R's defaults for 'alphap' and 'betap'
    # (ie: R's type 7)
    # See: 'http://docs.scipy.org/doc/scipy/reference/generated/
    # scipy.stats.mstats.mquantiles.html'
    quant = mquantiles(A, prob=ppoints(B), alphap=1., betap=1.)

    # Set order so the names of the axis when plotting are unchanged.
    if len(p_vals_f) >= len(p_vals_cl):
        quantiles = [sorted(B), sorted(quant.tolist())]
    else:
        quantiles = [sorted(quant.tolist()), sorted(B)]

    # Calculate CCC (Concordance correlation coefficient) for the quantiles.
    # https://en.wikipedia.org/wiki/Concordance_correlation_coefficient
    a = quantiles[0]
    b = quantiles[1]
    ccc = 2 * (np.std(a) * np.std(b)) / \
    (np.std(a) ** 2 + np.std(b) ** 2 + (np.mean(a) - np.mean(b)) ** 2)
    # Invert CCC for easier reading.
    ccc = 1 - ccc

    # Return parameters inside list.
    qq_params = [ccc, quantiles]

    return qq_params
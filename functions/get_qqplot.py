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
    see details at http://stat.ethz.ch/R-manual/R-patched/library/stats/html/ppoints.html 
    '''
    try:
        n = np.float(len(vector))
    except TypeError:
        n = np.float(vector)
    a = 3./8. if n <= 10 else 1./2
    
    return (np.arange(n) + 1 - a)/(n + 1 - 2*a)
    
    
def qqplot(p_vals_cl, p_vals_f):
    
    # Calculate the quantiles, using R's defaults for 'alphap' and 'betap'
    # (ie: R's type 7)    
    quant = mquantiles(p_vals_f, prob=ppoints(p_vals_cl), alphap=1., betap=1.)
    
    print quant

    raw_input()
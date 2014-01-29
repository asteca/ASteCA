# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 10:37:16 2014

@author: gabriel
"""

import numpy as np

def mag_completeness(mag_data, b_width):
    '''
    Calculate the completeness level in each magnitude bin beyond the one
    with the maximum count (ie: the assumed 100% completeness limit)
    '''
    # Max value of magnitude.
    max_mag = max(mag_data)    
    # Number of bins.
    bins = int((max(mag_data)-min(mag_data))/b_width)
    # Get histogram.
    mag_hist, bin_edges = np.histogram(mag_data, bins)
    # Index of the bin with the maximum number of stars.
    max_indx = mag_hist.argmax(axis=0)
    # Value of the magnitude where the peak starts.
#    mag_peak = bin_edges[max_indx]
    # Total number of stars beyond the peak bin (included).
    total = sum(mag_hist[max_indx:])
    # Get percentages per interval beyond the maximum interval (included).
    # These values tell me the percentages of stars beyond the magnitude peak
    # that are located inside each magnitude bin. The peak magnitude bin (the
    # first one) will have the biggest percentage.
    comp_perc = [(i*100.)/total for i in mag_hist[max_indx:]]
    
    # Store everything in a single list.
    completeness = [max_mag, bin_edges, max_indx, comp_perc]
    
    return completeness

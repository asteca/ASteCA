# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 15:51:34 2013

@author: gabriel
"""


def cont_indx(field_dens, a_clust, n_clust):
    '''
    Calculate the contamination index value. This parameter is defined as the
    ratio of field stars density over the density of stars in the cluster
    region.

    A small number (close to zero) means the field contamination in the
    cluster region is very small.
    If this number equals 1, it means that an equal number of field stars
    and cluster members are expected inside the cluster region. A number
    larger than 1 means there are more field stars on average than there
    are cluster members inside the cluster region (which isn't a good
    sign).
    '''

    # Star density in the cluster region.
    cl_dens = n_clust / a_clust

    # Final contamination index.
    cont_index = field_dens / cl_dens

    return cont_index
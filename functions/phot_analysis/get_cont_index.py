# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 15:51:34 2013

@author: gabriel
"""


def cont_indx(n_clust, cl_area, field_dens):
    '''
    Calculate the contamination index value. This parameter is defined as the
    ratio of field stars density over the density of stars in the cluster
    region.

    A small number (close to zero) means the field contamination in the
    cluster region is very small.
    If this number equals 0.5, it means that an equal number of field stars
    and cluster members are expected inside the cluster region. A value of
    1 means there are no expected cluster members inside the cluster region
    (which isn't a good sign).
    '''

    # Star density in the cluster region.
    cl_dens = n_clust / cl_area

    # Final contamination index.
    cont_index = field_dens / cl_dens

    if cont_index >= 1.:
        print "  WARNING: contamination index obtained is high: {:.2f}".format(
            cont_index)
    else:
        print 'Contamination index obtained ({:.2f}).'.format(cont_index)

    return cont_index
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 15:51:34 2013

@author: gabriel
"""

'''
Calculate and return the contamination index value cont_index.
The contamination index is defined as the ratio of field stars that should 
be present in the cluster region (calculated by means of the background value
'backg_val' and the cluster's radious 'r') to the totalnumber of stars in the
cluster region 'n_tot'.

CI = [backg_val*PI*r**2]/n_tot

If this number equals 1, it means that all of the stars in the cluster region
are expected to be field stars. A small number (close to zero) means the
field contamination in the cluster region is very small. A number larger
than 1 means there are more stars in average in the background than there are
inside the cluster region (which isn't a good sign).
'''

import numpy as np

def cont_indx(backg_val, clust_rad, stars_in, stars_in_rjct):
    
    # Cluster's area.
    a_c = np.pi*(clust_rad**2)
    
    # Total number of stars in cluster region.
    n_tot = len(stars_in) + len(stars_in_rjct)
    
    # Calculate contamination index.
    cont_index = backg_val*a_c/n_tot
    
    return cont_index
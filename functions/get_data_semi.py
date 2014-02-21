# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 11:12:55 2013

@author: gabriel
"""

from os.path import join

def get_semi(mypath, clust_name):
    '''
    Get center, radius and flags for semi automatic mode.
    '''
    # Flag to indicate if cluster was found in file.
    flag_clust_found = False
    myfile = 'clusters_input.dat'
    with open(join(mypath, myfile), mode="r") as f_cl_dt:
        for line in f_cl_dt:
            li=line.strip()
            # Skip comments.
            if not li.startswith("#"):
                reader = li.split()            
            
                # If cluster is found in file.
                if reader[0] == clust_name:
                    cl_cent_semi = [float(reader[1]), float(reader[2])]
                    cl_rad_semi = float(reader[3])
                    cent_flag_semi, rad_flag_semi, err_flag_semi = \
                    int(reader[4]), int(reader[5]), int(reader[6])
                    # Set flag to True if the cluster was found.                    
                    flag_clust_found = True
                    
    # If cluster was found.
    if flag_clust_found:
        semi_return = [cl_cent_semi, cl_rad_semi, cent_flag_semi, rad_flag_semi,
                       err_flag_semi]
    else:
        semi_return = []

    return semi_return
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 11:12:55 2013

@author: gabriel
"""

from os.path import join

def get_semi(mypath, clust_name):
    '''
    Get center, radius and errors flag for semi automatic mode.
    '''
    
    cent_cl_semi = [[],[]]
    # Get center, radius and error flag from file.
    myfile = 'data_input'
    with open(join(mypath, myfile), mode="r") as f_cl_dt:
       
        for line in f_cl_dt:

            li=line.strip()
            
            # Jump comments.
            if not li.startswith("#"):
                reader = li.split()            
            
                # reader[1]=c_x, reader[2]=c_y, [3]=radius, [8], [9], [10]=
                # center, radius and error flags.
                if reader[0] == clust_name:
                    cent_cl_semi[0], cent_cl_semi[1], cl_rad_semi, \
                    cent_flag_semi, rad_flag_semi, err_flag_semi = \
                    float(reader[1]), float(reader[2]), float(reader[3]), \
                    int(reader[8]), int(reader[9]), int(reader[10])

    return cent_cl_semi, cl_rad_semi, cent_flag_semi, rad_flag_semi, \
    err_flag_semi   
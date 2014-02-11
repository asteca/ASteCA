# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 14:03:44 2014

@author: gabriel
"""

import csv
from os.path import join


def get_in_params(mypath):
    '''
    This function reads the input data parameters stored in the 'ocaat_input.dat'
    file and returns them packaged for each function to use.
    '''
    
# Allows to work with columns data files.

    data_file = join(mypath, 'ocaat_input.dat')
    
    with open(data_file, mode="r") as f_dat:
        
        # Iterate through each line in the file.
        for line in f_dat:
            
            if not line.startswith("#") and line.strip() != '':
                reader = line.split()
                
                # Read folder paths where clusters are stored.
                if reader[0] == 'MO':
                    mode = str(reader[1])
                elif reader[0] == 'CP0':
                    mypath2 = str(reader[1])
                elif reader[0] == 'CP1':
                    mypath3 = str(reader[1])
                elif reader[0] == 'CP2':
                    output_dir = str(reader[1])
                elif reader[0] == 'PD':
                    gd_params = map(int, reader[1:])
                elif reader[0] == 'CC':
                    gc_params = map(float, reader[1:])
                    
                    
    in_dirs = [mypath2, mypath3, output_dir]
    
    return mode, in_dirs, gd_params, gc_params
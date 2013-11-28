"""
@author: gabriel
"""

import numpy as np

# Allows to work with columns data files.
import csv
from os.path import join

def get_data(mypath, sub_dir, myfile):
    '''
    Get data from the cluster's data files.
    '''

    data_file = join(mypath, sub_dir, myfile)

    # Obtain the number of columns in the data file
    with open(data_file) as f_cluster:
        # Read the f_cluster file using this call which separates each line in
        # columns automatically.
        reader = csv.reader(f_cluster, delimiter=' ', skipinitialspace=True)
        
        # Skip the first 20 lines of the file and then count the number of
        # columns. We do this to skip the header of the file if one exists.
        for _ in range(20):
            row_20 = next(reader)
        # Save number of columns in 'num_cols'.
        num_cols = len(row_20)
        
    # Headers
    #7 columnas = '# ID   x   y T1   sT1   C-T1   sCT1'
    #8 columnas = '# ID   x   y T1   sT1   C-T1   sCT1   No_obs'
    #9 columnas = '# ID   x   y T1   sT1   No_obs   C-T1   sCT1   No_obs'
    #10, 11, 12 columnas = etc 9
    
    # Loads the data in 'myfile' as a list of N lists where N is the number of
    # columns. Each of the N lists contains all the data for the column.
    data = np.loadtxt(data_file, unpack=True)
                    
    if num_cols <= 9:
        id_star, x_data, y_data, mag_data, e_mag, col1_data, e_col1 = data[0], \
        data[1], data[2], data[3], data[4], data[5], data[6]
    elif num_cols == 10:
        # MASSCLEAN input file
        id_star, x_data, y_data, mag_data, e_mag, col1_data, e_col1 = data[0], \
        data[8], data[9], data[2], data[0]/500., data[1]-data[2], data[0]/500.    
    else:
        id_star, x_data, y_data, mag_data, e_mag, col1_data, e_col1 = data[0], \
        data[1], data[2], data[3], data[4], data[6], data[7]
        
    return id_star, x_data, y_data, mag_data, e_mag, col1_data, e_col1
"""
@author: gabriel
"""

import numpy as np

from os.path import join

def get_data(mypath, sub_dir, myfile, gd_params):
    '''
    Get photometric data from the cluster's data file.
    '''

    data_file = join(mypath, sub_dir, myfile)

    # Loads the data in 'myfile' as a list of N lists where N is the number of
    # columns. Each of the N lists contains all the data for the column.
    # If any string is found (for example 'INDEF') it is converted to 99.999.
    data = np.genfromtxt(data_file, dtype=float, filling_values=99.999,
                         unpack=True)
                    
    # Read indexes from input params.
    id_inx, x_inx, y_inx, m_inx, em_inx, c_inx, ec_inx = gd_params
    
    # Read data colums.
    id_star, x_data, y_data, mag_data, e_mag, col1_data, e_col1 = \
    data[id_inx], data[x_inx], data[y_inx], data[m_inx], data[em_inx],\
    data[c_inx], data[ec_inx]

    return id_star, x_data, y_data, mag_data, e_mag, col1_data, e_col1
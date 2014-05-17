"""
@author: gabriel
"""

import numpy as np
from os.path import join


def rem_bad_stars(id_star, x_data, y_data, mag_data, e_mag, col1_data,
        e_col1):
    '''
    Remove stars from all lists that have too large magnitude or color
    values (or their errors) which indicates a bad photometry.
    '''
    # Store indexes of stars that should be removed.
    max_lim = 90.
    del_indexes = []
    for i, data_lsts in enumerate(zip(mag_data, e_mag, col1_data, e_col1)):
        if any(e > max_lim for e in data_lsts):
            del_indexes.append(i)

    # Remove stars from all lists simultaneously.
    big_array = np.array([id_star, x_data, y_data, mag_data, e_mag, col1_data,
        e_col1])
    clean_array = np.delete(big_array, del_indexes, axis=1)

    return clean_array


def get_data(mypath, sub_dir, myfile, gd_params):
    '''
    Get photometric data from the cluster's data file.
    '''
    data_file = join(mypath, sub_dir, myfile)

    # Read indexes from input params.
    id_inx, x_inx, y_inx, m_inx, em_inx, c_inx, ec_inx = gd_params

    # Loads the data in 'myfile' as a list of N lists where N is the number of
    # columns. Each of the N lists contains all the data for the column.
    # If any string is found (for example 'INDEF') it is converted to 99.999.
    data = np.genfromtxt(data_file, dtype=float, filling_values=99.999,
                         unpack=True)

    # Read data colums, except IDs.
    x_data, y_data, mag_data, e_mag, col1_data, e_col1 = \
    data[x_inx], data[y_inx], data[m_inx], data[em_inx],\
    data[c_inx], data[ec_inx]

    # Now read IDs as strings.
    data = np.genfromtxt(data_file, dtype=str, unpack=True)
    id_star = data[id_inx]

    # If any mag or color value (or their errors) is too large, discard
    # that star.
    print len(id_star)
    id_star, x_data, y_data, mag_data, e_mag, col1_data, e_col1 = \
    rem_bad_stars(id_star, x_data, y_data, mag_data, e_mag, col1_data, e_col1)

    print 'Data obtained from input file (N stars: %d).' % len(id_star)

    return id_star, x_data, y_data, mag_data, e_mag, col1_data, e_col1
"""
@author: gabriel
"""

import numpy as np
import get_in_params as g


def rem_bad_stars(id_star, x_data, y_data, mag_data, e_mag, col_data,
        e_col):
    '''
    Remove stars from all lists that have too large magnitude or color
    values (or their errors) which indicates a bad photometry.
    '''
    # Set photometric range for accepted stars.
    min_lim, max_lim = -50., 50.

    # Store indexes of stars that should be removed.

    # For magnitudes.
    m_del_ind = []
    for i, mag in enumerate(mag_data):
        lists_arr = zip(mag, e_mag[i])
        m_del_ind.append([i for i, t in enumerate(lists_arr) if
            any(e > max_lim for e in t) or any(e < min_lim for e in t)])
    # For colors.
    c_del_ind = []
    for i, col in enumerate(col_data):
        lists_arr = zip(col, e_col[i])
        c_del_ind.append([i for i, t in enumerate(lists_arr) if
            any(e > max_lim for e in t) or any(e < min_lim for e in t)])

    # Combine magnitude and color indexes of lines to be removed.
    del_ind = m_del_ind + c_del_ind
    del_indexes = [_ for sublist in del_ind for _ in sublist]

    # Remove stars from id list.
    id_clean = np.delete(np.array(id_star), del_indexes)

    # Remove stars from the x, y coordinates.
    coord_clean_arr = np.delete(np.array([x_data, y_data]), del_indexes, axis=1)

    # Remove stars from the magnitude columns.
    m_clean_arr, em_clean_arr = [], []
    for i, mag in enumerate(mag_data):
        m_clean_arr.append(np.delete(np.array(mag), del_indexes))
        em_clean_arr.append(np.delete(np.array(e_mag[i]), del_indexes))

    # Remove stars from the color columns.
    c_clean_arr, ec_clean_arr = [], []
    for i, col in enumerate(col_data):
        c_clean_arr.append(np.delete(np.array(col), del_indexes))
        ec_clean_arr.append(np.delete(np.array(e_col[i]), del_indexes))

    return id_clean, coord_clean_arr, m_clean_arr, em_clean_arr, c_clean_arr,\
    ec_clean_arr


def get_data(data_file, phot_params):
    '''
    Get spatial and photometric data from the cluster's data file.
    '''

    # Read indexes from input params.
    id_inx, x_inx, y_inx = g.gd_params[0][:-1]
    # INdexes.
    magnitudes, colors, e_mags, e_cols = phot_params[0]

    # Loads the data in 'myfile' as a list of N lists where N is the number of
    # columns. Each of the N lists contains all the data for the column.
    # If any string is found (for example 'INDEF') it is converted to 99.999.
    data = np.genfromtxt(data_file, dtype=float, filling_values=99.999,
                         unpack=True)

    # Extract coordinates and photometric data colums, except IDs.
    x_data, y_data = data[x_inx], data[y_inx]
    mag_data = [data[i] for i in magnitudes]
    e_mag = [data[i] for i in e_mags]
    col_data = [data[i] for i in colors]
    e_col = [data[i] for i in e_cols]

    # Now read IDs as strings.
    data = np.genfromtxt(data_file, dtype=str, unpack=True)
    id_star = data[id_inx]
    n_old = len(id_star)

    # If any mag or color value (or their errors) is too large, discard
    # that star.
    id_star, [x_data, y_data], mag_data, e_mag, col_data, e_col = \
    rem_bad_stars(id_star, x_data, y_data, mag_data, e_mag, col_data, e_col)

    print 'Data obtained from input file (N_stars: %d).' % len(id_star)
    if (n_old - len(id_star)) > 0:
        print ' Entries rejected: %d.' % (n_old - len(id_star))

    return [id_star, x_data, y_data], [mag_data, e_mag, col_data, e_col]
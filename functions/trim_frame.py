# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 10:20:44 2013

@author: gabriel
"""

import numpy as np
import matplotlib.pyplot as plt
from display_frame import disp_frame as d_f
import get_in_params as g


def trim_frame(id_coords, phot_data):
    '''
    Trim frame according to given values of new center and side lengths.
    '''

    if g.mode == 'manual':

        # Unpack data.
        id_star, x_data, y_data = id_coords
        mag_data, e_mag, col_data, e_col = phot_data

        # Show full frame plot.
        d_f(x_data, y_data, mag_data)
        plt.show()

        wrong_answer = True
        while wrong_answer:

            temp_cent, temp_side = [], []
            answer_fra = raw_input('Trim frame? (y/n) ')

            if answer_fra == 'n':
                id_coords_t, phot_data_t = id_coords, phot_data
                wrong_answer = False

            elif answer_fra == 'y':
                print 'Input center of new frame (in px).'
                temp_cent.append(float(raw_input('x: ')))
                temp_cent.append(float(raw_input('y: ')))
                print 'Input side lenght for new frame (in px).'
                temp_side.append(float(raw_input('x_side: ')))
                temp_side.append(float(raw_input('y_side: ')))
                wrong_answer = False

                # Indexes of elements to remove.
                del_indexes = []
                # Iterate through all stars.
                for st_indx, star in enumerate(id_star):
                    # Check if star is outside new frame boudaries.
                    if abs(temp_cent[0] - x_data[st_indx]) > temp_side[0] / 2.\
                    or abs(temp_cent[1] - y_data[st_indx]) > temp_side[1] / 2.:
                        del_indexes.append(st_indx)

                # Remove stars from id list.
                id_clean = np.delete(np.array(id_star), del_indexes)
                # Remove stars from the x, y coordinates.
                coord_clean_arr = np.delete(np.array([x_data, y_data]),
                    del_indexes, axis=1)
                # Remove stars from the magnitude columns.
                m_clean_arr, em_clean_arr = [], []
                for i, mag in enumerate(mag_data):
                    m_clean_arr.append(np.delete(np.array(mag), del_indexes))
                    em_clean_arr.append(np.delete(np.array(e_mag[i]),
                    del_indexes))
                # Remove stars from the color columns.
                c_clean_arr, ec_clean_arr = [], []
                for i, col in enumerate(col_data):
                    c_clean_arr.append(np.delete(np.array(col), del_indexes))
                    ec_clean_arr.append(np.delete(np.array(e_col[i]),
                    del_indexes))

                id_coords_t = [id_clean, coord_clean_arr]
                phot_data_t = [m_clean_arr, em_clean_arr, c_clean_arr,
                    ec_clean_arr]
            else:
                print 'Wrong input. Try again.\n'
    else:
        id_coords_t, phot_data_t = id_coords, phot_data

    return id_coords_t, phot_data_t

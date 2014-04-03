# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 10:20:44 2013

@author: gabriel
"""

import numpy as np
import matplotlib.pyplot as plt
from display_frame import disp_frame as d_f


def trim_frame(phot_data, mode):
    '''
    Trim frame according to given values of new center and side lengths.
    '''

    if mode == 'm':

        # Unpack data.
        id_star, x_data, y_data, mag_data, e_mag, col1_data, e_col1 = phot_data

        # Show plot with center obtained.
        d_f(x_data, y_data, mag_data)
        plt.show()

        wrong_answer = True
        while wrong_answer:

            temp_cent, temp_side = [], []
            answer_fra = raw_input('Trim frame? (y/n) ')

            if answer_fra == 'n':
                phot_data_t = phot_data
                wrong_answer = False

            elif answer_fra == 'y':
                print 'Input center of new frame (in px).'
                temp_cent.append(float(raw_input('x: ')))
                temp_cent.append(float(raw_input('y: ')))
                print 'Input side lenght for new frame (in px).'
                temp_side.append(float(raw_input('x_side: ')))
                temp_side.append(float(raw_input('y_side: ')))
                wrong_answer = False

                # Empty new lists.
                id_star2, x_data2, y_data2, T1_data2, e_T12, CT1_data2, \
                e_CT12 = [], [], [], [], [], [], []

                # Iterate through all stars.
                for st_indx, star in enumerate(id_star):

                    # Check if star is inside new frame boudaries.
                    if abs(temp_cent[0] - x_data[st_indx]) < temp_side[0] / 2.\
                    and abs(temp_cent[1] - y_data[st_indx]) < temp_side[1] / 2.:

                        id_star2.append(star)
                        x_data2.append(x_data[st_indx])
                        y_data2.append(y_data[st_indx])
                        T1_data2.append(mag_data[st_indx])
                        e_T12.append(e_mag[st_indx])
                        CT1_data2.append(col1_data[st_indx])
                        e_CT12.append(e_col1[st_indx])

                phot_data_t = [np.array(id_star2), np.array(x_data2),
                    np.array(y_data2), np.array(T1_data2), np.array(e_T12),
                    np.array(CT1_data2), np.array(e_CT12)]
            else:
                print 'Wrong input. Try again.\n'

    return phot_data_t

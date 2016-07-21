
import numpy as np
import matplotlib.pyplot as plt
import display_frame
from ..inp import input_params as g


def main(phot_data):
    '''
    Trim frame according to given values of new center and side lengths.
    '''

    if g.mode == 'manual':

        # Unpack data.
        id_star, x_data, y_data, mag_data, e_mag, col1_data, e_col1 = phot_data

        # Show full frame plot.
        display_frame.main(x_data, y_data, mag_data)
        plt.show()

        # Ask to trim frame.
        while True:
            temp_cent, temp_side = [], []
            answer_fra = raw_input('Trim frame? (y/n) ')
            if answer_fra == 'n':
                phot_data_t = phot_data
                break
            elif answer_fra == 'y':
                try:
                    print('Input center of new frame (in px).')
                    temp_cent.append(float(raw_input('x: ')))
                    temp_cent.append(float(raw_input('y: ')))
                    print('Input full side length for new frame (in px).')
                    temp_side.append(float(raw_input('x_side: ')))
                    temp_side.append(float(raw_input('y_side: ')))
                    # Empty new lists.
                    id_star2, x_data2, y_data2, mag_data2, e_m2, col_data2, \
                        e_c2 = [], [], [], [], [], [], []

                    # Iterate through all stars.
                    for st_indx, star in enumerate(id_star):

                        # Check if star is inside new frame boundaries.
                        if abs(temp_cent[0] - x_data[st_indx]) < \
                            temp_side[0] / 2. and \
                            abs(temp_cent[1] - y_data[st_indx]) < \
                                temp_side[1] / 2.:

                            id_star2.append(star)
                            x_data2.append(x_data[st_indx])
                            y_data2.append(y_data[st_indx])
                            mag_data2.append(mag_data[st_indx])
                            e_m2.append(e_mag[st_indx])
                            col_data2.append(col1_data[st_indx])
                            e_c2.append(e_col1[st_indx])

                    phot_data_t = [np.array(id_star2), np.array(x_data2),
                                   np.array(y_data2), np.array(mag_data2),
                                   np.array(e_m2), np.array(col_data2),
                                   np.array(e_c2)]
                    break
                except:
                    print("Sorry, input is not valid. Try again.")
            else:
                print("Sorry, input is not valid. Try again.\n")
    else:
        phot_data_t = phot_data

    return phot_data_t

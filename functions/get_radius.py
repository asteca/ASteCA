"""
@author: gabriel
"""

import numpy as np
from functions.display_rad import disp_rad as d_r
import matplotlib.pyplot as plt


def main_center_algor(rdp_params, cr_params, backg_value, bin_width):
    '''
    This function holds the main algorithm that returns a radius value.
    '''

    radii, ring_density = rdp_params[:2]
    # Find maximum density value and assume this is the central density.
    # Do not use previous values.
    max_dens_ind = np.argmax(ring_density)
    ring_dens_c, radii_c = ring_density[max_dens_ind:], radii[max_dens_ind:]

    # Assign a value to the number of points that should be found below
    # the delta values around the background to attain the 'stabilized'
    # condition.
    mode_r = cr_params[0]
    if mode_r not in {'auto', 'manual'}:
        print "  WARNING: CR mode is not valid. Default to 'auto'."
        mode_r = 'auto'
    # Set params.
    if mode_r == 'manual':
        n_left = int(cr_params[1])
    elif mode_r == 'auto':
        n_left = max(int(round(len(radii_c) * 0.2)), 3)

    # Delta step is fixed to 5%
    delta_step = 5

    # Difference between max density value and the background value.
    delta_total = (max(ring_dens_c) - backg_value)

    # If the difference between the max density value and the background is
    # less than 3 times the value of the background, raise a flag.
    flag_delta_total = False
    if delta_total < 3 * backg_value:
        flag_delta_total = True

    # Initialize index_rad value.
    index_rad = []
    for i in range(4):

        # Store value for delta_percentage --> 20, 15, 10, 5
        delta_percentage = (4. - i) * delta_step

        # % of difference between max density value and background.
        delta_backg = delta_percentage * delta_total / 100.

        # Initialize density values counter for points that fall inside the
        # range determined by the delta value around the background.
        in_delta_val, dens_dist, index_rad_i = 0, 1.e10, 0

        # Iterate through all values of star density in this "square ring".
        for index, item in enumerate(ring_dens_c):

            # Condition to iterate until at least n_left points below the
            # delta + background value are found.
            if in_delta_val < (n_left - i):

                # If the density value is closer than 'delta_backg' to the
                # background value or lower --> add it.
                if item <= delta_backg + backg_value:
                    # Augment value of counter.
                    in_delta_val += 1
                    # Store first radius value that falls below the upper delta
                    # limit.
                    if abs(item - backg_value) < dens_dist:
                        # Store distance of point to field density value.
                        dens_dist = abs(item - backg_value)
                        index_rad_i = index
                else:
                    # Reset.
                    in_delta_val, index_rad_i, dens_dist = 0, 0, 1.e10
            else:
                index_rad.append(index_rad_i)
                break

    # Raise a flag if only two radius values were found under the 4 deltas.
    flag_delta, flag_not_stable = False, False
    if len(index_rad) < 3:
        flag_delta = True

    rad_found = []
    for ind in index_rad:
        rad_found.append(radii_c[ind])
    if rad_found:
        clust_rad, e_rad = np.mean(rad_found), max(np.std(rad_found), bin_width)
    else:
        flag_not_stable = True
        # No radius value found. Assign radius value as the middle element
        # in the radii list.
        clust_rad, e_rad = radii_c[int(len(radii_c) / 2.)], 0.
        print '  WARNING: no radius found, setting value to: %0.2f' % clust_rad

    return clust_rad, e_rad, flag_delta_total, flag_not_stable, flag_delta


def get_clust_rad(phot_data, backg_value, cr_params, center_params,
    rdp_params, semi_return, mode, bin_width):
    """
    Obtain the value for the cluster's radius by counting the number of points
    that fall within a given interval of the background or lower. If this number
    is equal to a fixed number of points n_left then assign the radius as the
    closest point to the background value among the first n_left points
    counting from the first one that fell below the backg + delta limit.
    Iterate increasing the interval around the background until n_left points
    are found or the delta interval reaches its maximum allowed.
    """

    # Call function that holds the radius finding algorithm.
    clust_rad, e_rad, flag_delta_total, flag_not_stable, flag_delta = \
    main_center_algor(rdp_params, cr_params, backg_value, bin_width)

    # Check if semi or manual mode are set.
    flag_radius_manual = False
    if mode == 'auto':
        print 'Auto radius found: %0.1f px.' % clust_rad

    elif mode == 'semi':
        # Unpack semi values.
        cent_cl_semi, cl_rad_semi, cent_flag_semi, rad_flag_semi, \
        err_flag_semi = semi_return

        if rad_flag_semi == 1:
            # Update values.
            clust_rad, e_rad = cl_rad_semi, 0.
            print 'Semi radius set: %0.1f px.' % clust_rad
        else:
            print 'Auto radius found: %0.1f px.' % clust_rad

    # If Manual mode is set, display radius and ask the user to accept it or
    # input new one.
    elif mode == 'manual':

        print 'Radius found: %0.1f px' % clust_rad
        d_r(phot_data, center_params, clust_rad, backg_value, rdp_params)
        plt.show()

        wrong_answer = True
        while wrong_answer:
            answer_rad = raw_input('Accept radius (otherwise input new one \
manually)? (y/n) ')
            if answer_rad == 'y':
                print 'Value accepted.'
                wrong_answer = False
            elif answer_rad == 'n':
                clust_rad_m = float(raw_input('Input new radius value (in \
px): '))
                # Update radius value.
                clust_rad = clust_rad_m
                wrong_answer = False
                flag_radius_manual = True
            else:
                print 'Wrong input. Try again.\n'

    radius_params = [clust_rad, e_rad, flag_delta_total, flag_not_stable,
        flag_delta, flag_radius_manual]

    return radius_params
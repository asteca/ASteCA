"""
@author: gabriel
"""

import numpy as np
from functions.display_rad import disp_rad as d_r
import matplotlib.pyplot as plt


def main_center_algor(rdp_params, cr_params, backg_value):
    '''
    This function holds the main algorithm that returns a radius value.
    '''

    radii, ring_density = rdp_params[:2]
    # Find maximum density value and assume this is the central density.
    max_dens_ind = np.argmax(ring_density)
    ring_dens_c, radii_c = ring_density[max_dens_ind:], radii[max_dens_ind:]

    # Assign a value to the number of points that should be found below
    # the delta value around the background to attain the 'stabilized'
    # condition.
    mode_r = cr_params[0]
    if mode_r not in {'auto', 'manual'}:
        print "  WARNING: CR mode is not valid. Default to 'auto'."
        mode_r = 'auto'
    # Set params.
    if mode_r == 'manual':
        n_left, delta_step = int(cr_params[1]), cr_params[2]
    elif mode_r == 'auto':
        n_left, delta_step = max(int(round(len(radii_c) * 0.1)), 2), 5

    # Difference between max density value and the background value.
    delta_total = (max(ring_dens_c) - backg_value)

    # If the difference between the max density value and the background is
    # less than 3 times the value of the background, raise a flag.
    flag_delta_total = False
    if delta_total < 3 * backg_value:
        flag_delta_total = True

    # Start i value to cap the number of iterations.
    i = 1
    # Initialize condition to break out of 'while'.
    flag_not_stable = True
    # Iterate until a stable condition is attained or for a maximum
    # values of delta_backg.
    while flag_not_stable and i < 6:

        # Store value for delta_percentage.
        delta_percentage = i * delta_step

        # % of difference between max density value and background.
        delta_backg = delta_percentage * delta_total / 100.
        # Increase value of i for next iteration (if it happens)
        i += 1

        # Initialize density values counter for points that fall inside the
        # range determined by the delta value around the background.
        in_delta_val = 0

        # Initialize index_rad value.
        index_rad = 0

        # Iterate through all values of star density in this "square ring".
        for index, item in enumerate(ring_dens_c):

            # Condition to iterate until at least n_left points below the
            # delta + background value are found.
            if in_delta_val < n_left:

                # If the density value is closer than 'delta_backg' to the
                # background value or lower --> add it.
                if (item - backg_value) <= delta_backg:
                    # Augment value of counter.
                    in_delta_val += 1
                    # Store first radius value that falls below the upper delta
                    # limit.
                    if in_delta_val == 1:
                        index_rad = index
            else:
                # Exit 'for' and 'while' loops if n_left consecutive values were
                # found == "stable" condition.
                flag_not_stable = False
                break

    # Raise a flag if the delta used to get the stable condition is greater
    # than 10%.
    flag_delta = False
    if delta_percentage > 10:
        flag_delta = True

    # The first condition is there in case that the stable condition was reached
    # with the last item.
    if (in_delta_val < n_left) and flag_not_stable is True:
        # No radius value found. Assign radius value as the middle element
        # in the radii list.
        clust_rad = radii_c[int(len(radii_c) / 2.)]
    else:
        # Stable condition was attained, assign radius value as the one with
        # the density value closest to the background value among the first
        # n_left points counting from the one determined by the index index_rad.
        radii_dens = [ring_dens_c[index_rad + i] for i in range(n_left)]
        clust_rad = radii_c[index_rad + min(range(len(radii_dens)), key=lambda
        i:abs(radii_dens[i] - backg_value))]

    return clust_rad, delta_backg, delta_percentage, flag_delta_total, \
    flag_not_stable, flag_delta


def get_clust_rad(phot_data, backg_value, cr_params, center_params,
    rdp_params, semi_return, mode):
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
    clust_rad, delta_backg, delta_percentage, flag_delta_total, \
    flag_not_stable, flag_delta = main_center_algor(rdp_params, cr_params,
        backg_value)

    # Check if semi or manual mode are set.
    flag_radius_manual = False
    if mode == 'auto':
        print 'Auto radius found: %0.1f px.' % clust_rad

    elif mode == 'semi':
        # Unpack semi values.
        cent_cl_semi, cl_rad_semi, cent_flag_semi, rad_flag_semi, \
        err_flag_semi = semi_return

        if rad_flag_semi == 1:
            # Update value.
            clust_rad = cl_rad_semi
            print 'Semi radius set: %0.1f px.' % clust_rad
        else:
            print 'Auto radius found: %0.1f px.' % clust_rad

    # If Manual mode is set, display radius and ask the user to accept it or
    # input new one.
    elif mode == 'manual':

        print 'Radius found: %0.1f px' % clust_rad
        d_r(phot_data, center_params, clust_rad, delta_backg, delta_percentage,
            backg_value, rdp_params)
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

    radius_params = [clust_rad, delta_backg, delta_percentage,
    flag_delta_total, flag_not_stable, flag_delta, flag_radius_manual]

    return radius_params
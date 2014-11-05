"""
@author: gabriel
"""

import numpy as np
from functions.display_rad import disp_rad as d_r
import matplotlib.pyplot as plt
import get_in_params as g


def main_rad_algor(rdp_params, field_dens, bin_width):
    '''
    This function holds the main algorithm that returns a radius value.
    '''

    radii, rdp_points = rdp_params[:2]
    # Find maximum density value and assume this is the central density.
    # Do not use previous values.
    max_dens_ind = np.argmax(rdp_points)
    rdp_points_c, radii_c = rdp_points[max_dens_ind:], radii[max_dens_ind:]

    # Assign a value to the number of points that should be found below
    # the delta values around the field density to attain the 'stabilized'
    # condition.
    mode_r = g.cr_params[0]
    if mode_r not in {'auto', 'manual'}:
        print "  WARNING: CR mode is not valid. Default to 'auto'."
        mode_r = 'auto'
    # Set params.
    if mode_r == 'manual':
        # Read the value from input file.
        n_left = int(g.cr_params[1])
    elif mode_r == 'auto':
        # Calculate the value automatically --> 25% of the points in the RDP
        # (min 3 points)
        n_left = max(int(round(len(radii_c) * 0.25)), 3)

    # Fix to a minimum value of 4 to avoid conclict when '(n_left - i)'
    # happens below.
    n_left = max(n_left, 4)

    # Delta step is fixed to 5%.
    delta_step = 5

    # Difference between max RDP density value and the field density value.
    delta_total = (max(rdp_points_c) - field_dens)

    # If the difference between the max density value and the field density is
    # less than 3 times the value of the field density, raise a flag.
    flag_delta_total = False
    if delta_total < 3 * field_dens:
        flag_delta_total = True

    # Initialize index_rad value.
    index_rad = []
    for i in range(4):

        # Store value for delta_percentage --> 20, 15, 10, 5
        delta_percentage = (4. - i) * delta_step / 100.
        # % of difference between max density value and field density.
        delta_field = delta_percentage * delta_total

        # Initialize density values counter for points that fall inside the
        # range determined by the delta value around the field density.
        in_delta_val, index_rad_i, dens_dist = 0, 0, 1.e10

        # Iterate through all values of star density in the RDP.
        for index, item in enumerate(rdp_points_c):

            # Condition to iterate until at least n_left points below the
            # (delta + field density) value are found.
            if in_delta_val < (n_left - i):

                # If the density value is closer than 'delta_field' to the
                # field density value or lower --> add it.
                if item <= delta_field + field_dens:
                    # Augment value of counter.
                    in_delta_val += 1
                    # Store radius value closer to the field density.
                    if abs(item - field_dens) < dens_dist:
                        dens_dist = abs(item - field_dens)
                        index_rad_i = index
                # If the RDP point is outside the (delta + field density) range
                # reset all values.
                else:
                    # Reset.
                    in_delta_val, index_rad_i, dens_dist = 0, 0, 1.e10
            # If enough RDP points have been found within the fiel density
            # range, store the radius value closer to the field density value
            # and break out of the for loop.
            else:
                index_rad.append(index_rad_i)
                break

    # Raise a flag if only two radius values were found under the 4 deltas.
    flag_delta, flag_not_stable = False, False
    if len(index_rad) < 3:
        flag_delta = True

    rad_found = []
    for ind in index_rad:
        # Use the stored indexes to obtain the actual radius values.
        rad_found.append(radii_c[ind])
    # If at least one radius value was found.
    if rad_found:
        clust_rad, e_rad = np.mean(rad_found), max(np.std(rad_found), bin_width)
    else:
        flag_not_stable = True
        # No radius value found. Assign radius value as the middle element
        # in the radii list.
        clust_rad, e_rad = radii_c[int(len(radii_c) / 2.)], 0.
        print '  WARNING: no radius found, setting value to: {:g} {}'.format(
            clust_rad, g.gd_params[0][-1])

    return clust_rad, e_rad, flag_delta_total, flag_not_stable, flag_delta


def get_clust_rad(phot_data, field_dens, center_params, semi_return, bin_width,
    rdp_params):
    """
    Obtain the value for the cluster's radius by counting the number of points
    that fall within a given interval of the field density or lower. If this
    number is equal to a minimum fixed number of points 'n_left', then assign
    the radius as the point closest to the field density value among those
    first n_left points counting from the first one that fell below the
    (field dens + delta limit) range.
    Iterate increasing the interval around the field density and finally
    average all the radius values found for each interval.
    """

    # Call function that holds the radius finding algorithm.
    clust_rad, e_rad, flag_delta_total, flag_not_stable, flag_delta = \
    main_rad_algor(rdp_params, field_dens, bin_width)

    # Check if semi or manual mode are set.
    flag_radius_manual = False
    if g.mode == 'auto':
        print 'Auto radius found: {:g} {}.'.format(clust_rad,
            g.gd_params[0][-1])

    elif g.mode == 'semi':
        # Unpack semi values.
        cent_cl_semi, cl_rad_semi, cent_flag_semi, rad_flag_semi, \
        err_flag_semi = semi_return

        if rad_flag_semi == 1:
            # Update values.
            clust_rad, e_rad = cl_rad_semi, 0.
            print 'Semi radius set: {:g} {}.'.format(clust_rad,
                g.gd_params[0][-1])
        else:
            print 'Auto radius found: {:g} {}.'.format(clust_rad,
                g.gd_params[0][-1])

    # If Manual mode is set, display radius and ask the user to accept it or
    # input new one.
    elif g.mode == 'manual':

        print 'Radius found: {:g} {}.'.format(clust_rad, g.gd_params[0][-1])
        d_r(phot_data, bin_width, center_params, clust_rad, field_dens,
            rdp_params)
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
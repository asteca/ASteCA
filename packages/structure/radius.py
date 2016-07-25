
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import warnings
from ..inp import input_params as g
import display_rad
from ..out import prep_plots


def radius_algor(rdp_params, field_dens, bin_width, coord):
    '''
    This function holds the main algorithm that returns a radius value.
    '''

    radii, rdp_points = rdp_params[:2]
    # Find maximum density value and assume this is the central density.
    # Do not use previous values.
    max_dens_ind = np.argmax(rdp_points)
    rdp_points_m, radii_m = rdp_points[max_dens_ind:], radii[max_dens_ind:]

    # Interpolate extra RDP points between those calculated.
    N = 1000
    t, xp = np.linspace(0, 1, N), np.linspace(0, 1, len(rdp_points_m))
    # Store interpolated RDP and radii values.
    rdp_points_c = np.interp(t, xp, rdp_points_m)
    radii_c = np.interp(t, xp, radii_m)

    # Assign a value to the number of points that should be found below
    # the delta values around the field density to attain the 'stabilized'
    # condition.
    mode_r = g.cr_params[0]
    if mode_r == 'low':
        # Fix to 5% of the total number of interpolated points in the RDP.
        n_left = int(0.05 * N)
    elif mode_r == 'mid':
        # Fix to 10% of the total number of interpolated points in the RDP.
        n_left = int(0.1 * N)
    elif mode_r == 'high':
        # Fix to 20% of the total number of interpolated points in the RDP.
        n_left = int(0.2 * N)

    # Difference between max RDP density value and the field density value.
    delta_total = (max(rdp_points_c) - field_dens)

    # If the difference between the max density value and the field density is
    # less than 3 times the value of the field density, raise a flag.
    flag_delta_total = False
    if delta_total < 3 * field_dens:
        flag_delta_total = True

    # Initialize index_rad value.
    rad_found = []
    for k, delta_percentage in enumerate(np.arange(0.2, 0.1, -0.01)):
        # The 'k' param relaxes the condition that requires a certain number of
        # points to be located below the 'delta_field + field_dens' value
        # to establish that the RDP has stabilized.
        # The smaller the value of 'delta_field', the fewer the number of
        # consecutive points that are required.

        # delta_field = % of difference between max density value and field
        # density.
        delta_field = delta_percentage * delta_total

        # Initialize density values counter for points that fall inside the
        # range determined by the delta value around the field density.
        in_delta_val, index_rad_i, dens_dist = 0, 0, 1.e10

        # Iterate through all values of star density in the RDP.
        for index, dens in enumerate(rdp_points_c):

            # Condition to iterate until at least in_delta_val *consecutive*
            # points below the (delta + field density) value are found.
            if in_delta_val < (n_left - (4 * k)):

                # If the density value is closer than 'delta_field' to the
                # field density value or lower --> add it.
                if dens <= delta_field + field_dens:
                    # Augment value of counter.
                    in_delta_val += 1
                    # Store radius value closer to the field density.
                    if abs(dens - field_dens) < dens_dist:
                        dens_dist = abs(dens - field_dens)
                        index_rad_i = index
                # If the RDP point is outside the (delta + field density) range
                # reset all values. I.e.: the points should be located below
                # the 'delta_field + field_dens' *consecutively*.
                else:
                    # Reset.
                    in_delta_val, index_rad_i, dens_dist = 0, 0, 1.e10
            # If enough RDP points have been found within the field density
            # range, store the radius value closer to the field density value
            # and break out of the for loop.
            else:
                rad_found.append(radii_c[index_rad_i])
                break

    # Raise a flag if only two radius values were found under the deltas.
    flag_delta, flag_not_stable = False, False
    if len(rad_found) < 3:
        flag_delta = True

    # If at least one radius value was found.
    if rad_found:
        # Use the median to avoid outliers.
        clust_rad = np.median(rad_found)

        # Obtain error as the 1 sigma confidence interval (68.27%).
        _std = np.std(rad_found)
        # Catch warning if stats fails to obtain confidence interval.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            conf_int = stats.norm.interval(0.6827, loc=clust_rad, scale=_std /
                                           np.sqrt(len(rad_found)))
        # If stats returns a confidence interval with a NaN, discard it.
        if np.any(np.isnan(conf_int)):
            conf_dif = 0.
        else:
            conf_dif = abs(clust_rad - max(conf_int))
        e_rad = max(conf_dif, bin_width)
        # Prevent too small radius by fixing the minimum value to the second
        # RDP point.
        if clust_rad < radii[1]:
            clust_rad = radii[1]
    else:
        flag_not_stable = True
        # No radius value found. Assign radius value as the middle element
        # in the radii list.
        clust_rad, e_rad = radii_c[int(len(radii_c) / 2.)], 0.
        print('  WARNING: no radius found, setting value to: {:g} {}'.format(
            clust_rad, coord))

    return clust_rad, e_rad, flag_delta_total, flag_not_stable, flag_delta


def main(phot_data, field_dens, center_params, rdp_params,
         semi_return, bin_width):
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

    coord = prep_plots.coord_syst()[0]
    # Call function that holds the radius finding algorithm.
    clust_rad, e_rad, flag_delta_total, flag_not_stable, flag_delta = \
        radius_algor(rdp_params, field_dens, bin_width, coord)

    # Check if semi or manual mode are set.
    flag_radius_manual = False
    if g.mode == 'auto':
        print('Auto radius found: {:g} {}.'.format(clust_rad, coord))

    elif g.mode == 'semi':
        # Unpack semi values.
        cl_rad_semi, rad_flag_semi = semi_return[1], semi_return[4]

        if rad_flag_semi == 1:
            # Update values.
            clust_rad, e_rad = cl_rad_semi, 0.
            print('Semi radius set: {:g} {}.'.format(clust_rad, coord))
        else:
            print('Auto radius found: {:g} {}.'.format(clust_rad, coord))

    # If Manual mode is set, display radius and ask the user to accept it or
    # input new one.
    elif g.mode == 'manual':

        print 'Radius found: {:g} {}.'.format(clust_rad, coord)
        display_rad.main(phot_data, bin_width, center_params, clust_rad,
                         e_rad, field_dens, rdp_params)
        plt.show()

        # Ask if the radius is accepted, or a if a another one should be used.
        while True:
            answer_rad = raw_input('Input new radius value? (y/n) ')
            if answer_rad == 'n':
                print('Value accepted.')
                break
            elif answer_rad == 'y':
                try:
                    clust_rad_m = float(raw_input('cluster_rad: '))
                    # Update radius value.
                    clust_rad = clust_rad_m
                    flag_radius_manual = True
                    break
                except:
                    print("Sorry, input is not valid. Try again.")
            else:
                print("Sorry, input is not valid. Try again.\n")

    radius_params = [clust_rad, e_rad, flag_delta_total, flag_not_stable,
                     flag_delta, flag_radius_manual]

    return radius_params

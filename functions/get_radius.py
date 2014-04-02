"""
@author: gabriel
"""


def get_clust_rad(backg_value, radii, ring_density, cr_params):
    """
    Obtain the value for the cluster's radius by counting the number of points
    that fall within a given interval of the background or lower. If this number
    is equal to a fixed number of points n_left then assign the radius as the
    closest point to the background value among the first n_left points
    counting from the first one that fell below the backg + delta limit.
    Iterate increasing the interval around the background until n_left points
    are found or the delta interval reaches its maximum allowed.
    """

    # Assign a value to the number of points that should be found below
    # the delta value around the background to attain the 'stabilized'
    # condition.
    if cr_params[0] == 'auto':
        n_left = max(int(round(len(radii) * 0.1)), 2)
    else:
        n_left = cr_params[1]

    # Difference between max density value and the background value.
    delta_total = (max(ring_density) - backg_value)

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
        delta_percentage = i * 5

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
        for index, item in enumerate(ring_density):

            # Condition to iterate until at least n_left points below the
            # delta + background value are found.
            if in_delta_val < n_left:

                # If the density value is closer than 'delta_backg' to the
                # background value or lower --> add it.
                # The condition index!=0 is there to avoid counting the first
                # point which sometimes presents a very low density value due
                # to the small width of the bin used.
                if (item - backg_value) <= delta_backg and index != 0:
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
        clust_rad = radii[int(len(radii) / 2.)]
    else:
        # Stable condition was attained, assign radius value as the one with
        # the density value closest to the background value among the first
        # n_left points counting from the one determined by the index index_rad.
        radii_dens = [ring_density[index_rad + i] for i in range(n_left)]
        clust_rad = radii[index_rad + min(range(len(radii_dens)), key=lambda
        i:abs(radii_dens[i] - backg_value))]

    radius_params = [clust_rad, delta_backg, delta_percentage,
    flag_delta_total, flag_not_stable, flag_delta]

    return radius_params
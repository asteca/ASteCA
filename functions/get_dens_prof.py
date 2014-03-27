"""
@author: gabriel
"""

import math


def get_dens_prof(histo, x_c_b, y_c_b, width_bins, inner_ring):
    """
    Calculate the density profile by counting the number of stars in the center
    bin first (r aprox width_bins/2 px), then moving to the 8 adyacent bins
    (r aprox width_bins + (width_bins/2) px), then the following 16 (r aprox
    2*width_bins + (width_bins/2) px), then the next 24 (r aprox 3*width_bins +
    (width_bins/2) px) and so forth. These "square rings" have consecutive
    multiples of 8 bins each (if no bin falls outside the range of the x,y
    chart) and an area of N*width_bin^2 [px^2], where N is the number of bins.
    Dividing the number of stars in each "square ring" by this area, we get all
    the "ring densities" for those approximate values of r.
    """

    # Initialize all lists with len(width_bins) items all set to zero.
    radii, ring_density, ring_count, bin_count, poisson_error = \
    [0] * len(width_bins), [0] * len(width_bins), [0] * len(width_bins), \
    [0] * len(width_bins), [0] * len(width_bins)

    # Iterate through all 2D histograms obtained with different bin widths
    for hist_index, hist_2d_item in enumerate(histo):

        # Number of points in lists (stops at inner_ring px)
        points_list = int(inner_ring / width_bins[hist_index]) + 1
        # Initialize all empty lists inside each list belonging to each
        # iteration.
        ring_density[hist_index], ring_count[hist_index], \
        bin_count[hist_index], poisson_error[hist_index] = [0] * points_list,\
        [0] * points_list, [0] * points_list, [0] * points_list

        # Iterate through all the "square rings" up to inner_ring px (where the
        # background suposedly starts)
        for i in range(points_list):
            # Iterate through items in the x dimension for the 2D hist
            for xindex, xitem in enumerate(hist_2d_item):
                # Iterate through items in the y dimension for the 2D hist
                for yindex, yitem in enumerate(xitem):
                    # Bin is in the top row
                    if yindex == (y_c_b[hist_index] + i) and \
                    (x_c_b[hist_index] - i) <= xindex and \
                    xindex <= (x_c_b[hist_index] + i):
                        # Add stars in bin to corresponding ring
                        ring_count[hist_index][i] = ring_count[hist_index][i] +\
                        yitem
                        # Add 1 more bin to the "square ring".
                        bin_count[hist_index][i] += 1
                    # Bin is in the bottom row
                    elif yindex == (y_c_b[hist_index] - i) and \
                    (x_c_b[hist_index] - i) <= xindex and \
                    xindex <= (x_c_b[hist_index] + i):
                        # Add stars in bin to corresponding ring
                        ring_count[hist_index][i] = ring_count[hist_index][i] +\
                        yitem
                        # Add 1 more bin to the "square ring".
                        bin_count[hist_index][i] += 1
                    # Bin is in the left column
                    elif xindex == (x_c_b[hist_index] - i) and \
                    (y_c_b[hist_index] - (i - 1)) <= yindex and \
                    yindex <= (y_c_b[hist_index] + (i - 1)):
                        # Add stars in bin to corresponding ring
                        ring_count[hist_index][i] = ring_count[hist_index][i] +\
                        yitem
                        # Add 1 more bin to the "square ring"
                        bin_count[hist_index][i] += 1
                    # Bin is in the right column
                    elif xindex == (x_c_b[hist_index] + i) and \
                    (y_c_b[hist_index] - (i - 1)) <= yindex and \
                    yindex <= (y_c_b[hist_index] + (i - 1)):
                        # Add stars in bin to corresponding ring
                        ring_count[hist_index][i] = ring_count[hist_index][i] +\
                        yitem
                        # Add 1 more bin to the "square ring"
                        bin_count[hist_index][i] += 1

            # If no stars are inside the bins, set value to 1 to avoid
            # division by zero.
            if bin_count[hist_index][i] == 0:
                bin_count[hist_index][i] = 1
            # The number of bins times the area of each bin gives the density
            # area.
            dens_area = bin_count[hist_index][i] * (width_bins[hist_index] ** 2)

            # Calculate density corresponding to "square ring" i
            ring_density[hist_index][i] = ring_count[hist_index][i] / dens_area
            # Obtain the Poisson error bar for each value
            poisson_error[hist_index][i] = \
            math.sqrt(ring_count[hist_index][i]) / dens_area

        # Store values for radii to go with the densities obtained above
        # and stored in 'ring_density'
        radii[hist_index] = [width_bins[hist_index] / 2. +
        (width_bins[hist_index] * i) for i in range(points_list)]

    # Return 4 lists of lists, ie: radii = [[12.5, 37.5, ..], [25.0, 75.0, ..],\
    # [37.5, ..], [50.0, ..]]
    return radii, ring_density, poisson_error
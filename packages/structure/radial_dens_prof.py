
import math


def main(clp):
    """
    Calculate the density profile by counting the number of stars in the center
    bin first (r aprox width_bin/2 px), then moving to the 8 adjacent bins
    (r aprox width_bin + (width_bin/2) px), then the following 16 (r aprox
    2*width_bin + (width_bin/2) px), then the next 24 (r aprox 3*width_bin +
    (width_bin/2) px) and so forth. These "square rings" have consecutive
    multiples of 8 bins each (if no bin falls outside the range of the x,y
    chart) and an area of N*width_bin^2 [px^2], where N is the number of bins.
    Dividing the number of stars in each "square ring" by this area, we get all
    the "ring densities" for those approximate values of r.
    """

    x_c_b, y_c_b = clp['bin_cent']

    square_rings, radii, rdp_points, poisson_error = [], [], [], []
    # Use max x,y length defined in the 2D histogram.
    rdp_length = max(len(clp['hist_2d']), len(clp['hist_2d'][0]))
    # Iterate through all the bins in the largest dimension.
    for i in range(rdp_length):
        # Store here the coordinates of the bins.
        bins_coords = []

        # Initialize bin_count for this square ring.
        ring_count, bin_count = 0, 0

        # Iterate through bins in the x dimension for the 2D hist.
        for xindex, xitem in enumerate(clp['hist_2d']):
            # Iterate through bins in the y dimension for the 2D hist.
            for yindex, st_in_bin in enumerate(xitem):

                # Bin is in the top row
                if yindex == (y_c_b + i) and abs(xindex - x_c_b) <= i:
                    # Add stars in bin to corresponding ring.
                    ring_count += st_in_bin
                    # Add 1 more bin to the "square ring".
                    bin_count += 1
                    bins_coords.append([xindex - x_c_b, i])
                # Bin is in the bottom row
                elif yindex == (y_c_b - i) and abs(xindex - x_c_b) <= i:
                    ring_count += st_in_bin
                    bin_count += 1
                    bins_coords.append([xindex - x_c_b, -i])
                # Bin is in the left column
                elif xindex == (x_c_b - i) and abs(yindex - y_c_b) <= (i - 1):
                    ring_count += st_in_bin
                    bin_count += 1
                    bins_coords.append([-i, yindex - y_c_b])
                # Bin is in the right column
                elif xindex == (x_c_b + i) and abs(yindex - y_c_b) <= (i - 1):
                    ring_count += st_in_bin
                    bin_count += 1
                    bins_coords.append([i, yindex - y_c_b])

        # Break when no more bins are stored in this square ring. This means
        # we reached the border of the frame.
        if bin_count == 0:
            break

        # Store bin coordinates in each square ring.
        square_rings.append(bins_coords)
        # If no stars are inside this square ring, set value to 1 to avoid a
        # division by zero.
        bin_count = 1 if bin_count == 0 else bin_count
        # The number of bins times the area of each bin gives the area of
        # this square ring.
        area = bin_count * (clp['bin_width'] ** 2)

        # Calculate density corresponding to "square ring" i
        rdp_points.append(ring_count / area)
        # Obtain the Poisson error bar for each value
        poisson_error.append(math.sqrt(ring_count) / area)

        # Store values for radii to go with the densities obtained above
        # and stored in 'rdp_points'
        radii.append(clp['bin_width'] / 2. + (clp['bin_width'] * i))

    # Transform from bin units to coordinate units before passing.
    rdp_length = rdp_length * clp['bin_width']

    if rdp_points:
        print("Radial density profile (RDP) calculated")
    else:
        raise ValueError("ERROR: RDP is empty. Check the center coordinates")

    # Add data to dictionary.
    clp['radii'], clp['rdp_points'], clp['poisson_error'],\
        clp['square_rings'], clp['rdp_length'] = radii, rdp_points,\
        poisson_error, square_rings, rdp_length
    return clp

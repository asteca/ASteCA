"""
@author: gabriel
"""


def get_background(x_data, y_data, x_c_b, y_c_b, h_not_filt, width_bins,
                   inner_ring, outer_ring):
    """
    Get background level of stars by calculating the density
    of (num of stars)/area for a "square ring" centered at the cluster's center
    and whose minor and major radii are: inner_ring and outer_ring pixels
    respectively. We use "square rings" instead of circles because that way we
    can avoid using areas that are too big in the case where the cluster might
    be near a border and the raddi fall outside the boundaries of the frame.
    """

    # Use only histogram made wiht the smallest bin width.
    index, item = 0, h_not_filt[0]

    # If the bin is at least inner_ring px away from the cluster's center
    # and less than outer_ring px, then use it to calculate the background
    # value. We set an outer boundary to try to prevent empty areas of the
    # frame from artificially lowering the final density value.
    inner_bound, outer_bound = int(inner_ring / width_bins[index]), \
    int(outer_ring / width_bins[index])

    # Initialize bin and star counters
    bin_count, star_count = 0, 0

    # Iterate through items in the x dimension for the non-filtered 2D hist.
    for xindex, xitem in enumerate(item):
        # Iterate through items in the y dimension for the non-filtered
        # 2D hist.
        for yindex, yitem in enumerate(xitem):

            # Two columns bounded in y axis and in the x axis by the outer
            # limit.
            if (abs(yindex - y_c_b[index]) < outer_bound) and \
            (abs(yindex - y_c_b[index]) > inner_bound)\
            and (abs(xindex - x_c_b[index]) < outer_bound):
                # Add stars in bin to total star count
                star_count = star_count + yitem
                # Increase number of bins (ie: increase area covered)
                bin_count += 1

            # Two columns bounded in x axis and in the y axis by the inner
            # limit.
            if (abs(xindex - x_c_b[index]) > inner_bound) and \
            (abs(xindex - x_c_b[index]) < outer_bound)\
            and (abs(yindex - y_c_b[index]) < inner_bound):
                # Add stars in bin to total star count
                star_count = star_count + yitem
                # Increase number of bins (ie: increase area covered)
                bin_count += 1

    # Calculate the background value by dividing the number of stars found
    # away from the cluster's center by the total area
    # (where width_bins[index]**2 is the area of a single bin)
    flag_bin_count = False
    if bin_count != 0:
        area = bin_count * (width_bins[index] ** 2)
        backg_value = star_count / area
    else:
        # A value of bin_count=0 for a given 2D histogram means that the
        # inner limit pushed the bins outside the frame. This could be
        # indicative of a small field given the size of the cluster.
        flag_bin_count = True
        backg_value = 0.

    # Return background value obtained.
    return backg_value, flag_bin_count


from ..inp import input_params as g
import spiral as sp
import manual_histo


def spiral_index(spiral, sp_indx, histo, x_c_b, y_c_b, num_bins_area):
    '''
    Take the fixed x,y coordinates for a squared spiral of bins (spiral)
    centered at [0,0] and an index (sp_indx) that points to a given
    coordinate (where 0 means center coords: [0,0]).

    Center the spiral at the coordinates (x_c_b, y_c_b) for a 2D histogram
    (histo) and loop through the spiral storing the coordinates of each
    histogram bin until a given total area (num_bins_area) is obtained.
    '''

    # Initialize the bin counter that indicates how many bins are already
    # added to the region.
    bin_count = 0
    # Initialize empty list that will hold the coordinates of the spiral bins.
    sp_coords = [[], []]

    # Loop spiral.
    for indx, sp_item in enumerate(spiral[sp_indx:]):
        # Loop until region is filled.
        if bin_count <= num_bins_area:

            # Check if the bin exists in the 2D histogram.
            try:
                histo[x_c_b + sp_item[0]][y_c_b + sp_item[1]]
            except IndexError:
                pass  # Item out of histogram range.
            else:
                # Check that the index is not negative because python
                # will assign items in lists even if they are pointed
                # as negative values, ie: list1[1][-2]; which in this
                # case makes no sense because it would fall out of the
                # 2D histogram.
                if (x_c_b + sp_item[0]) >= 0 and (y_c_b + sp_item[1]) >= 0:
                    # If the item exists, we store the coordinates of
                    # that bin in this region in both coordinates.
                    sp_coords[0].append(x_c_b + sp_item[0])
                    sp_coords[1].append(y_c_b + sp_item[1])
                    # Increase the bin count.
                    bin_count += 1
                    # Store the index of the last bin.
                    sp_indx2 = indx + sp_indx

    return sp_indx2, sp_coords


def spiral_region(h_manual, sp_coords):
    '''
     At this point we have the list 'sp_coords' composed of two lists, the
     first one containing the x coordinates for every bin that corresponds
     to the region and the second list the y coordinates. We need to
     obtain the stars located inside each of those bins. To do this
     we use 'h_manual' which is a list that already contains the stars that
     fall inside each of the bins in the 2D histogram along with their
     relevant data (ID, x, y, mag, etc..)
    '''

    # Initialize empty field region.
    f_region = []

    # Iterate through all the bins in the spiral list defined.
    for xbin, ybin in zip(*sp_coords):

        if h_manual[xbin][ybin][0] > 0:
            # Add all the stars inside this bin to the region.
            # We use [1:] to skip the first item which holds the
            # number of stars in the bin.
            f_region.extend(h_manual[xbin][ybin][1:])

    return f_region


def main(semi_return, hist_lst, cent_bin, clust_rad, cl_area,
         stars_out):
    '''
    Define empty region around the cluster via a spiral centered on it
    and of area a bit larger than that defined by the cluster's radius.

    Define as many field regions as requested (or permitted) around this empty
    region.
    '''

    hist_2d, xedges, yedges, bin_width = hist_lst
    x_c_b, y_c_b = cent_bin

    # Number of field regions defined in 'params_input.dat' file.
    f_regs_num = g.fr_number

    # Check if semi is set.
    if g.mode == 'semi':
        # Unpack semi values.
        cl_f_regs_semi, freg_flag_semi = semi_return[2], semi_return[5]

        if freg_flag_semi == 1:
            # Update value.
            f_regs_num = cl_f_regs_semi

    # Get area as total number of bins in 2D hist times the area of each bin.
    total_area = len(hist_2d[0]) * len(hist_2d) * (bin_width ** 2)
    # Define empty area around the cluster region.
    sq_area = 2. * cl_area
    # Maximum number of field regions possible.
    f_regs_max = int((total_area - sq_area) / cl_area)

    # If the remaining area in the frame after subtracting the cluster region
    # is smaller than the cluster's area, this means that the cluster is either
    # too large or the frame too small and no field region of equal area than
    # that of the cluster can be obtained.
    # Raise a flag.
    flag_no_fl_regs = False
    if f_regs_max < 1:
        print ("  WARNING: cluster region is too large or frame\n"
               "  is too small. No field regions available.")
        flag_no_fl_regs = True
    else:
        # If the number of field regions defined is larger than the maximum
        # allowed, use the maximum.
        if f_regs_num == 'max':
            f_regions = f_regs_max
            print 'Using max number of field regions ({}).'.format(f_regions)
        elif f_regs_num > f_regs_max:
            f_regions = f_regs_max
            print ("  WARNING: Number of FR defined ({}) is larger than\n"
                   "  the maximum allowed ({}). "
                   "Using max number.").format(f_regs_num, f_regs_max)
        elif f_regs_num <= 0:
            f_regions = f_regs_max
            print ("  WARNING: Number of FR ({}) is less than or equal\n"
                   "  to zero. No field region will be defined.").format(
                f_regs_num)
            flag_no_fl_regs = True
        else:
            print ("Using defined number of field regions ({}).".format(
                f_regs_num))
            f_regions = f_regs_num

    # Get list that contains the spiral as a list of x,y coordinates (also
    # stored as lists) starting from the initial bin [0, 0].
    spiral = sp.main()

    # Calculate number of bins such that their combined area equals the
    # larger area around the cluster defined above.
    num_bins_area = int(round(sq_area / (bin_width ** 2), 0))
    # Obtain index of spiral bin where field regions should begin to be formed.
    # dummy_lst is not important here.
    sp_indx, dummy_lst = spiral_index(spiral, 0, hist_2d, x_c_b, y_c_b,
                                      num_bins_area)

    # Obtain field regions only if it is possible.
    # This list holds all the field regions.
    field_regions = []
    if flag_no_fl_regs is False:

        # Obtain filled 2D histogram for the field with star's values attached
        # to each bin.
        h_manual = manual_histo.main(stars_out, xedges, yedges)

        # This ensures that the areas of the field regions are equal
        # to the cluster area.
        num_bins_area = int(cl_area / (bin_width ** 2))

        for _ in range(f_regions):
            # Retrieve spiral index where this field region should end and
            # coordinates of its bins.
            sp_indx, sp_coords = spiral_index(spiral, sp_indx, hist_2d, x_c_b,
                                              y_c_b, num_bins_area)
            # Fill spiral section for this field region with all the stars
            # that fall inside of it.
            f_region = spiral_region(h_manual, sp_coords)
            field_regions.append(f_region)

        # If any of the field regions has less than 4 stars then we remove it
        # from the list otherwise the decontamination or the p-value algorithms
        # will fail. This list stores the indexes of the empty regions.
        field_regs_del = []
        for indx, s_lst in enumerate(field_regions):
            if len(s_lst) < 4:
                field_regs_del.append(indx)
        # Delete empty regions this way to avoid messing with the indexes.
        for index in sorted(field_regs_del, reverse=True):
            del field_regions[index]
        if field_regs_del:
            print ('  {} field regions with less than 4 stars each were'
                   ' removed.').format(len(field_regs_del))

        # If after removing the empty regions no regions are left, raise the
        # flag.
        if f_regions > 0 and not(field_regions):
            print ('  WARNING: no field regions left after the removal of\n' +
                   '  those containing less than 4 stars.')
            flag_no_fl_regs = True

    return flag_no_fl_regs, field_regions

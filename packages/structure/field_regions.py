
import numpy as np
import spiral as sp
import field_manual_histo


def spiral_index(spiral, sp_indx, histo, cent_bin, num_bins_area):
    '''
    Take the fixed x,y coordinates for a squared spiral of bins (spiral)
    centered at [0,0] and an index (sp_indx) that points to a given
    coordinate (where 0 means center coords: [0,0]).

    Center the spiral at the coordinates (cent_bin[0], cent_bin[1]) for a
    2D histogram (histo) and loop through the spiral storing the coordinates
    of each histogram bin until a given total area (num_bins_area) is obtained.
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
                histo[cent_bin[0] + sp_item[0]][cent_bin[1] + sp_item[1]]
            except IndexError:
                pass  # Item out of histogram range.
            else:
                # Check that the index is not negative because python
                # will assign items in lists even if they are pointed
                # as negative values, ie: list1[1][-2]; which in this
                # case makes no sense because it would fall out of the
                # 2D histogram.
                if (cent_bin[0] + sp_item[0]) >= 0 and\
                        (cent_bin[1] + sp_item[1]) >= 0:
                    # If the item exists, we store the coordinates of
                    # that bin in this region in both coordinates.
                    sp_coords[0].append(cent_bin[0] + sp_item[0])
                    sp_coords[1].append(cent_bin[1] + sp_item[1])
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


def fieldRegs(hist_2d, bin_width, cl_area):
    """
    Estimate the maximum number of field regions that can possibly be defined,
    and the number of bins whose combined area equals the cluster region
    plus the empty region around it.
    """
    # Number of bins in x and y.
    x_bins, y_bins = len(hist_2d[0][0]), len(hist_2d[0])
    # Total area: total number of bins in 2D hist times the area of each bin.
    total_area = x_bins * y_bins * (bin_width ** 2)

    # Several empty areas, in decreasing order of size.
    sq_areas = np.arange(2., 1.2, -.05) * cl_area
    # All the possible total number of field regions.
    f_regs_all = [int(_) for _ in (total_area - sq_areas) / cl_area]

    # Find the first index where the number of field regions is >= 1. If
    # none is found, then the index '0' is returned, as desired.
    i = np.argmax(np.array(f_regs_all) >= 1)
    sq_area, f_regs_max = sq_areas[i], f_regs_all[i]

    # Number of bins such that their combined area equals the
    # larger 'sq_area' area around the cluster.
    num_bins_sqarea = int(round(sq_area / (bin_width ** 2), 0))

    return num_bins_sqarea, f_regs_max


def main(i_c, clp, run_mode, fr_number, cl_f_regs_semi, freg_flag_semi,
         **kwargs):
    '''
    Define empty region around the cluster via a spiral centered on it
    and of area a bit larger than that defined by the cluster's radius.

    Define as many field regions as requested (or permitted) around this empty
    region.
    '''
    # Number of field regions defined in 'params_input.dat' file.
    f_regs_num = fr_number

    # Check if semi is set.
    if run_mode == 'semi':
        if freg_flag_semi == 1:
            # Update value.
            f_regs_num = cl_f_regs_semi

    num_bins_sqarea, f_regs_max = fieldRegs(
        clp['hist_2d'], clp['bin_width'], clp['cl_area'])

    # If the maximum number of field regions that can be formed is zero, it
    # means that the cluster is either too large or the frame too small.
    flag_no_fl_regs = False
    if f_regs_max < 1:
        print ("    WARNING: cluster region is too large.\n"
               "    No field regions available.")
        flag_no_fl_regs = True
    else:
        # If the number of field regions defined is larger than the maximum
        # allowed, use the maximum.
        if f_regs_num == 'max':
            f_regions = f_regs_max
            print('  Using maximum number of field regions ({}).'.format(
                f_regions))
        elif f_regs_num > f_regs_max:
            f_regions = f_regs_max
            print("    WARNING: Number of FR defined ({}) is larger than\n"
                  "    the maximum allowed ({}). Using max number.").format(
                      f_regs_num, f_regs_max)
        elif f_regs_num <= 0:
            f_regions = f_regs_max
            print("    WARNING: Number of FR ({}) is less than or equal\n"
                  "    to zero. No field region will be defined.").format(
                f_regs_num)
            flag_no_fl_regs = True
        else:
            print("  Using defined number of field regions ({}).".format(
                f_regs_num))
            f_regions = f_regs_num

    # Obtain field regions only if it is possible.
    field_regions = []
    if flag_no_fl_regs is False:

        # List that contains the spiral as a list of x,y coordinates (also
        # stored as lists) starting from the initial bin [0, 0].
        spiral = sp.main()

        # Index of spiral bin where field regions should begin to be
        # formed. dummy list is not important here.
        sp_indx, dummy = spiral_index(
            spiral, 0, clp['hist_2d'][0], clp['bin_cent'], num_bins_sqarea)

        # Obtain filled 2D histogram for the field with star's values attached
        # to each bin.
        h_manual = field_manual_histo.main(
            clp['stars_out_' + i_c[0]], clp['xedges'], clp['yedges'])

        # This ensures that the areas of the field regions are equal
        # to the cluster area.
        num_bins_area = int(clp['cl_area'] / (clp['bin_width'] ** 2))

        for _ in range(f_regions):
            # Retrieve spiral index where this field region should end and
            # coordinates of its bins.
            sp_indx, sp_coords = spiral_index(
                spiral, sp_indx, clp['hist_2d'][0], clp['bin_cent'],
                num_bins_area)
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
            print('    {} field regions with less than 4 stars each were'
                  ' removed.').format(len(field_regs_del))

        # If after removing the empty regions no regions are left, raise the
        # flag.
        if f_regions > 0 and not(field_regions):
            print('    WARNING: no field regions left after the removal of\n' +
                  '    those containing less than 4 stars.')
            flag_no_fl_regs = True

    if i_c == 'incomp':
        clp['flag_no_fl_regs_i'], clp['field_regions_i'] = flag_no_fl_regs,\
            field_regions
    elif i_c == 'comp':
        clp['flag_no_fl_regs_c'], clp['field_regions_c'] = flag_no_fl_regs,\
            field_regions

    return clp

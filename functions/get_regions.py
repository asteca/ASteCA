# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 10:06:31 2013

@author: gabriel
"""

from functions.get_spiral import spiral as gs
import numpy as np


def spiral_region(histo, h_manual, stars_in, stars_out, x_c_b, y_c_b, spiral,
                  num_bins_area, sp_indx):
    '''
    Define a spiral region starting from a certain coordinate until a
    given area is obtained.
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

    # At this point we have the list 'sp_coords' composed of two lists, the
    # first one containing the x coordinates for every bin that corresponds
    # to the region and the second list the y coordinates. We need to
    # obtain the stars located insiden each of those bins. To do this
    # we use 'h_manual' wich is a list that already contains the stars that
    # fall inside each of the bins in the 2D histogram along with their
    # relevant data (ID, x, y, T1, etc..)

    # Initialize empty region.
    region = []

    # Iterate through all the bins in the spiral list defined.
    for xbin_index, xbin in enumerate(sp_coords[0]):
        ybin = sp_coords[1][xbin_index]

        for star in range(h_manual[xbin][ybin][0]):
            # Add all the stars inside this bin to the region.
            # We use 'star+1' to skip the first item which holds the
            # number of stars in the bin.
            # Check to see if star belongs to the group of stars
            # that were NOT removed because of its large error.
            # h_manual[xbin][ybin][star+1][0] is the ID of the star.
            st_id = h_manual[xbin][ybin][star + 1][0]
            st_x = h_manual[xbin][ybin][star + 1][1]
            if any(st_id == i[0] for i in stars_in) and \
            any(st_x == i[1] for i in stars_in) or \
            any(st_id == i[0] for i in stars_out) and \
            any(st_x == i[1] for i in stars_out):
                region.append(h_manual[xbin][ybin][star + 1])

    return region, sp_indx2


def get_regions(bin_center, bin_width, histo, clust_rad, h_manual, stars_in,
    stars_out, gr_params):
    '''
    Define cluster and field regions around the cluster's center.
    '''

    x_c_b, y_c_b = bin_center
    # Maximum number of field regions to attempt to fill.
    f_regions = gr_params[0]

    # Define region around the cluster as a spiral centered in it
    # and of area a bit larger than that defined by the cluster's radius.

    # Get area as total number of bins in 2D hist times the area of each bin.
    area = len(histo[0]) * len(histo) * (bin_width ** 2)

    # Length of the side of the square that contains the cluster. Increase it
    # if the radius is comparable to the bin width used to make sure the
    # region covers the entire cluster.
    length = 2.5 if clust_rad > 2 * bin_width else 5.

    # If the remaining area in the frame after substracting the cluster region
    # is smaller than the cluster's area, this means that the cluster is either
    # too large or the frame too small and no field region of equal area than
    # that of the cluster can be obtained.
    # Raise a flag.
    flag_area_stronger = False
    if (area - (length * clust_rad) ** 2) < np.pi * clust_rad ** 2:
        print '  WARNING: cluster region too large, no field region available.'
        flag_area_stronger = True
    else:
        # Calculate maximum number of field regions possible. The minus 1
        # at the end is there so the last field region is defined a bit
        # removed from the edges of the frame.
        f_regs_max = int((area - (length * clust_rad) ** 2.) /
            (np.pi * clust_rad ** 2)) - 1
        # If the number of field regions defined is larger than the maximum
        # allowed, use the maximum.
        if f_regions > f_regs_max:
            print ('  WARNING: Number of FR defined (%d) larger\n' +
            '  than the maximum allowed (%d). Using max number.') % (f_regions,
            f_regs_max)
            f_regions = f_regs_max

    # Get list that contains the spiral as a list of x,y coordinates (also
    # stored as lists) starting from the initial bin [0, 0].
    spiral = gs()

    # Obtain cluster region.
    # Calculate number of bins such that their combined area is
    # approximately (l*r)^2. See that: num_bins_area * width_bins[0]^2 =
    # (l*r)^2.
    num_bins_area = int(round(((length * clust_rad / bin_width) ** 2), 0))
    # Store 'big' cluster region as the slightly larger area around the
    # cluster's radius. Used for plotting.
    sp_indx = 0
    cl_reg_big, sp_indx = spiral_region(histo, h_manual, stars_in,
                                            stars_out, x_c_b, y_c_b, spiral,
                                            num_bins_area, sp_indx)

    # Obtain field regions.
    # This list holds all the field regions.
    field_regions = []
    if not flag_area_stronger and f_regions > 0:

        # This ensures that the decontamination algorithm uses CMD's
        # of areas equal to the cluster area for the field regions since
        # only stars inside the cluster's radius are used to obtain
        # the cluster's CMD.
        num_bins_area = int(np.pi * ((clust_rad / bin_width) ** 2))

        for f_reg in range(f_regions):
            f_region, sp_indx = spiral_region(histo, h_manual, stars_in,
                                              stars_out, x_c_b, y_c_b, spiral,
                                              num_bins_area, sp_indx)
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

        # If after removing the empty regions no regions are left, raise the
        # flag.
        if not(field_regions):
            print ('  WARNING: no field regions left after removal of those\n' +
            '  with less than 4 stars.')
            flag_area_stronger = True
    else:
        print ('  WARNING: number of field regions is 0.\n' +
        '  No field region defined.')
        flag_area_stronger = True

    return flag_area_stronger, cl_reg_big, field_regions
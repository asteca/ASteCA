# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 10:06:31 2013

@author: gabriel
"""

from get_spiral import spiral as gs
from get_histo_manual import manual_histo as mh
import numpy as np
import get_in_params as g


def spiral_index(spiral, sp_indx, histo, x_c_b, y_c_b, num_bins_area):
    '''
    Take the fixed x,y coordinates for a squared spiral of bins (spiral)
    centered at [0,0] and an index (sp_indx) that points to a given
    coordinate (where 0 means centger coords: [0,0]).

    Center the spiral at the coordinates (x_c_b, y_c_b) for a 2D histogram
    (histo) and loop through the spiralt storing the coordinates of each
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
     we use 'h_manual' wich is a list that already contains the stars that
     fall inside each of the bins in the 2D histogram along with their
     relevant data (ID, x, y, T1, etc..)
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


def get_regions(hist_lst, cent_bin, clust_rad, stars_out):
    '''
    Define empty region around the cluster via a spiral centered on it
    and of area a bit larger than that defined by the cluster's radius.

    Define as many field regions as requested (or permitted) around this empty
    region.
    '''

    hist_2d, xedges, yedges, bin_width = hist_lst
    x_c_b, y_c_b = cent_bin

    # Get area as total number of bins in 2D hist times the area of each bin.
    total_area = len(hist_2d[0]) * len(hist_2d) * (bin_width ** 2)
    cl_area = np.pi * clust_rad ** 2

    # Begin with a length factor of 2.5 and decrease it to 2.0 if the maximu
    # number of field regions that can be defined is smaller than 1.
    for l_factor in np.arange(2.5, 1.99, -0.1):
        # l_factor: Length of the side of the square that contains the cluster.

        # Increase length if the radius is comparable to the bin width used,
        # to make sure the region covers the entire cluster.
        length = l_factor if clust_rad > 2 * bin_width else (2. * l_factor)
        # Area of the square around the cluster area.
        sq_area = (length * clust_rad) ** 2.
        # Maximum number of field regions possible.
        f_regs_max = int((total_area - sq_area) / cl_area)

        if f_regs_max > 0:
            break

    # If the remaining area in the frame after substracting the cluster region
    # is smaller than the cluster's area, this means that the cluster is either
    # too large or the frame too small and no field region of equal area than
    # that of the cluster can be obtained.
    # Raise a flag.
    flag_area_stronger = False
    if f_regs_max < 1:
        print ("  WARNING: cluster region is too large or frame\n"
        "  is too small. No field regions available.")
        flag_area_stronger = True
    else:
        # If the number of field regions defined is larger than the maximum
        # allowed, use the maximum.
        if g.fr_number == 'max':
            f_regions = f_regs_max
            print 'Using max number of field regions ({}).'.format(f_regions)
        elif g.fr_number > f_regs_max:
            f_regions = f_regs_max
            print ("  WARNING: Number of FR defined ({}) larger than\n"
            "  the maximum allowed ({}). "
            "Using max number.").format(g.fr_number, f_regs_max)
        elif g.fr_number < 0:
            f_regions = f_regs_max
            print ("  WARNING: Number of FR defined ({}) is less than\n"
            "  zero. Using max number ({}).").format(g.fr_number, f_regs_max)
        else:
            print ("Using defined number of field regions ({}).".format(
                g.fr_number))
            f_regions = g.fr_number

    # Get list that contains the spiral as a list of x,y coordinates (also
    # stored as lists) starting from the initial bin [0, 0].
    spiral = gs()

    # Obtain empty area around cluster region.
    # Calculate number of bins such that their combined area is
    # approximately (l*r)^2. See that: num_bins_area * width_bins[0]^2 =
    # (l*r)^2.
    num_bins_area = int(round(((length * clust_rad / bin_width) ** 2), 0))
    # Obtain index of spiral bin where field regions should begin to be formed.
    # no_lst is not important here.
    sp_indx, no_lst = spiral_index(spiral, 0, hist_2d, x_c_b, y_c_b,
        num_bins_area)

    # Obtain field regions.
    # This list holds all the field regions.
    field_regions = []
    if not flag_area_stronger:

        # Obtain filled 2D histogram for the field with star's values attached
        # to each bin.
        h_manual = mh(stars_out, xedges, yedges)

        # This ensures that the decontamination algorithm uses CMD's
        # of areas equal to the cluster area for the field regions since
        # only stars inside the cluster's radius are used to obtain
        # the cluster's CMD.
        num_bins_area = int(np.pi * ((clust_rad / bin_width) ** 2))

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

        # If after removing the empty regions no regions are left, raise the
        # flag.
        if not(field_regions):
            print ('  WARNING: no field regions left after removal of those\n' +
            '  with less than 4 stars.')
            flag_area_stronger = True
        else:
            print "Field regions obtained."

    return flag_area_stronger, field_regions
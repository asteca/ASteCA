
import numpy as np


def main(clp, x, y, **kwargs):
    """
    Obtain the cluster's area. If the cluster is cropped, calculate the
    correct area by counting the bins that compose it's region.

    x, y are coordinates form the incomplete data set, i.e.: all stars prior to
    photometric removal.
    """
    x_max, x_min = max(x), min(x)
    y_max, y_min = max(y), min(y)

    # Used to normalize the LF of the entire frame.
    frame_area = (x_max - x_min) * (y_max - y_min)

    # Check if a portion of the cluster's region falls outside the frame.
    if clp['kde_cent'][0] + clp['clust_rad'] <= x_max and \
        clp['kde_cent'][0] - clp['clust_rad'] >= x_min and \
        clp['kde_cent'][1] + clp['clust_rad'] <= y_max and \
            clp['kde_cent'][1] - clp['clust_rad'] >= y_min:

        # Cluster's area.
        cl_area = np.pi * clp['clust_rad'] ** 2
        # Fraction of cluster's area present in frame.
        frac_cl_area = 1.

    else:
        # Index of point in RDP closer to the calculated cluster radius.
        sq_indx = int(
            round(((clp['clust_rad'] - (clp['bin_width'] / 2.)) /
                   clp['bin_width']) + 1))

        # Length of bin diagonal.
        bin_diag = np.sqrt(2 * clp['bin_width'] ** 2)

        sum_bins_in_rad = 0.
        # For each square ring until reaching the limit imposed by the
        # cluster's radius.
        for sq_ring in clp['square_rings'][:sq_indx]:
            # For each bin within this square ring.
            for bin_coords in sq_ring:
                # Coordinates of bin with center in (0, 0)
                x_b, y_b = bin_coords
                # Sign of each coordinate.
                x_s = x_b / abs(float(x_b)) if x_b != 0 else 1.
                y_s = y_b / abs(float(y_b)) if y_b != 0 else 1.
                # Coordinates of farthest corner of bin.
                x_cor, y_cor = x_b + x_s / 2., y_b + y_s / 2.
                # Coordinates of corner of bin with center in (0., 0.)
                x_c = (x_cor * clp['bin_width']) + (clp['bin_width'] / 2.)
                y_c = (y_cor * clp['bin_width']) + (clp['bin_width'] / 2.)
                # Distance of center to corner of bin.
                bin_dist = np.sqrt(x_c ** 2 + y_c ** 2)
                if bin_dist - clp['clust_rad'] > bin_diag:
                    # The entire bin is outside of the cluster radius range.
                    pass
                elif 0. < bin_dist - clp['clust_rad'] <= bin_diag:
                    # Add a portion of the bin to the total sum of bins.
                    sum_bins_in_rad += min(
                        1., (bin_dist - clp['clust_rad']) / clp['bin_width'])
                else:
                    # Add entire bin.
                    sum_bins_in_rad += 1.

        # Cluster's area.
        cl_area = (sum_bins_in_rad * clp['bin_width'] ** 2)
        # Fraction of cluster's area present in frame.
        frac_cl_area = cl_area / (np.pi * clp['clust_rad'] ** 2)

        print("  WARNING: only a portion of the cluster\n  is present "
              "in the observed frame ({:.2f})".format(frac_cl_area))

    # Normalization of frame to cluster area, used by LF plotting.
    frame_norm = frame_area / cl_area

    clp['cl_area'], clp['frac_cl_area'], clp['frame_norm'] =\
        cl_area, frac_cl_area, frame_norm
    print("Area of cluster obtained")

    return clp


import numpy as np


def main(kde_center, clust_rad, x_data, y_data, square_rings, bw):
    """
    Obtain the cluster's area. If the cluster is cropped, calculate the
    correct area by counting the bins that compose it's region.'
    """

    x_max, x_min = max(x_data), min(x_data)
    y_max, y_min = max(y_data), min(y_data)
    # Check if a portion of the cluster's region falls outside the frame.
    if kde_center[0] + clust_rad <= x_max and \
        kde_center[0] - clust_rad >= x_min and \
        kde_center[1] + clust_rad <= y_max and \
            kde_center[1] - clust_rad >= y_min:

        # Cluster's area.
        cl_area = np.pi * clust_rad ** 2
        # Fraction of cluster's area present in frame.
        frac_cl_area = 1.

    else:
        # Index of point in RDP closer to the calculated cluster radius.
        sq_indx = int(round(((clust_rad - (bw / 2.)) / bw) + 1))

        sum_bins_in_rad = 0.
        # For each square ring until reaching the limit imposed by the
        # cluster's radius.
        for sq_ring in square_rings[:sq_indx]:
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
                x_c, y_c = (x_cor * bw) + (bw / 2.), (y_cor * bw) + (bw / 2.)
                # Distance of center to corner of bin.
                bin_dist = np.sqrt(x_c ** 2 + y_c ** 2)
                # Length of bin diagonal.
                bin_diag = np.sqrt(2 * bw ** 2)
                if bin_dist - clust_rad > bin_diag:
                    # The entire bin is outside of the cluster radius range.
                    pass
                elif 0. < bin_dist - clust_rad <= bin_diag:
                    # Add a portion of the bin to the total sum of bins.
                    sum_bins_in_rad += min(1., (bin_dist - clust_rad) / bw)
                else:
                    # Add entire bin.
                    sum_bins_in_rad += 1.

        # Cluster's area.
        cl_area = (sum_bins_in_rad * bw ** 2)
        # Fraction of cluster's area present in frame.
        frac_cl_area = cl_area / (np.pi * clust_rad ** 2)

    print "Area of cluster obtained."

    return cl_area, frac_cl_area

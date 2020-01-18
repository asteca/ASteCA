
import numpy as np
from ..aux_funcs import circFrac


def main(clp, x, y, **kwargs):
    """
    Obtain the cluster's area. If the cluster is cropped, calculate the
    correct area via Monte Carlo.

    x, y are coordinates form the incomplete data set, i.e.: all stars prior to
    photometric removal.
    """
    x_max, x_min = max(x), min(x)
    y_max, y_min = max(y), min(y)

    # Used to normalize the LF of the entire frame.
    frame_area = (x_max - x_min) * (y_max - y_min)

    # Cluster's area.
    cl_area = np.pi * clp['clust_rad'] ** 2
    # Fraction of cluster's area present in frame.
    frac_cl_area = 1.

    # Check if a portion of the cluster's region falls outside the frame.
    if clp['kde_cent'][0] + clp['clust_rad'] <= x_max and \
        clp['kde_cent'][0] - clp['clust_rad'] >= x_min and \
        clp['kde_cent'][1] + clp['clust_rad'] <= y_max and \
            clp['kde_cent'][1] - clp['clust_rad'] >= y_min:
        # Nothing to do here.
        pass

    else:
        frac_cl_area = circFrac(
            clp['kde_cent'], clp['clust_rad'], x_min, x_max, y_min, y_max,
            clp['N_MC'], clp['rr'], clp['cos_t'], clp['sin_t'])
        cl_area = cl_area * frac_cl_area

        print("  WARNING: only a portion of the cluster\n  is present "
              "in the observed frame ({:.2f})".format(frac_cl_area))

    # Normalization of frame to cluster area, used by LF plotting.
    frame_norm = frame_area / cl_area

    clp['cl_area'], clp['frac_cl_area'], clp['frame_norm'] =\
        cl_area, frac_cl_area, frame_norm
    print("Area of cluster obtained")

    return clp

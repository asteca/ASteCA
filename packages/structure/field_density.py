
import numpy as np
from ..out import prep_plots


def main(clp, coords, fdens_method, **kwargs):
    """
    Get field density level of stars through an iterative process. Start with
    the complete set of radial density points and obtain its median and
    standard deviation. Then reject the point located the farthest beyond
    the 1 sigma range around the median and obtain the new median. Repeat
    the process until no points are left beyond the 1 sigma level.
    """

    # Copy list.
    reduced_rd = list(clp['rdp_points'])

    if fdens_method == 'auto':

        stable_cond = False
        while stable_cond is False:

            # Obtain median and standard deviation.
            median, sigma = np.median(reduced_rd), np.std(reduced_rd)

            # Check if at least one element in the list is beyond the 1 sigma
            # level.
            rm_elem = False
            dist_r = -1.
            for indx, elem in enumerate(reduced_rd):
                dist = abs(elem - median)

                if dist > sigma and dist > dist_r:
                    # Update distance removal value.
                    dist_r = dist
                    # Store index of element.
                    rm_index = indx
                    # Raise flag.
                    rm_elem = True

            if rm_elem is True:
                # Remove element from list and iterate again.
                del reduced_rd[rm_index]
            else:
                stable_cond = True

            field_dens = median

    else:
        field_dens = float(fdens_method)

    coord = prep_plots.coord_syst(coords)[0]
    print("Field density ({:.1E} stars/{c}^2)".format(field_dens, c=coord))
    clp['field_dens'] = field_dens

    return clp


import numpy as np
from scipy import spatial
from ..aux_funcs import circFrac
from ..out import prep_plots


def main(clp, cld_i, coords, NN_dd, fdens_method, **kwargs):
    """
    Get field density level of stars through.
    """

    # DEPRECATED 11/19
    # # Copy list.
    # reduced_rd = list(clp['rdp_points'])
    # field_dens = iterativeRDP(reduced_rd)

    clp['xy_dens'], clp['NN_dist'], clp['fr_dist'], clp['fr_dens'] = distDens(
        cld_i['x'], cld_i['y'], clp['kde_cent'], NN_dd)

    clp['fdens_min_d'], clp['fdens_lst'], clp['fdens_std_lst'],\
        clp['field_dens_d'], clp['field_dens'], clp['field_dens_std'] = kNNRDP(
        clp['fr_dist'], clp['fr_dens'], fdens_method)

    print("Field density ({:.1E} stars/{c}^2)".format(
        clp['field_dens'], c=prep_plots.coord_syst(coords)[0]))

    return clp


def distDens(x, y, kde_cent, NN_dd):
    """
    Filter stars in frame to reject those too close to the frame's borders.
    For the remaining stars, obtain their radius (distance to farthest
    neighbor), distance to the given center, and density value.
    """

    # HARDCODED outer frame limits
    lb, rt = 1., 99.
    # Filter out stars too close to the frame's borders.
    xl, xh = np.percentile(x, (lb, rt))
    yl, yh = np.percentile(y, (lb, rt))
    msk_in_frame = (x >= xl) & (x <= xh) & (y >= yl) & (y <= yh)
    x, y = x[msk_in_frame], y[msk_in_frame]
    # Frame limits
    x0, x1, y0, y1 = min(x), max(x), min(y), max(y)

    # Find NN_dd nearest neighbors.
    xy_dens = np.array((x, y)).T
    tree = spatial.cKDTree(xy_dens)
    inx = tree.query(xy_dens, k=NN_dd + 1)
    # Keep the distance to the most distant neighbor, ie: the radius.
    NN_dist = inx[0][:, NN_dd]

    # For stars close to the border frames, these areas will be overestimated.
    # That's why we filter some of them out above and use Monte Carlo here for
    # the rest of them.
    areas = np.pi * NN_dist**2
    for i, rad in enumerate(NN_dist):
        fr_area = 1.
        # If the area of the star lies outside of the frame.
        if (x[i] - rad < x0) or (x[i] + rad > x1) or (y[i] - rad < y0) or\
                (y[i] + rad > y1):
            fr_area = circFrac((x[i], y[i]), rad, x0, x1, y0, y1)
        areas[i] *= fr_area

    # Per star densities.
    dens = NN_dd / areas

    # Distances of (almost) all stars to the center of the cluster.
    dist = spatial.distance.cdist([kde_cent], xy_dens)[0]
    ds = np.argsort(dist)
    xy_dens, NN_dist, dens, dist = xy_dens[ds], NN_dist[ds], dens[ds], dist[ds]

    return xy_dens, NN_dist, dist, dens


def kNNRDP(dist, dens, fdens_method):
    """
    Use the previously obtained densities and distances to estimate the
    field density, as the mean of the X percent of the stars with the largest
    distances to the cluster's center.
    """

    fdens_min_d, fdens_lst, fdens_std_lst = [], [], []

    # HARDCODED: distance to center range.
    # dens_in_p, N_in_p = [], 0
    for perc in range(10, 99, 5):
        dist_min = np.percentile(dist, perc)
        msk_dens = dist > dist_min

        # N_in_p += msk_dens.sum()
        # dens_in_p += list(dens[msk_dens])

        # HARDCODED minimum number of stars in percentile.
        # # Only use values estimated with a reasonable number of stars.
        # if N_in_p > 100:
        fdens_min_d.append(dist_min)
        fdens_lst.append(np.median(dens[msk_dens]))
        fdens_std_lst.append(np.std(dens[msk_dens]))
        # # Reset
        # dens_in_p, N_in_p = [], 0

    if fdens_method == 'min':
        # Use the minimum value.
        i = np.argmin(fdens_lst)
        field_dens_d, field_dens, field_dens_std = fdens_min_d[i],\
            fdens_lst[i], fdens_std_lst[i]
    elif fdens_method == 'last':
        # Use the last value in the series.
        field_dens_d, field_dens, field_dens_std = fdens_min_d[-1],\
            fdens_lst[-1], fdens_std_lst[-1]
    elif fdens_method == 'iter':
        field_dens_d = iterativeRDP(fdens_lst)
        idx = np.argmin(abs(field_dens_d - np.array(fdens_lst)))
        field_dens, field_dens_std = fdens_lst[idx], fdens_std_lst[idx]
    elif fdens_method[-1] == '%':
        # Use the x% percentile.
        idx = int(len(fdens_min_d) * float(fdens_method[:-1]) / 100.)
        field_dens_d, field_dens, field_dens_std = fdens_min_d[idx],\
            fdens_lst[idx], fdens_std_lst[idx]
    else:
        # The 'nan' value identifies this as a manual value.
        field_dens_d, field_dens, field_dens_std = np.nan,\
            float(fdens_method), 0.

    return fdens_min_d, fdens_lst, fdens_std_lst, field_dens_d, field_dens,\
        field_dens_std


def iterativeRDP(fdens_lst):
    """
    Iterative process:

    1. Start with the complete set of radial density points
    2. Obtain its median and standard deviation.
    3. Reject the point located the farthest away from the 1 sigma range around
       the median
    4. Obtain new median and standard deviation.
    5. Repeat the process until no points are left beyond the 1 sigma level.
    6. Return the field density as the median value.
    """
    stable_cond = False
    while stable_cond is False:

        # Obtain median and standard deviation.
        median, sigma = np.median(fdens_lst), np.std(fdens_lst)

        # Check if at least one element in the list is beyond the 1 sigma
        # level.
        rm_elem = False
        dist_r = -1.
        for indx, elem in enumerate(fdens_lst):
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
            del fdens_lst[rm_index]
        else:
            stable_cond = True

        field_dens = median

    return field_dens

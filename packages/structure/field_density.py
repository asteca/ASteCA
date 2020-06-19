
import numpy as np
from scipy import spatial
from ..aux_funcs import circFrac
from ..out import prep_plots


def main(clp, cld_i, coords, NN_dd, fdens_method, **kwargs):
    """
    Get field density level of frame.
    """

    # Parameters that are obtained here and used in the functions below
    # and other structural analysis functions later on.
    clp['xy_filtered'], clp['xy_cent_dist'], clp['N_MC'], clp['rand_01_MC'],\
        clp['cos_t'], clp['sin_t'] = fixedParams(cld_i['x'], cld_i['y'], **clp)

    clp['NN_dist'], clp['fr_dens'] = distDens(NN_dd, **clp)

    clp['fdens_min_d'], clp['fdens_lst'], clp['fdens_std_lst'],\
        clp['field_dens_d'], clp['field_dens'], clp['field_dens_std'] = kNNRDP(
        clp['xy_cent_dist'], clp['fr_dens'], fdens_method)

    print("Field density ({:.1E} stars/{c}^2)".format(
        clp['field_dens'], c=prep_plots.coord_syst(coords)[0]))

    return clp


def fixedParams(x, y, kde_cent, lb=1., rt=99., N_MC=100000, **kwargs):
    """
    Estimate here parameters that need to be obtained only once.

    * HARDCODED:

    lb, rt: lower bottom / right top
       Frame limits in percentile values.
    N_MC: number of Monte Carlo points to use when estimating circular areas.
    """

    # Filter out stars that are too close to the frame's borders.
    xl, xh = np.percentile(x, (lb, rt))
    yl, yh = np.percentile(y, (lb, rt))
    msk_in_frame = (x >= xl) & (x <= xh) & (y >= yl) & (y <= yh)
    x, y = x[msk_in_frame], y[msk_in_frame]
    xy_filtered = np.array((x, y)).T

    # Distances of (almost) all stars to the center of the cluster.
    xy_cent_dist = spatial.distance.cdist([kde_cent], xy_filtered)[0]

    # For stars close to the border frames, the areas will be overestimated.
    # That's why we filter some of them out above and use Monte Carlo here for
    # the rest of them.
    # Obtain these values here so they are not estimated every time the MC
    # in circFrac() is called.
    rand_01_MC = np.sqrt(np.random.uniform(0., 1., N_MC))
    theta = np.random.uniform(0., 1., N_MC) * 2 * np.pi
    cos_t, sin_t = np.cos(theta), np.sin(theta)

    return xy_filtered, xy_cent_dist, N_MC, rand_01_MC, cos_t, sin_t


def distDens(NN_dd, xy_filtered, N_MC, rand_01_MC, cos_t, sin_t, **kwargs):
    """
    Obtain the NN densities and radius for each star in the frame.
    """
    # Frame limits
    x0, x1 = min(xy_filtered.T[0]), max(xy_filtered.T[0])
    y0, y1 = min(xy_filtered.T[1]), max(xy_filtered.T[1])

    # Find NN_dd nearest neighbors.
    tree = spatial.cKDTree(xy_filtered)
    inx = tree.query(xy_filtered, k=NN_dd + 1)
    # Keep the distance to the most distant neighbor, ie: the radius.
    NN_dist = inx[0][:, NN_dd]

    x, y = xy_filtered.T
    areas = np.pi * NN_dist**2
    for i, rad in enumerate(NN_dist):
        fr_area = 1.
        # If the area of the star lies outside of the frame.
        if (x[i] - rad < x0) or (x[i] + rad > x1) or (y[i] - rad < y0) or\
                (y[i] + rad > y1):
            # Use Monte Carlo to estimate its area
            fr_area = circFrac(
                (x[i], y[i]), rad, x0, x1, y0, y1, N_MC, rand_01_MC, cos_t,
                sin_t)
        areas[i] *= fr_area

    # Per star densities.
    fr_dens = NN_dd / areas

    return NN_dist, fr_dens


def kNNRDP(xy_cent_dist, dens, fdens_method):
    """
    Use the previously obtained densities and distances to estimate the
    field density.
    """

    # Obtain field density estimates using all stars beyond increasingly large
    # radius values (percentages of the distances to the pre-defined center)
    fdens_min_d, fdens_lst, fdens_std_lst = [], [], []
    # HARDCODED: distance to center range.
    # dens_in_p, N_in_p = [], 0
    for perc in range(10, 99, 5):
        dist_min = np.percentile(xy_cent_dist, perc)
        msk_dens = xy_cent_dist > dist_min

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
        field_dens = iterativeRDP(fdens_lst)
        idx = np.argmin(abs(field_dens - np.array(fdens_lst)))
        field_dens_d, field_dens_std = fdens_min_d[idx], fdens_std_lst[idx]
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
    fdens_copy = fdens_lst.copy()

    stable_cond = False
    while stable_cond is False:

        # Obtain median and standard deviation.
        median, sigma = np.median(fdens_copy), np.std(fdens_copy)

        # Check if at least one element in the list is beyond the 1 sigma
        # level.
        rm_elem = False
        dist_r = -1.
        for indx, elem in enumerate(fdens_copy):
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
            del fdens_copy[rm_index]
        else:
            stable_cond = True

        field_dens = median

    return field_dens

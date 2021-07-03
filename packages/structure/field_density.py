
import numpy as np
from scipy import spatial
from ..aux_funcs import monteCarloPars, circFrac


def main(clp, cld_i, fdens_method, **kwargs):
    """
    Get field density level of frame.
    """
    print("Estimating the field density")

    # Parameters that are obtained here and used in the functions below
    # and other structural analysis functions later on.
    clp['xy_filtered'], clp['xy_cent_dist'] = fixedParams(
        cld_i['x'], cld_i['y'], **clp)

    clp['NN_dd'], clp['NN_dist'], clp['fr_dens'] = distDens(**clp)

    clp['fdens_min_d'], clp['fdens_lst'], clp['fdens_std_lst'],\
        clp['field_dens_d'], clp['field_dens'], clp['field_dens_std'] = kNNRDP(
        clp['xy_cent_dist'], clp['fr_dens'], fdens_method)

    print("Field density ({:.3E} stars/deg^2)".format(clp['field_dens']))

    return clp


def fixedParams(x, y, kde_cent, lb=1., rt=99., **kwargs):
    """
    Estimate here parameters that need to be obtained only once.

    * HARDCODED:

    lb, rt: lower bottom / right top
       Frame limits in percentile values.
    """

    # Filter out stars that are too close to the frame's borders.
    xl, xh = np.percentile(x, (lb, rt))
    yl, yh = np.percentile(y, (lb, rt))
    msk_in_frame = (x >= xl) & (x <= xh) & (y >= yl) & (y <= yh)
    x, y = x[msk_in_frame], y[msk_in_frame]
    xy_filtered = np.array((x, y)).T

    # Distances of (almost) all stars to the center of the cluster.
    xy_cent_dist = spatial.distance.cdist([kde_cent], xy_filtered)[0]

    return xy_filtered, xy_cent_dist


def distDens(bw_list, xy_filtered, **kwargs):
    """
    Obtain the NN densities and radius for each star in the frame.

    The 'NN_dd' parameter is estimated using the bandwidth employed in the
    KDE analysis. We calculate the number of neighbors within a radius equal to
    the bandwidth, and take the median of the number of neighbors as NN_dd.
    """

    # Frame limits
    x0, x1 = min(xy_filtered.T[0]), max(xy_filtered.T[0])
    y0, y1 = min(xy_filtered.T[1]), max(xy_filtered.T[1])

    # Find NN_dd nearest neighbors.
    tree = spatial.cKDTree(xy_filtered)
    # Estimate NN_dd
    N_in_vol = tree.query_ball_point(xy_filtered, r=bw_list[1])
    NN_dd = int(np.median([len(_) for _ in N_in_vol]))

    inx = tree.query(xy_filtered, k=NN_dd + 1)
    # Keep the distance to the most distant neighbor, ie: the radius.
    NN_dist = inx[0][:, NN_dd]

    # Parameters for the Monte Carlo function
    rand_01_MC, cos_t, sin_t = monteCarloPars()

    x, y = xy_filtered.T
    areas = np.pi * NN_dist**2
    for i, rad in enumerate(NN_dist):
        fr_area = 1.
        # If the area of the star lies outside of the frame.
        if (x[i] - rad < x0) or (x[i] + rad > x1) or (y[i] - rad < y0) or\
                (y[i] + rad > y1):
            # Use Monte Carlo to estimate its area
            fr_area = circFrac(
                (x[i], y[i]), rad, x0, x1, y0, y1, rand_01_MC, cos_t,
                sin_t)
        areas[i] *= fr_area

    # Per star densities.
    fr_dens = NN_dd / areas

    return NN_dd, NN_dist, fr_dens


def kNNRDP(xy_cent_dist, fr_dens, fdens_method, pmin=1, pmax=99, Nrings=100):
    """
    Use the previously obtained densities and distances to estimate the
    field density.
    """
    dmin, dmax = np.percentile(xy_cent_dist, (pmin, pmax))
    rad_range = np.linspace(dmin, dmax, Nrings)

    # # Used for testing the King profile fit
    # dmin, dmax = np.percentile(xy_cent_dist, (.1, 25))
    # rad_range0 = np.linspace(dmin, dmax, 50)
    # dmin, dmax = np.percentile(xy_cent_dist, (25, 50))
    # rad_range1 = np.linspace(dmin, dmax, 25)
    # dmin, dmax = np.percentile(xy_cent_dist, (50, 99))
    # rad_range2 = np.linspace(dmin, dmax, 25)
    # rad_range = np.array(list(rad_range0) + list(rad_range1[1:])  + list(rad_range2[1:]))

    # Obtain field density estimates using circular rings of
    # increasingly large radius (percentages of the distances to the
    # pre-defined center)
    fdens_min_d, fdens_lst, fdens_std_lst = [], [], []
    rad_old = 0.
    for rad in rad_range:

        # Define ring with a minimum of 5 stars
        msk_dens = (xy_cent_dist > rad_old) & (xy_cent_dist <= rad)
        if msk_dens.sum() < 5:
            continue

        fdens_min_d.append((rad + rad_old) * .5)
        fdens_lst.append(np.median(fr_dens[msk_dens]))
        fdens_std_lst.append(np.std(fr_dens[msk_dens]))

        rad_old = rad

    if fdens_method == 'a':
        field_dens = iterativeRDP(fdens_lst)
    else:
        field_dens = float(fdens_method)
    idx = np.argmin(abs(field_dens - np.array(fdens_lst)))
    field_dens_d, field_dens_std = fdens_min_d[idx], fdens_std_lst[idx]

    # # Used for testing the King profile fit
    # from astropy.io import ascii
    # ascii.write(np.array([fdens_min_d, fdens_lst, fdens_std_lst]).T, "KP_test.dat")

    return fdens_min_d, fdens_lst, fdens_std_lst, field_dens_d, field_dens,\
        field_dens_std


def iterativeRDP(fdens_lst):
    """
    Iterative process:

    1. Start with the complete set of field density values
    2. Obtain its median and standard deviation.
    3. Reject the point located the farthest away from the 1 sigma range around
       the median
    4. Obtain new median and standard deviation.
    5. Repeat the process until no points are left beyond the 1 sigma level.
    6. Return the field density as the median value.
    """

    # field_dens_old = np.inf
    # for i in range(len(fdens_lst)):
    #     field_dens = np.median(fdens_lst[i + 1:])
    #     if field_dens > field_dens_old:
    #         break
    #     field_dens_old = field_dens

    fdens_copy = fdens_lst.copy()

    stable_cond = False
    while stable_cond is False:

        # Obtain median and standard deviation.
        median, sigma = np.median(fdens_copy), np.std(fdens_copy)

        # Check if at least one element in the list is beyond the 1 sigma
        # level.
        rm_elem = False
        dist_r = -1.
        for indx, fd in enumerate(fdens_copy):
            dist = abs(fd - median)

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


import numpy as np
from scipy import spatial
import scipy.integrate as integrate
from ..aux_funcs import monteCarloPars, circFrac


def main(clp, cld, fdens_method, **kwargs):
    """
    Get field density level of frame.
    """
    print("Estimating the field density")

    # Parameters that are obtained here and used in the functions below
    # and other structural analysis functions later on.
    clp['xy_filtered'], clp['xy_cent_dist'] = fixedParams(
        cld['x'], cld['y'], **clp)

    clp['pts_dens'] = distDens(clp['bw_list'], clp['xy_filtered'])

    clp['fdens_min_d'], clp['fdens_lst'], clp['fdens_std_lst'] = kNNRDP(
        clp['xy_cent_dist'], clp['pts_dens'])

    if fdens_method == 'a':
        field_dens, field_dens_std = integFieldDens(
            clp['fdens_min_d'], clp['fdens_lst'])
    else:
        field_dens, field_dens_std = float(fdens_method), np.nan
    clp['field_dens'], clp['field_dens_std'] = field_dens, field_dens_std

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


def distDens(bw_list, xy_filtered, Npts_max=150000, Nbins=100, NN_dd_max=200):
    """
    Obtain the NN densities and radius for each star in the frame.

    The 'NN_dd' parameter is estimated using the bandwidth employed in the
    KDE analysis. We calculate the number of neighbors within a radius equal to
    the bandwidth, and take the median of the number of neighbors as NN_dd.

    HARDCODED to avoid eating up all the memory
    Npts_max: maximum number of stars before using the coarse method to
    estimate per-star densities
    Nbins: number of bins used for the 2D histogram in case of too many stars
    NN_dd_max: maximum for 'NN_dd'
    """
    # Use this approach for large fields
    if xy_filtered.shape[0] > Npts_max:

        points = xy_filtered.T

        def histoDens(Nbins):
            H, xedges, yedges = np.histogram2d(*points, Nbins)
            bin_area = (xedges[1] - xedges[0]) * (yedges[1] - yedges[0])

            # Find indexes of points within edges
            xi = np.digitize(points[0], xedges, right=True)
            yi = np.digitize(points[1], yedges, right=True)

            # Handle border cases
            xi_e = np.clip(xi - 1, a_min=0, a_max=np.inf).astype(int)
            yi_e = np.clip(yi - 1, a_min=0, a_max=np.inf).astype(int)

            pts_dens = np.zeros(points.shape[-1])
            for i in range(points.shape[-1]):
                pts_dens[i] = H[xi_e[i], yi_e[i]]

            return pts_dens / bin_area

        pts_dens = histoDens(Nbins)
        msk = pts_dens == 0.
        pts_dens[msk] = np.median(pts_dens)

        return pts_dens

    else:
        # Frame limits
        x0, x1 = min(xy_filtered.T[0]), max(xy_filtered.T[0])
        y0, y1 = min(xy_filtered.T[1]), max(xy_filtered.T[1])

        # Find NN_dd nearest neighbors.
        tree = spatial.cKDTree(xy_filtered)
        # Estimate NN_dd
        N_in_vol = tree.query_ball_point(xy_filtered, r=bw_list[1])
        NN_dd = int(np.median([len(_) for _ in N_in_vol]))
        NN_dd = min(NN_dd, NN_dd_max)

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
        pts_dens = NN_dd / areas

        return pts_dens


def kNNRDP(xy_cent_dist, pts_dens, pmin=1, pmax=99, Nrings=100):
    """
    Use the previously obtained densities and distances to estimate the
    field density.
    """
    dmin, dmax = np.percentile(xy_cent_dist, (pmin, pmax))
    rad_range = np.linspace(dmin, dmax, Nrings)

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
        fdens_lst.append(np.median(pts_dens[msk_dens]))
        # Only for plotting
        fdens_std_lst.append(np.std(pts_dens[msk_dens]))

        rad_old = rad

    return fdens_min_d, fdens_lst, fdens_std_lst


def integFieldDens(fdens_min_d, fdens_lst, step=5, perc=.95):
    """
    Obtain the field density integrating the (fdens_min_d, fdens_lst) arrays
    with in steps. When the ratio of integrals to the first integral reaches
    a fixed percentage (perc), select that as the cut radius. Beyond this
    value all densities are assumed to represent the field density.

    The final value and STDDEV are the median and STDDEV of the post cut radius
    density values.
    """

    xx = np.linspace(0., 1., 1000)
    xp = np.linspace(0, 1, len(fdens_min_d))
    interp_lst = []
    for lst in (fdens_min_d, fdens_lst):
        interp_lst.append(np.interp(xx, xp, lst))
    rads, dens = interp_lst

    def curveInt(N1, N2, x, y):
        xx, yy = x[N1:N2], y[N1:N2]
        area_T = integrate.simpson(yy, xx)
        area = area_T
        return area

    i_old = 0
    rads_v, vals = [], []
    for i in range(step, len(rads), step):
        integ = curveInt(i_old, i, rads, dens)
        vals.append(integ)
        rads_v.append(np.mean(rads[i_old:i]))
        i_old = i - int(step * .5)
    rads_v = np.array(rads_v)
    areas = 1 - np.array(vals) / vals[0]

    membs_ratio = areas / areas.max()
    xp = np.linspace(0, 1, len(areas))
    rads_v_interp = np.interp(xx, xp, rads_v)
    membs_ratio_interp = np.interp(xx, xp, membs_ratio)

    idx = int(len(membs_ratio_interp) * .5)
    for i, v in enumerate(membs_ratio_interp):
        if v > perc:
            idx = i
            break
    f_rad = rads_v_interp[idx]
    j = np.argmin(abs(f_rad - rads))
    field_dens, field_dens_std = np.median(dens[j:]), np.std(dens[j:])

    return field_dens, field_dens_std


# DEPRECATED Feb 2022
# def iterativeRDP(fdens_lst):
#     """
#     Iterative process:

#     1. Start with the complete set of field density values
#     2. Obtain its median and standard deviation.
#     3. Reject the point located the farthest away from the 1 sigma range around
#        the median
#     4. Obtain new median and standard deviation.
#     5. Repeat the process until no points are left beyond the 1 sigma level.
#     6. Return the field density as the median value.
#     """

#     # field_dens_old = np.inf
#     # for i in range(len(fdens_lst)):
#     #     field_dens = np.median(fdens_lst[i + 1:])
#     #     if field_dens > field_dens_old:
#     #         break
#     #     field_dens_old = field_dens

#     fdens_copy = fdens_lst.copy()

#     stable_cond = False
#     while stable_cond is False:

#         # Obtain median and standard deviation.
#         median, sigma = np.median(fdens_copy), np.std(fdens_copy)

#         # Check if at least one element in the list is beyond the 1 sigma
#         # level.
#         rm_elem = False
#         dist_r = -1.
#         for indx, fd in enumerate(fdens_copy):
#             dist = abs(fd - median)

#             if dist > sigma and dist > dist_r:
#                 # Update distance removal value.
#                 dist_r = dist
#                 # Store index of element.
#                 rm_index = indx
#                 # Raise flag.
#                 rm_elem = True

#         if rm_elem is True:
#             # Remove element from list and iterate again.
#             del fdens_copy[rm_index]
#         else:
#             stable_cond = True

#         field_dens = median

#     return field_dens

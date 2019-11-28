
import numpy as np
from scipy import spatial
import warnings
from .contamination_index import CIfunc
from ..aux_funcs import circFrac
from ..out import prep_plots
from .. import update_progress


def main(cld_i, clp, coords, rad_manual, nsteps_rad, **kwargs):
    """
    Estimate the radius through the optimization of the #cluster-members vs
    #field-stars values. Assign the uncertainty through a bootstrap process
    based on the field density's uncertainty.
    """

    # RDP. For plotting purposes only, except the 'rdp_radii' values which are
    # used below.
    clp['rdp_radii'], clp['rdp_points'], clp['rdp_stddev'] = kNNRDP(
        clp['fr_dist'], clp['fr_dens'])

    radii, n_in_cl_reg, areas = radsAreas(
        clp['rdp_radii'], cld_i['x'], cld_i['y'], clp['kde_cent'])

    # Obtain the optimal radius and arrays for plotting.
    clp['clust_rad'], clp['rad_rads'], clp['rad_N_membs'],\
        clp['rad_N_field'], clp['rad_CI'] = optimalRadius(
            clp['field_dens'], radii, n_in_cl_reg, areas)

    coord = prep_plots.coord_syst(coords)[0]
    if rad_manual == 'n':
        print("Estimating radius")
        # Use bootstrap to estimate the uncertainty.
        clp['e_rad'] = radError(
            clp['field_dens'], clp['field_dens_std'], radii, n_in_cl_reg,
            areas, nsteps_rad)
        print("Radius found: {:g} {}".format(clp['clust_rad'], coord))

    elif rad_manual != 'n':
        # Update radius, assign zero error.
        clp['clust_rad'], clp['e_rad'] = rad_manual, 0.
        print("Manual radius set: {:g} {}".format(clp['clust_rad'], coord))

    return clp


def radsAreas(rdp_radii, x, y, kde_cent):
    """
    Helper function. Determine the radii values, and estimate the number of
    stars and areas for each value.
    """

    # Frame limits
    x0, x1, y0, y1 = min(x), max(x), min(y), max(y)
    dx0, dx1 = abs(kde_cent[0] - x0), abs(kde_cent[0] - x1)
    dy0, dy1 = abs(kde_cent[1] - y0), abs(kde_cent[1] - y1)
    dxy = min(dx0, dx1, dy0, dy1)
    # Distances to center
    dist = spatial.distance.cdist([kde_cent], np.array([x, y]).T)[0]

    # HARDCODED this '200' value gives reasonable results with reasonable
    # performance
    radii = np.linspace(rdp_radii[0], rdp_radii[-1], 200)

    # Areas and #stars for all rad values.
    areas, n_in_cl_reg = np.pi * radii**2, []
    for i, rad in enumerate(radii):
        fr_area = 1.
        if rad > dxy:
            fr_area = circFrac((kde_cent), rad, x0, x1, y0, y1)
        areas[i] *= fr_area
        # Stars within radius.
        n_in_cl_reg.append((dist < rad).sum())

    return radii, n_in_cl_reg, areas


def radiiParams(radii, n_in_cl_reg, areas, fd):
    """
    Given an array if radii, estimate the number of members, the number of
    field stars, and the contamination index within it.
    """
    rads, N_membs, N_field, cont_index, break_counter = [], [], [], [], 0
    for i, rad in enumerate(radii):

        if n_in_cl_reg[i] == 0:
            continue

        CI, n_memb, n_fl = CIfunc(n_in_cl_reg[i], fd, areas[i])

        # Break check
        if n_memb <= 0. and i > .5 * len(radii):
            break_counter += 1
            if break_counter > 3:
                break

        rads.append(rad)
        N_field.append(max(0., n_fl))
        N_membs.append(max(0., n_memb))
        cont_index.append(max(0, CI))

    return rads, N_membs, N_field, cont_index


def optimalRadius(field_dens, radii, n_in_cl_reg, areas):
    """
    Estimate the optimal radius as the value that maximizes the (normalized)
    number of members stars versus the number of field stars. The rationale
    is that we want the maximum possible of member stars, ans the minimum
    possible of field stars within the region.
    """

    rads, N_membs, N_field, CI = radiiParams(
        radii, n_in_cl_reg, areas, field_dens)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # Normalizing separately is important.
        N_membs = np.array(N_membs) / max(N_membs)
        N_field = np.array(N_field) / max(N_field)
        clust_rad = rads[np.argmax(N_membs - N_field)]

    return clust_rad, rads, N_membs, N_field, CI


def radError(field_dens, field_dens_std, radii, n_in_cl_reg, areas, N_btstrp):
    """
    Bootstrap the field density to estimate the radius uncertainty.
    """

    steps = int(.1 * N_btstrp)
    all_rads = np.empty(N_btstrp)
    for i, fd in enumerate(
            np.random.normal(field_dens, field_dens_std, N_btstrp)):
        all_rads[i] = optimalRadius(fd, radii, n_in_cl_reg, areas)[0]

        if i + 1 == N_btstrp:
            pass
        elif (i + 1) % steps:
            continue
        update_progress.updt(N_btstrp, i + 1)

    # From the bootstrap samples we estimate the mean standard error.
    return np.std(all_rads)


def kNNRDP(fr_dist, fr_dens):
    """
    Use the kNN's per-star densities. Average these values for several circular
    rings to obtain the RDP.
    """

    # HARDCODED: remove the more conflicting last ~10% of radii values.
    radii = np.linspace(0., fr_dist.max(), 55)[:-5]

    rdp_radii, rdp_NN, rdp_stddev = [], [], []
    for l, h in zip(*[radii[:-1], radii[1:]]):
        # Stars within this ring.
        msk_in = (fr_dist >= l) & (fr_dist < h)
        if sum(msk_in) > 0:
            rdp_radii.append(.5 * (l + h))
            rdp_NN.append(np.median(fr_dens[msk_in]))
            rdp_stddev.append(np.std(fr_dens[msk_in]))

    if not rdp_radii:
        raise ValueError("ERROR: RDP is empty. Check the center coordinates")

    return rdp_radii, rdp_NN, rdp_stddev


# DEPRECATED 11/19
# def radius_algor(clp, coord, radius_method):
#     """
#     Obtain the value for the cluster's radius by counting the number of points
#     that fall within a given interval of the field density or lower. If this
#     number is equal to a minimum fixed number of points 'n_left', then assign
#     the radius as the point closest to the field density value among those
#     first n_left points counting from the first one that fell below the
#     (field dens + delta limit) range.
#     Iterate increasing the interval around the field density and finally
#     average all the radius values found for each interval.
#     """

#     radii, rdp_points, bin_width, field_dens = clp['rdp_radii'],\
#         clp['rdp_points'], clp['bin_width'], clp['field_dens']
#     # Find maximum density value and assume this is the central density.
#     # Do not use previous values.
#     max_dens_ind = np.argmax(rdp_points)
#     rdp_points_m, radii_m = rdp_points[max_dens_ind:], radii[max_dens_ind:]

#     # Interpolate extra RDP points between those calculated.
#     N = 1000
#     t, xp = np.linspace(0, 1, N), np.linspace(0, 1, len(rdp_points_m))
#     # Store interpolated RDP and radii values.
#     rdp_points_c = np.interp(t, xp, rdp_points_m)
#     radii_c = np.interp(t, xp, radii_m)

#     # Assign a value to the number of points that should be found below
#     # the delta values around the field density to attain the 'stabilized'
#     # condition.
#     # Fix to X% of the total number of interpolated points in the RDP.
#     lmh = {'low': 0.05, 'mid': 0.1, 'high': 0.2}
#     n_left = int(lmh[radius_method] * N)

#     # Difference between max RDP density value and the field density value.
#     delta_total = (max(rdp_points_c) - field_dens)

#     # If the difference between the max density value and the field density is
#     # less than 3 times the value of the field density, raise a flag.
#     flag_delta_total = False
#     if delta_total < 3 * field_dens:
#         flag_delta_total = True

#     # Initialize index_rad value.
#     rad_found = []
#     for k, delta_percentage in enumerate(np.arange(0.2, 0.1, -0.01)):
#         # The 'k' param relaxes the condition that requires a certain number of
#         # points to be located below the 'delta_field + field_dens' value
#         # to establish that the RDP has stabilized.
#         # The smaller the value of 'delta_field', the fewer the number of
#         # consecutive points that are required.

#         # delta_field = % of difference between max density value and field
#         # density.
#         delta_field = delta_percentage * delta_total

#         # Initialize density values counter for points that fall inside the
#         # range determined by the delta value around the field density.
#         in_delta_val, index_rad_i, dens_dist = 0, 0, np.inf

#         # Iterate through all values of star density in the RDP.
#         for index, dens in enumerate(rdp_points_c):

#             # Condition to iterate until at least in_delta_val *consecutive*
#             # points below the (delta + field density) value are found.
#             if in_delta_val < (n_left - (4 * k)):

#                 # If the density value is closer than 'delta_field' to the
#                 # field density value or lower --> add it.
#                 if dens <= delta_field + field_dens:
#                     # Augment value of counter.
#                     in_delta_val += 1
#                     # Store radius value closer to the field density.
#                     if abs(dens - field_dens) < dens_dist:
#                         dens_dist = abs(dens - field_dens)
#                         index_rad_i = index
#                 # If the RDP point is outside the (delta + field density) range
#                 # reset all values. I.e.: the points should be located below
#                 # the 'delta_field + field_dens' *consecutively*.
#                 else:
#                     # Reset.
#                     in_delta_val, index_rad_i, dens_dist = 0, 0, np.inf

#             # If enough RDP points have been found within the field density
#             # range, store the radius value closer to the field density value
#             # and break out of the for loop.
#             else:
#                 rad_found.append(radii_c[index_rad_i])
#                 break

#     # Raise a flag if only two radius values were found under the deltas.
#     flag_delta, flag_not_stable = False, False
#     if len(rad_found) < 3:
#         flag_delta = True

#     # If at least one radius value was found.
#     if rad_found:
#         # Use the median to avoid outliers.
#         clust_rad = np.median(rad_found)

#         # Obtain error as the 1 sigma confidence interval (68.27%).
#         _std = np.std(rad_found)
#         # Catch warning if stats fails to obtain confidence interval.
#         with warnings.catch_warnings():
#             warnings.simplefilter("ignore")
#             conf_int = stats.norm.interval(0.6827, loc=clust_rad, scale=_std /
#                                            np.sqrt(len(rad_found)))
#         # If stats returns a confidence interval with a NaN, discard it.
#         if np.any(np.isnan(conf_int)):
#             conf_dif = 0.
#         else:
#             conf_dif = abs(clust_rad - max(conf_int))
#         e_rad = max(conf_dif, bin_width)
#         # Prevent too small radius by fixing the minimum value to the second
#         # RDP point.
#         if clust_rad < radii[1]:
#             clust_rad = radii[1]
#     else:
#         flag_not_stable = True
#         # No radius value found. Assign radius value as the middle element
#         # in the radii list.
#         clust_rad, e_rad = radii_c[int(len(radii_c) / 2.)], 0.
#         print('  WARNING: no radius found, setting value to: {:g} {}'.format(
#             clust_rad, coord))

#     return clust_rad, e_rad, flag_delta_total, flag_not_stable, flag_delta


import numpy as np
from scipy import spatial
from .contamination_index import CIfunc
from ..aux_funcs import circFrac
from ..out import prep_plots
from .. import update_progress


def main(cld_i, clp, coords, rad_manual, nsteps_rad, NN_dd, **kwargs):
    """
    Estimate the radius through the optimization of the #cluster-members vs
    #field-stars values. Assign the uncertainty through a bootstrap process
    based on the field density's uncertainty.
    """

    # Parameters used internally only.
    rad_radii, rad_areas, st_dists_cent = rdpAreasDists(cld_i, **clp)

    # Obtain the optimal radius and arrays for plotting.
    clp['clust_rad'], clp['rad_rads'], clp['rad_N_membs'], clp['rad_N_field'],\
        clp['rad_CI'] = optimalRadius(
            rad_radii, rad_areas, st_dists_cent, **clp)

    coord = prep_plots.coord_syst(coords)[0]
    if rad_manual == 'n':
        print("Estimating radius")

        if nsteps_rad > 2:
            # Use bootstrap to estimate the uncertainty.
            clp['e_rad'] = radError(
                rad_radii, rad_areas, st_dists_cent, nsteps_rad, **clp)
            print("Radius found: {:g} {}".format(clp['clust_rad'], coord))
        else:
            clp['e_rad'] = np.array([np.nan, np.nan])
            print("Radius found (no bootstrap): {:g} {}".format(
                clp['clust_rad'], coord))

    elif rad_manual != 'n':
        # Update radius, assign zero error.
        clp['clust_rad'], clp['e_rad'] = rad_manual, np.array([np.nan, np.nan])
        print("Manual radius set: {:g} {}".format(clp['clust_rad'], coord))

    return clp


def rdpAreasDists(
        cld_i, kde_cent, rdp_radii, N_MC, rand_01_MC, cos_t, sin_t, **kwargs):
    """
    The areas for each radius value in 'rdp_radii' are obtained here once.
    We also calculate here the distance of each star to the defined center.
    """

    # Proper array of coordinates.
    xy = np.array([cld_i['x'], cld_i['y']]).T

    # Frame limits
    x0, x1 = min(xy.T[0]), max(xy.T[0])
    y0, y1 = min(xy.T[1]), max(xy.T[1])

    dx0, dx1 = abs(kde_cent[0] - x0), abs(kde_cent[0] - x1)
    dy0, dy1 = abs(kde_cent[1] - y0), abs(kde_cent[1] - y1)
    dxy = min(dx0, dx1, dy0, dy1)

    # We use the 'rdp_radii' here since it is already processed to contain
    # reasonable starting and finishing values.
    # Defining this array here gives the radius finding function more
    # resolution,almost independently of the number of RDP rings used.
    rad_radii = np.linspace(rdp_radii[0], rdp_radii[-1], 100)

    # Areas associated to the radii defined in 'rdp_radii'.
    rad_areas = np.pi * np.array(rad_radii)**2
    for i, rad in enumerate(rad_radii):
        fr_area = 1.
        if rad > dxy:
            fr_area = circFrac(
                (kde_cent), rad, x0, x1, y0, y1, N_MC, rand_01_MC, cos_t,
                sin_t)
        rad_areas[i] *= fr_area

    # Distances of all stars to the center of the cluster.
    st_dists_cent = spatial.distance.cdist([kde_cent], xy)[0]

    return rad_radii, rad_areas, st_dists_cent


def optimalRadius(
    rad_radii, rad_areas, st_dists_cent, field_dens, bttstrp_flag=False,
        **kwargs):
    """
    Estimate the optimal radius as the value that maximizes the (normalized)
    number of members stars versus the number of field stars. The rationale
    is that we want the maximum possible of member stars, ans the minimum
    possible of field stars within the region.
    """

    data, break_counter = [], 0
    for i, rad in enumerate(rad_radii):

        # Stars within radius.
        n_in_cl_reg = (st_dists_cent <= rad).sum()
        if n_in_cl_reg == 0:
            continue

        # Tried (21/01/20) the variable field density method, but it failed.
        # It is too fragile and ends up selecting very large radius values.
        # Field density
        # (dist > rad).sum(): #stars outside the cluster region
        # (area_frame - rdp_areas[i]): are outside the cluster region.
        # field_dens = (N_tot - n_in_cl_reg) / (area_tot - rad_areas[i])

        if bttstrp_flag:
            n_fl = field_dens * rad_areas[i]
            n_memb = n_in_cl_reg - n_fl
            CI = np.nan
        else:
            CI, n_memb, n_fl = CIfunc(n_in_cl_reg, field_dens, rad_areas[i])

        # Break check
        if (n_memb <= 0. or CI > 1.1) and i > .5 * len(rad_radii):
            break_counter += 1
            if break_counter > 3:
                break
        data.append([rad, n_memb, n_fl, CI])

    # rads, N_membs, N_field, CI
    data = np.clip(data, a_min=0., a_max=None).T

    # Normalizing separately is important. Otherwise it selects low radii
    # values.
    N_membs = data[1] / data[1].max()
    N_field = data[2] / data[2].max()
    idx = np.argmax(N_membs - N_field)
    clust_rad = data[0][idx]

    if bttstrp_flag:
        return clust_rad

    return clust_rad, data[0], N_membs, N_field, data[3]


def radError(
        rad_radii, rad_areas, st_dists_cent, N_btstrp, field_dens, **kwargs):
    """
    Bootstrap the distances to estimate the radius uncertainty.
    """

    steps = int(.1 * N_btstrp)
    all_rads = np.empty(N_btstrp)
    for i in range(N_btstrp):

        # Sample distances with replacement.
        st_dists_btstrp = np.random.choice(st_dists_cent, st_dists_cent.size)

        # Obtain the optimal radius and arrays for plotting.
        all_rads[i] = optimalRadius(
            rad_radii, rad_areas, st_dists_btstrp, field_dens, True)

        if i + 1 == N_btstrp:
            pass
        elif (i + 1) % steps:
            continue
        update_progress.updt(N_btstrp, i + 1)

    return np.percentile(all_rads, (16, 84))


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

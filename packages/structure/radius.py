
import numpy as np
import warnings
from ..out import prep_plots

from scipy.spatial.distance import cdist
from scipy.signal import savgol_filter

import scipy.stats as stats
from scipy.stats import median_abs_deviation as MAD
from .contamination_index import CIfunc
from ..aux_funcs import monteCarloPars, circFrac


def main(cld_i, clp, coords, rad_manual, clust_rad_mode, **kwargs):
    """
    Estimate the radius through the optimization of the #cluster-members vs
    #field-stars values. Assign the uncertainty through a bootstrap process
    based on the field density's uncertainty.

    Nboot: number of bootstrap runs
    """
    print("Estimating the radius")

    clp['all_rads'], clp['e_rad'] = np.array([]), np.array([np.nan, np.nan])

    coord = prep_plots.coord_syst(coords)[0]
    if rad_manual == 'n':
        # Use bootstrap to estimate the uncertainty.
        rad_radii, rad_areas, N_in_cl_rad, N_in_ring = rdpAreasDists(
            cld_i['x'], cld_i['y'], clp['kde_cent'], clp['xy_cent_dist'])
        clp['clust_rad'] = optimalRadius(
            rad_radii, rad_areas, N_in_cl_rad, N_in_ring, clp['field_dens'])

        # if not np.isnan(clp['field_dens_std']):
        #     clp['all_rads'], clp['e_rad'] = radError(
        #         rad_radii, rad_areas, N_in_cl_rad, N_in_ring,
        #         clp['field_dens'], clp['field_dens_std'])

        # clp['clust_rad'], clp['all_rads'], clp['e_rad'] = radBootstrp(
        #     clp['xy_cent_dist'], clp['fr_dens'], clp['field_dens'],
        #     clp['field_dens_std'], Nboot)

        if clust_rad_mode == 'auto':
            print("Radius found: {:g} {}".format(clp['clust_rad'], coord))

        elif clust_rad_mode == 'max':
            clp['clust_rad'] = maxRadius(
                cld_i['x'], cld_i['y'], clp['kde_cent'])
            print("Large radius selected: {:g} {}".format(
                clp['clust_rad'], coord))

    elif rad_manual != 'n':
        clp['clust_rad'] = rad_manual
        print("Manual radius set: {:g} {}".format(clp['clust_rad'], coord))

    return clp


def optimalRadius(rad_radii, rad_areas, N_in_cl_rad, N_in_ring, field_dens):
    """
    Estimate the optimal radius as the value that...
    """
    # data, break_counter, N_membs = [], 0, []
    # for i, rad in enumerate(rad_radii):

    #     # # Stars within radius.
    #     # n_in_cl_reg = (xy_cent_dist <= rad).sum()
    #     # if n_in_cl_reg == 0:
    #     #     continue

    #     # Ring's area
    #     # area_ring = rad_areas[i] - area_old
    #     # area_old = rad_areas[i]

    #     # # Stars within ring
    #     # N_in_ring = n_in_cl_reg - N_in_old
    #     # N_in_old = n_in_cl_reg

    #     # # Estimated number of field stars in ring
    #     # n_fl_r = field_dens * area_ring

    #     # # Estimated number of members in ring
    #     # n_memb_r = N_in_r - n_fl_r
    #     # # Cumulative number of members
    #     # n_memb = n_memb_r + n_memb_old
    #     # n_memb_old = n_memb

    #     n_memb = N_in_cl_rad[i] - field_dens * rad_areas[i]

    #     # Break check
    #     if n_memb <= 0. and i > .5 * len(rad_radii):
    #         break_counter += 1
    #         if break_counter > 3:
    #             break
    #     N_membs.append(n_memb)

    N_membs = N_in_cl_rad - field_dens * rad_areas

    # rads, N_membs, N_field, CI
    # data = np.clip(data, a_min=0., a_max=None).T

    # # Normalizing separately is important. Otherwise it selects low radii
    # # values.
    # membs_max = data[1].max()
    # with warnings.catch_warnings():
    #     warnings.simplefilter("ignore")
    #     N_membs = data[1] / membs_max
    # N_field = data[2] / data[2].max()

    # # Smooth the curve. Tried smoothing before the normalization but it
    # # increases the uncertainty region enormously.
    # # Catch error in savgol_filter() when there is a 'nan' or 'inf'
    # if np.isnan(N_membs).any() or np.isinf(N_membs).any():
    #     N_membs_smooth = N_membs
    # else:
    #     # ws: window size, pol: polynomial order
    #     ws, pol = int(len(N_membs) / 5.), 3
    #     # must be odd
    #     ws = ws + 1 if ws % 2 == 0 else ws
    #     N_membs_smooth = savgol_filter(N_membs, ws, pol)

    # Optimal radius selection
    # idx = np.argmax(N_membs_smooth - N_field)
    idx = np.argmax(N_membs / N_in_ring)
    clust_rad = rad_radii[idx]

    # # if plotflag:
    # print(field_dens, idx, clust_rad)
    # import matplotlib.pyplot as plt
    # plt.subplot(121)
    # plt.scatter(rad_radii, N_membs / N_in_ring)
    # plt.subplot(122)
    # plt.scatter(rad_radii, N_membs, c='g')
    # plt.scatter(rad_radii, N_in_ring, c='r')
    # plt.show()

    # return data[0], data[3], membs_max, N_membs_smooth, N_field, clust_rad
    return clust_rad


def rdpAreasDists(
    x, y, kde_cent, xy_cent_dist, pmin=3, pmax=90, Nrads=300, Nmax=50000,
        N_MC=1000000):
    """
    The areas for each radius value in 'rad_radii' are obtained here once.
    We also calculate here the distance of each star to the defined center.

    HARDCODED
    pmin, pmax: minimum and maximum percentiles used to define the radii range
    Nrads: number of values used to generate the 'rad_radii' array.
    Nmax: maximum number of stars used in the process
    N_MC: points in the Monte Carlo area estimation. Use 1e6 for stability.
    """

    rand_01_MC, cos_t, sin_t = monteCarloPars(N_MC)

    # Define the radii values
    dmin, dmax = np.percentile(xy_cent_dist, (pmin, pmax))
    all_rads = np.linspace(dmin, dmax, Nrads)

    # Frame limits
    x0, x1 = min(x), max(x)
    y0, y1 = min(y), max(y)

    # Estimate the minimum distance from the center of the cluster to any
    # border of the frame
    dx0, dx1 = abs(kde_cent[0] - x0), abs(kde_cent[0] - x1)
    dy0, dy1 = abs(kde_cent[1] - y0), abs(kde_cent[1] - y1)
    dxy = min(dx0, dx1, dy0, dy1)

    # Areas associated to the radii defined in 'all_rads'.
    rad_areas, rad_radii, N_in_ring, N_in_cl_rad =\
        np.pi * np.array(all_rads)**2, [], [], []
    N_in_old = 0
    for i, rad in enumerate(all_rads):

        # Stars within radius.
        n_in_cl_reg = (xy_cent_dist <= rad).sum()
        if n_in_cl_reg == 0:
            continue

        rad_radii.append(rad)
        # Stars within radius
        N_in_cl_rad.append(n_in_cl_reg)

        fr_area = 1.
        if rad > dxy:
            fr_area = circFrac(
                (kde_cent), rad, x0, x1, y0, y1, rand_01_MC, cos_t, sin_t)
        rad_areas[i] *= fr_area

        # Stars within ring
        N_in_ring.append(n_in_cl_reg - N_in_old)
        N_in_old = n_in_cl_reg

    xx = np.linspace(0., 1., 1000)
    xp = np.linspace(0, 1, len(rad_radii))
    interp_lst = []
    for lst in (rad_radii, rad_areas, N_in_cl_rad, N_in_ring):
        interp_lst.append(np.interp(xx, xp, lst))
    rad_radii, rad_areas, N_in_cl_rad, N_in_ring = interp_lst

    return rad_radii, rad_areas, N_in_cl_rad, N_in_ring


def maxRadius(x, y, kde_cent):
    """
    Estimate a large radius that encloses the entire frame
    """
    xmin, xmax, ymin, ymax = min(x), max(x), min(y), max(y)
    pts = np.array([(xmin, ymin), (xmin, ymax), (xmax, ymin), (xmax, ymax)])
    # Fast euclidean distance: https://stackoverflow.com/a/47775357/1391441
    a_min_b = kde_cent - pts
    dist = np.sqrt(np.einsum('ij,ij->i', a_min_b, a_min_b))
    # Use the maximum distance as the radius
    clust_rad = max(dist)

    return clust_rad


# # DEPRECATED Nov 2020; RE-IMPLEMENTED April 2021
# def radError(
#     rad_radii, rad_areas, N_in_cl_rad, N_in_ring, field_dens, field_dens_std,
#         Nboot=1000):
#     """
#     Bootstrap the distances to estimate the radius uncertainty.
#     """
#     field_dens_s = np.random.normal(field_dens, field_dens_std, Nboot)

#     # mu, sigma = field_dens, field_dens_std
#     # lower, upper = mu - sigma, mu + sigma
#     # X = stats.truncnorm(
#     #     (lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma)
#     # field_dens_s = X.rvs(Nboot)

#     field_dens_s[field_dens_s <= 0.] = field_dens

#     # steps = int(.1 * Nboot)
#     all_rads = []
#     for field_dens_s_i in field_dens_s:

#         # Obtain the optimal radius and arrays for plotting.
#         rad = optimalRadius(
#             rad_radii, rad_areas, N_in_cl_rad, N_in_ring, field_dens_s_i)
#         if not np.isnan(rad):
#             all_rads.append(rad)

#         # if i + 1 == Nboot:
#         #     pass
#         # elif (i + 1) % steps:
#         #     continue
#         # update_progress.updt(Nboot, i + 1)

#     all_rads = np.array(all_rads)
#     msk = all_rads < np.median(all_rads) + MAD(all_rads) * 5
#     if msk.sum() > 10:
#         all_rads = all_rads[msk]

#     return all_rads, np.percentile(all_rads, (16, 84))


# def radBootstrp(xy_cent_dist, fr_dens, field_dens, field_dens_std, Nboot):
#     """
#     """
#     steps = int(.1 * Nboot)

#     count_max = np.random.randint(1, 4, Nboot)
#     Nrings = np.random.randint(50, 100, Nboot)
#     pmin = np.random.randint(2, 5, Nboot)
#     pmax = np.random.randint(90, 99, Nboot)

#     all_rads = []
#     for i in range(Nboot):
#         all_rads.append(
#             densRad(xy_cent_dist, fr_dens, field_dens, field_dens_std,
#                     Nrings[i], count_max[i], pmin[i], pmax[i]))
#         if i + 1 == Nboot:
#             pass
#         elif (i + 1) % steps:
#             continue
#         update_progress.updt(Nboot, i + 1)

#     return np.median(all_rads), np.array(all_rads),\
#         np.percentile(all_rads, (16, 84))


# def densRad(
#     xy_cent_dist, fr_dens, field_dens, field_dens_std, Nrings, count_max,
#         pmin, pmax, stars_min=10):
#     """
#     Estimate the radius as the value where the density (from the center of the
#     cluster) reaches the region of the field density:
#     (field_dens - field_dens_std)
#     """
#     dmin, dmax = np.percentile(xy_cent_dist, (pmin, pmax))
#     rad_radii = np.linspace(dmin, dmax, Nrings)

#     rad_old, count = 0, 0
#     clust_rad = rad_radii[int(Nrings / 2.)]
#     for rad in rad_radii:
#         # Define ring with a minimum of 10 stars
#         msk_dens = (xy_cent_dist > rad_old) & (xy_cent_dist <= rad)
#         if msk_dens.sum() < stars_min:
#             continue

#         fdens_r = np.median(fr_dens[msk_dens])
#         if abs(fdens_r - field_dens) <= field_dens_std:
#             count += 1
#             if count >= count_max:
#                 clust_rad = rad
#                 break
#         rad_old = rad

#     return clust_rad

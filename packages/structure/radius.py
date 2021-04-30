
import numpy as np
from scipy.signal import savgol_filter
from ..out import prep_plots
from ..aux_funcs import monteCarloPars, circFrac


def main(cld_i, clp, coords, rad_method, **kwargs):
    """
    Estimate the radius through the optimization of the #cluster-members vs
    #field-stars values. Assign the uncertainty through a bootstrap process
    based on the field density's uncertainty.

    Nboot: number of bootstrap runs
    """
    print("Estimating the radius")
    clp['rad_radii'], clp['rad_areas'], clp['N_in_cl_rad'], clp['N_in_ring'] =\
        rdpAreasDists(cld_i['x'], cld_i['y'], clp['kde_cent'],
                      clp['xy_cent_dist'], clp['field_dens'])

    coord = prep_plots.coord_syst(coords)[0]
    if rad_method == 'a':

        clp['clust_rad'] = optimalRadius(
            clp['rad_radii'], clp['rad_areas'], clp['N_in_cl_rad'],
            clp['N_in_ring'], clp['field_dens'])
        print("Radius found: {:g} {}".format(clp['clust_rad'], coord))

        # if not np.isnan(clp['field_dens_std']):
        #     clp['all_rads'], clp['e_rad'] = radError(
        #         rad_radii, rad_areas, N_in_cl_rad, N_in_ring,
        #         clp['field_dens'], clp['field_dens_std'])

        # clp['clust_rad'], clp['all_rads'], clp['e_rad'] = radBootstrp(
        #     clp['xy_cent_dist'], clp['fr_dens'], clp['field_dens'],
        #     clp['field_dens_std'], Nboot)

    elif rad_method == 'max':
        clp['clust_rad'] = maxRadius(
            cld_i['x'], cld_i['y'], clp['kde_cent'])
        print("Large radius selected: {:g} {}".format(
            clp['clust_rad'], coord))

    else:
        clp['clust_rad'] = rad_method
        print("Manual radius set: {:g} {}".format(clp['clust_rad'], coord))

    return clp


def optimalRadius(rad_radii, rad_areas, N_in_cl_rad, N_in_ring, field_dens):
    """
    Estimate the optimal radius as the rad value where the 'N_membs/N_in_ring'
    ratio is maximized.
    """
    N_membs = N_in_cl_rad - field_dens * rad_areas
    membs_ratio = N_membs / N_in_ring

    if np.isnan(membs_ratio).any() or np.isinf(membs_ratio).any():
        pass
    else:
        # window size (must be odd). There are ~1000 values by default
        ws = max(3, int(len(membs_ratio) / 10.))
        ws = ws + 1 if ws % 2 == 0 else ws
        # polynomial order
        pol = 3
        membs_ratio = savgol_filter(membs_ratio, ws, pol)

    idx = np.argmax(membs_ratio)
    clust_rad = rad_radii[idx]

    return clust_rad


def rdpAreasDists(
    x, y, kde_cent, xy_cent_dist, field_dens, pmin=2, pmax=90, Nrads=300,
        N_MC=1000000, Ninterp=1000):
    """
    The areas for each radius value in 'rad_radii' are obtained here once.
    We also calculate here the distance of each star to the defined center.

    HARDCODED
    pmin, pmax: minimum and maximum percentiles used to define the radii range
    Nrads: number of values used to generate the 'rad_radii' array.
    N_MC: points in the Monte Carlo area estimation. Use 1e6 for stability.
    Ninterp: number of interpolated points in the final arrays
    """

    rand_01_MC, cos_t, sin_t = monteCarloPars(N_MC)

    # # Define the radii values
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

        # Total stars within ring
        N_in_ring.append(n_in_cl_reg - N_in_old)
        N_in_old = n_in_cl_reg

    # INterpolate extra points
    xx = np.linspace(0., 1., Ninterp)
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


# # DEPRECATED Nov 2020; RE-IMPLEMENTED April 2021 (it had almost no effect
# # on the radius uncertainty)
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

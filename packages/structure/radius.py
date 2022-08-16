
import numpy as np
# import scipy.integrate as integrate
from scipy.spatial.distance import cdist
from scipy import stats
from ..out.prep_plots import RDPCurve


def main(cld, clp, rad_method, **kwargs):
    """
    Estimate the radius through the integrals of the field density array.
    Assign the uncertainty based on the field density's uncertainty.
    """
    print("Estimating the radius")
    clp['rad_uncert'], clp['rads_fit_line'] = np.array([np.nan, np.nan]),\
        np.array([])

    if rad_method == 'a':
        clp['clust_rad'], clp['rad_uncert'], clp['rads_fit_line'] =\
            optimalRadius(
                clp['xy_filtered'], clp['xy_cent_dist'], clp['kde_cent'],
                clp['field_dens'], clp['field_dens_std'])

        print("Radius found: {:g} deg".format(clp['clust_rad']))

    elif rad_method == 'max':
        clp['clust_rad'] = maxRadius(
            cld['x'], cld['y'], clp['kde_cent'], clp['xy_cent_dist'])
        print("Max radius selected: {:g} deg".format(
            clp['clust_rad']))
    else:
        clp['clust_rad'] = rad_method
        print("Manual radius set: {:g} deg".format(clp['clust_rad']))

    return clp


def optimalRadius(
    xy_filtered, xy_cent_dist, kde_cent, field_dens, field_dens_std,
        rvalue2_min=0.9):
    """
    """
    # Obtain the RDP using the concentric rings method.
    rdp_radii, rdp_points, _ = RDPCurve(
        0, xy_filtered, xy_cent_dist, kde_cent, 0, 0)

    # Logarithm values. The first three points are removed because they tend
    # to disturb the fit
    x, y = np.log(rdp_radii[3:]), np.log(rdp_points[3:])
    xx = np.linspace(x.min(), x.max(), 100)
    N = len(rdp_radii)

    def getCLRad(fdens):
        for i in range(int(len(rdp_radii) * .5)):
            x_j, y_j = x[:N - i], y[:N - i]
            res = stats.linregress(x_j, y_j)
            idx = np.argmin(abs(
                res.intercept + res.slope * xx - np.log(field_dens + fdens)))
            clust_rad = np.exp(xx[idx])
            if res.rvalue**2 > rvalue2_min:
                break
        return clust_rad, res

    # Estimate value
    rad, res = getCLRad(0)
    cl_rads = [rad]
    rads_fit_line = np.array([
        res.intercept, res.slope, res.intercept_stderr, res.stderr])
    # Estimate STDDEV range
    for N_std in (-1, 1):
        if not np.isnan(field_dens_std):
            fdens = N_std * field_dens_std
            cl_rads.append(getCLRad(fdens)[0])
        else:
            cl_rads.append(np.nan)

    clust_rad, rmax, rmin = cl_rads
    rmin = min(rmin, clust_rad)
    rmax = max(rmax, clust_rad)
    rad_uncert = np.array([rmin, rmax])

    return clust_rad, rad_uncert, rads_fit_line


# def optimalRadius(x, y, kde_cent, field_dens, field_dens_std, xy_cent_dist):
#     """
#     """
#     max_lx = min(max(x) - kde_cent[0], kde_cent[0] - min(x))
#     max_ly = min(max(y) - kde_cent[1], kde_cent[1] - min(y))
#     max_l = min(max_lx, max_ly)
#     d_idx = np.argsort(xy_cent_dist)
#     min_l = xy_cent_dist[d_idx[25]]  # HARDCODED

#     N_tot = (xy_cent_dist < max_ly).sum()
#     N_rings = max(5, int(N_tot / 25))

#     rads = np.linspace(min_l, max_l, N_rings)
#     all_vals = []
#     r_old = 0
#     for r in rads:
#         # Stars inside this ring
#         Nr = ((xy_cent_dist >= r_old) & (xy_cent_dist < r)).sum()
#         fd = np.random.normal(field_dens, field_dens_std, 1000)
#         # Estimated members stars inside this ring
#         Nm_r = Nr - np.pi * (r**2 - r_old**2) * fd
#         Nc = np.median(Nm_r)
#         Nc_min = Nc + np.std(Nm_r)
#         Nc_max = Nc - np.std(Nm_r)
#         all_vals.append([r, Nc_min, Nc_max])
#         r_old = r
#     all_vals = np.array(all_vals).T

#     i_low, i_high = 0, len(rads)
#     for i, _ in enumerate(all_vals[1]):
#         if _ <= 0:
#             i_low = i
#             break
#     for i, _ in enumerate(all_vals[2]):
#         if _ <= 0:
#             i_high = i
#             break

#     rad_uncert = np.array([rads[i_low], rads[i_high]])
#     clust_rad = (rads[i_high] + rads[i_low]) * .5

#     return clust_rad, rad_uncert


# def optimalRadius(field_dens, field_dens_std, fdens_min_d, fdens_lst):
#     """
#     """
#     # Interpolate extra points into de radii and field density arrays.
#     xx = np.linspace(0., 1., 1000)
#     xp = np.linspace(0, 1, len(fdens_min_d))
#     interp_lst = []
#     for lst in (fdens_min_d, fdens_lst):
#         interp_lst.append(np.interp(xx, xp, lst))
#     rads, dens = interp_lst

#     clust_rad = integRad(field_dens, xx, xp, rads, dens)
#     if np.isnan(field_dens_std):
#         r_16, r_84 = np.nan, np.nan
#     else:
#         r_16 = integRad(field_dens + field_dens_std, xx, xp, rads, dens)
#         r_84 = integRad(field_dens - field_dens_std, xx, xp, rads, dens)
#     rad_uncert = np.array([r_16, r_84])

#     rads_interp, integ_interp = integRad(field_dens, xx, xp, rads, dens, True)

#     return clust_rad, rad_uncert, rads_interp, integ_interp


# def integRad(
#         field_dens, xx, xp, rads, dens, memb_flag=False, step=5, perc=.975):
#     """
#     The (fdens_min_d, fdens_lst) arrays are integrated in steps and the areas
#     normalized (removing the area for the field density). When these ratios
#     reach the 'perc' value, select that as the cluster radius.

#     The hard-coded values for 'step' and 'perc' were selected after much
#     testing.
#     """

#     dens = np.clip(dens, a_min=field_dens, a_max=np.inf)

#     def curveInt(N1, N2, x, y):
#         xx, yy = x[N1:N2], y[N1:N2]
#         area_T = integrate.simpson(yy, xx)
#         area_F = field_dens * (xx.max() - xx.min())
#         area = area_T - area_F
#         return max(1.e-6, area)

#     i_old = 0
#     rads_v, vals = [], []
#     for i in range(step, len(rads), step):
#         vals.append(curveInt(i_old, i, rads, dens))
#         rads_v.append(np.mean(rads[i_old:i]))
#         i_old = i - int(step * .5)
#     rads_v = np.array(rads_v)

#     # Avoid an error adding a small value here
#     vals[0] += 1.e-6

#     areas = 1 - np.array(vals) / vals[0]
#     areas_ratio = areas / areas.max()
#     xp = np.linspace(0, 1, len(areas_ratio))
#     rads_v_interp = np.interp(xx, xp, rads_v)
#     areas_ratio_interp = np.interp(xx, xp, areas_ratio)

#     if memb_flag is True:
#         return rads_v_interp, areas_ratio_interp

#     # Default cluster value in case the method fails
#     clust_rad = rads_v_interp[int(len(rads_v_interp) * .9)]
#     for i, v in enumerate(areas_ratio_interp):
#         if v > perc:
#             clust_rad = rads_v_interp[i]
#             break
#     return clust_rad


def maxRadius(x, y, kde_cent, xy_cent_dist):
    """
    Estimate a radius large enough for every last star
    """
    clust_rad = cdist([kde_cent], np.array((x, y)).T)[0].max()

    return clust_rad

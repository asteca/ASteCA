
import numpy as np
import scipy.integrate as integrate
from scipy.spatial.distance import cdist


def main(cld, clp, rad_method, **kwargs):
    """
    Estimate the radius through the integrals of the field density array.
    Assign the uncertainty based on the field density's uncertainty.
    """
    print("Estimating the radius")

    rad_uncert, rads_interp, integ_interp = np.array([np.nan, np.nan]),\
        np.array([]), np.array([])
    if rad_method == 'a':
        clp['clust_rad'], rad_uncert, rads_interp, integ_interp =\
            optimalRadius(
                clp['field_dens'], clp['field_dens_std'],
                np.array(clp['fdens_min_d']), np.array(clp['fdens_lst']))
        print("Radius found: {:g} deg".format(clp['clust_rad']))

    elif rad_method == 'max':
        clp['clust_rad'] = maxRadius(
            cld['x'], cld['y'], clp['kde_cent'], clp['xy_cent_dist'])
        print("Max radius selected: {:g} deg".format(
            clp['clust_rad']))

    else:
        clp['clust_rad'] = rad_method
        print("Manual radius set: {:g} deg".format(clp['clust_rad']))

    clp['rad_uncert'], clp['rads_interp'], clp['integ_interp'] =\
        rad_uncert, rads_interp, integ_interp

    return clp


def optimalRadius(field_dens, field_dens_std, fdens_min_d, fdens_lst):
    """
    """
    # Interpolate extra points into de radii and field density arrays.
    xx = np.linspace(0., 1., 1000)
    xp = np.linspace(0, 1, len(fdens_min_d))
    interp_lst = []
    for lst in (fdens_min_d, fdens_lst):
        interp_lst.append(np.interp(xx, xp, lst))
    rads, dens = interp_lst

    clust_rad = integRad(field_dens, xx, xp, rads, dens)
    if np.isnan(field_dens_std):
        r_16, r_84 = np.nan, np.nan
    else:
        r_16 = integRad(field_dens + field_dens_std, xx, xp, rads, dens)
        r_84 = integRad(field_dens - field_dens_std, xx, xp, rads, dens)
    rad_uncert = np.array([r_16, r_84])

    rads_interp, integ_interp = integRad(field_dens, xx, xp, rads, dens, True)

    return clust_rad, rad_uncert, rads_interp, integ_interp


def integRad(
        field_dens, xx, xp, rads, dens, memb_flag=False, step=5, perc=.975):
    """
    The (fdens_min_d, fdens_lst) arrays are integrated in steps and the areas
    normalized (removing the area for the field density). When these ratios
    reach the 'perc' value, select that as the cluster radius.

    The hard-coded values for 'step' and 'perc' were selected after much
    testing.
    """

    dens = np.clip(dens, a_min=field_dens, a_max=np.inf)

    def curveInt(N1, N2, x, y):
        xx, yy = x[N1:N2], y[N1:N2]
        area_T = integrate.simpson(yy, xx)
        area_F = field_dens * (xx.max() - xx.min())
        area = area_T - area_F
        return max(1.e-6, area)

    i_old = 0
    rads_v, vals = [], []
    for i in range(step, len(rads), step):
        vals.append(curveInt(i_old, i, rads, dens))
        rads_v.append(np.mean(rads[i_old:i]))
        i_old = i - int(step * .5)
    rads_v = np.array(rads_v)

    # Avoid an error adding a small value here
    vals[0] += 1.e-6

    areas = 1 - np.array(vals) / vals[0]
    areas_ratio = areas / areas.max()
    xp = np.linspace(0, 1, len(areas_ratio))
    rads_v_interp = np.interp(xx, xp, rads_v)
    areas_ratio_interp = np.interp(xx, xp, areas_ratio)

    if memb_flag is True:
        return rads_v_interp, areas_ratio_interp

    # Default cluster value in case the method fails
    clust_rad = rads_v_interp[int(len(rads_v_interp) * .9)]
    for i, v in enumerate(areas_ratio_interp):
        if v > perc:
            clust_rad = rads_v_interp[i]
            break
    return clust_rad


def maxRadius(x, y, kde_cent, xy_cent_dist):
    """
    Estimate a radius large enough for every last star
    """
    clust_rad = cdist([kde_cent], np.array((x, y)).T)[0].max()

    return clust_rad

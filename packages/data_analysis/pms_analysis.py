
import numpy as np
from scipy.spatial.distance import cdist
from scipy import stats
from scipy import spatial
from astropy.stats import sigma_clipped_stats


def main(
    clp, cld_i, coords, project, flag_make_plot, flag_PM_coord, PM_nnmax,
        PM_KDE_std, **kwargs):
    """
    """
    PM_flag, PM_flag_all, plx_pm_flag = False, False, False
    pmMP, pmRA_DE, e_pmRA_DE, pmDE, e_pmDE, mmag_pm, PM_cl_x, PM_cl_y,\
        PM_cl_z, PM_fl_x, PM_fl_y, pm_mag_fl, pmDE_fl, e_pmDE_fl, pmRA_fl_DE,\
        e_pmRA_fl_DE, pm_dist_max, pm_Plx_cl, pm_ePlx_cl, pm_Plx_fr,\
        pm_ePlx_fr, pmRA_all, pmDE_all, pmMag_all, xRA_all, yDE_all,\
        PM_kde_all, PMs_d_median = [[] for _ in range(28)]
    PM_fl_z = np.array([])
    PMs_cl_cx, PMs_cl_cy, PMs_fl_cx, PMs_fl_cy = [np.nan] * 4

    # Check PMs data.
    PM_flag = checkPMs(clp, 1)
    if PM_flag:

        print("Processing proper motions.")

        if 'C3' in flag_make_plot:
            PM_flag_all, pmRA_all, pmDE_all, pmMag_all, xRA_all, yDE_all,\
                PM_kde_all = PMsKDEAll(coords, project, cld_i, clp, PM_KDE_std)
            PMs_d_median = PMsNNAll(
                PM_flag_all, flag_PM_coord, xRA_all, yDE_all, pmRA_all,
                pmDE_all, PM_nnmax)

        # Cluster region data.
        pmMP, pmRA, e_pmRA, pmDE, e_pmDE =\
            np.array(list(zip(*clp['cl_reg_fit']))[9]),\
            np.array(list(zip(*list(zip(*clp['cl_reg_fit']))[7]))[1]),\
            np.array(list(zip(*list(zip(*clp['cl_reg_fit']))[8]))[1]),\
            np.array(list(zip(*list(zip(*clp['cl_reg_fit']))[7]))[2]),\
            np.array(list(zip(*list(zip(*clp['cl_reg_fit']))[8]))[2])
        mmag_pm = np.array(list(zip(*list(zip(*clp['cl_reg_fit']))[3]))[0])

        # Check Plx data.
        plx_pm_flag = checkPMs(clp, 0)

        if plx_pm_flag:
            pm_Plx_cl = np.array(
                list(zip(*list(zip(*clp['cl_reg_fit']))[7]))[0])
            pm_ePlx_cl = np.array(
                list(zip(*list(zip(*clp['cl_reg_fit']))[8]))[0])

        # Apply cos(delta) correction to pmRA if possible.
        if coords == 'deg':
            if project:
                DE_pm = np.array(list(zip(*clp['cl_reg_fit']))[2]) +\
                    clp['y_offset']
            else:
                DE_pm = np.array(list(zip(*clp['cl_reg_fit']))[2])
        else:
            print("  WARNING: can not apply cos(dec) factor to pmRA.")
            DE_pm = np.zeros(pmRA.size)

        # Remove nan PM values from cluster region
        msk = ~np.isnan(pmRA) & ~np.isnan(e_pmRA) & ~np.isnan(pmDE) &\
            ~np.isnan(e_pmDE)
        pmMP, pmRA, e_pmRA, pmDE, e_pmDE, DE_pm, mmag_pm =\
            pmMP[msk], pmRA[msk], e_pmRA[msk], pmDE[msk], e_pmDE[msk],\
            DE_pm[msk], mmag_pm[msk]
        if plx_pm_flag:
            pm_Plx_cl, pm_ePlx_cl = pm_Plx_cl[msk], pm_ePlx_cl[msk]

        pmRA_DE = pmRA * np.cos(np.deg2rad(DE_pm))
        # Propagate error in RA*cos(delta)
        e_pmRA_DE = np.sqrt(
            (e_pmRA * pmRA * np.sin(np.deg2rad(DE_pm)))**2 +
            (e_pmDE * np.cos(np.deg2rad(DE_pm)))**2)

        # Error weighted 2D KDE for cluster region
        PM_cl_x, PM_cl_y, PM_cl_z, PMs_cl_cx, PMs_cl_cy, _ = kde_2d(
            pmRA_DE, e_pmRA_DE, pmDE, e_pmDE)

        # PM distances to the KDE's center.
        pm_dist_max = cdist(
            np.array([
                [PM_cl_x[PMs_cl_cx][PMs_cl_cy],
                 PM_cl_y[PMs_cl_cx][PMs_cl_cy]]]), np.array([pmRA_DE, pmDE]).T)

        # Process field regions
        if not clp['flag_no_fl_regs_i']:
            pmRA_fl, e_pmRA_fl, DE_fl_pm = [], [], []
            # Field region(s) data.
            for fl_rg in clp['field_regions_i']:
                pm_mag_fl += list(zip(*list(zip(*fl_rg))[3]))[0]
                pmRA_fl += list(zip(*list(zip(*fl_rg))[7]))[1]
                e_pmRA_fl += list(zip(*list(zip(*fl_rg))[8]))[1]
                pmDE_fl += list(zip(*list(zip(*fl_rg))[7]))[2]
                e_pmDE_fl += list(zip(*list(zip(*fl_rg))[8]))[2]
                DE_fl_pm += list(list(zip(*fl_rg)))[2]
                if plx_pm_flag:
                    pm_Plx_fr += list(zip(*list(zip(*fl_rg))[7]))[0]
                    pm_ePlx_fr += list(zip(*list(zip(*fl_rg))[8]))[0]

            pm_mag_fl, pmRA_fl, e_pmRA_fl, pmDE_fl, e_pmDE_fl, DE_fl_pm,\
                pm_Plx_fr, pm_ePlx_fr = [np.asarray(_) for _ in (
                    pm_mag_fl, pmRA_fl, e_pmRA_fl, pmDE_fl, e_pmDE_fl,
                    DE_fl_pm, pm_Plx_fr, pm_ePlx_fr)]

            # Remove nan values from field region(s)
            msk = ~np.isnan(pmRA_fl) & ~np.isnan(e_pmRA_fl) &\
                ~np.isnan(pmDE_fl) & ~np.isnan(e_pmDE_fl)
            pm_mag_fl, pmRA_fl, e_pmRA_fl, pmDE_fl, e_pmDE_fl, DE_fl_pm = \
                pm_mag_fl[msk], pmRA_fl[msk], e_pmRA_fl[msk], pmDE_fl[msk],\
                e_pmDE_fl[msk], DE_fl_pm[msk]
            if plx_pm_flag:
                pm_Plx_fr, pm_ePlx_fr = pm_Plx_fr[msk], pm_ePlx_fr[msk]

            pmRA_fl_DE = pmRA_fl * np.cos(np.deg2rad(DE_fl_pm))
            # Propagate error in RA*cos(delta)
            e_pmRA_fl_DE = np.sqrt(
                (e_pmRA_fl * pmRA_fl * np.sin(np.deg2rad(DE_fl_pm)))**2 +
                (e_pmDE_fl * np.cos(np.deg2rad(DE_fl_pm)))**2)

            # Error weighted 2D KDE for field region
            PM_fl_x, PM_fl_y, PM_fl_z, PMs_fl_cx, PMs_fl_cy, _ = kde_2d(
                pmRA_fl_DE, e_pmRA_fl_DE, pmDE_fl, e_pmDE_fl)

    clp.update({
        'PM_flag': PM_flag, 'pmMP': pmMP, 'pmRA_DE': pmRA_DE,
        'e_pmRA_DE': e_pmRA_DE, 'pmDE': pmDE, 'e_pmDE': e_pmDE,
        'mmag_pm': mmag_pm, 'PM_cl_x': PM_cl_x, 'PM_cl_y': PM_cl_y,
        'PM_cl_z': PM_cl_z, 'PMs_cl_cx': PMs_cl_cx, 'PMs_cl_cy': PMs_cl_cy,
        'pmRA_fl_DE': pmRA_fl_DE, 'e_pmRA_fl_DE': e_pmRA_fl_DE,
        'pmDE_fl': pmDE_fl, 'e_pmDE_fl': e_pmDE_fl, 'pm_mag_fl': pm_mag_fl,
        'PM_fl_x': PM_fl_x, 'PM_fl_y': PM_fl_y, 'PM_fl_z': PM_fl_z,
        'PMs_fl_cx': PMs_fl_cx, 'PMs_fl_cy': PMs_fl_cy,
        'pm_dist_max': pm_dist_max, 'pm_Plx_cl': pm_Plx_cl,
        'pm_ePlx_cl': pm_ePlx_cl, 'pm_Plx_fr': pm_Plx_fr,
        'pm_ePlx_fr': pm_ePlx_fr, 'plx_pm_flag': plx_pm_flag,
        'PM_flag_all': PM_flag_all, 'PM_kde_all': PM_kde_all,
        'pmRA_all': pmRA_all, 'pmDE_all': pmDE_all,
        'PMs_d_median': PMs_d_median, 'pmMag_all': pmMag_all,
        'xRA_all': xRA_all, 'yDE_all': yDE_all})
    return clp


def checkPMs(clp, _id):
    """
    Check that Plx/PMs were defined within the cluster region.

    id=0 --> plx
    id=1 --> PMs (RA)
    """
    # Extract cluster region Plx / RA PMs data.
    data = np.array(list(zip(*list(zip(*clp['cl_reg_fit']))[7]))[_id])
    # Array with no nan values
    clrg = data[~np.isnan(data)]
    if clrg.any() and np.min(clrg) < np.max(clrg):
        return True
    else:
        return False


def PMsKDEAll(coords, project, cld_i, clp, PM_KDE_std):
    """
    Process *all* the PMs in the frame.
    """
    # Extract and mask nan values
    pmRA_all, e_pmRA_all = cld_i['kine'][1], cld_i['ek'][1]
    pmDE_all, e_pmDE_all = cld_i['kine'][2], cld_i['ek'][2]
    msk = (~np.isnan(pmRA_all)) & (~np.isnan(pmDE_all))
    pmRA_all, e_pmRA_all = pmRA_all[msk], e_pmRA_all[msk]
    pmDE_all, e_pmDE_all = pmDE_all[msk], e_pmDE_all[msk]
    xRA_all, yDE_all, pmMag_all = cld_i['x'][msk], cld_i['y'][msk],\
        cld_i['mags'][0][msk]

    # If there are values left, process.
    if pmRA_all.any():
        PM_flag_all = True
        # Apply cos(delta) correction to pmRA if possible.
        if coords == 'deg':
            if project:
                DE_pm = cld_i['y'][msk] + clp['y_offset']
            else:
                DE_pm = cld_i['y'][msk]
        else:
            DE_pm = np.zeros(pmRA_all.size)
        pmRA_all = pmRA_all * np.cos(np.deg2rad(DE_pm))

        # Error weighted 2D KDE for cluster+field region.
        PMx_kdeall, PMy_kdeall, PMz_kdeall, _, _, extent =\
            kde_2d(pmRA_all, e_pmRA_all, pmDE_all, e_pmDE_all, PM_KDE_std, 100)
        PM_kde_all = [PMx_kdeall, PMy_kdeall, PMz_kdeall, extent]
    else:
        PM_flag_all, PM_kde_all = True, []

    return PM_flag_all, pmRA_all, pmDE_all, pmMag_all, xRA_all, yDE_all,\
        PM_kde_all


def PMsNNAll(
    PM_flag_all, flag_PM_coord, xRA_all, yDE_all, pmRA_all, pmDE_all,
        nnmax):
    """
    Nearest-neighbor analysis that combines coordinates and PMs data.
    """
    if PM_flag_all:

        # Normalize data
        def norm(data):
            data_norm = []
            for arr in data:
                dmin = np.min(arr)
                arr_n = arr - dmin
                dmax = np.max(arr_n)
                arr_n /= dmax
                data_norm.append(arr_n)
            return np.array(data_norm)

        if flag_PM_coord is True:
            data = norm([xRA_all, yDE_all, pmRA_all, pmDE_all]).T
        else:
            data = norm([pmRA_all, pmDE_all]).T

        # Create the tree
        tree = spatial.cKDTree(data)
        # Find the closest nnmax-1 neighbors (first entry is the point itself)
        dists = tree.query(data, nnmax)
        # Extract distances
        nn_dist = dists[0][:, 1:]
        # Median values
        PMs_d_median = np.median(nn_dist, 1)

        return PMs_d_median


def kde_2d(xarr, xsigma, yarr, ysigma, Nstd=3, grid_dens=50):
    '''
    Take an array of x,y data with their errors, create a grid of points in x,y
    and return the 2D KDE density map.
    '''
    _, x_median, x_std = sigma_clipped_stats(xarr)
    _, y_median, y_std = sigma_clipped_stats(yarr)

    # Mask zoomed region
    x_range, y_range = 2. * x_std * Nstd, 2. * y_std * Nstd
    xyrange = .5 * max(x_range, y_range)
    xmin, xmax = x_median - xyrange, x_median + xyrange
    ymin, ymax = y_median - xyrange, y_median + xyrange
    msk = (xarr > xmin) & (xarr < xmax) & (yarr > ymin) & (yarr < ymax)
    xarr, xsigma, yarr, ysigma = xarr[msk], xsigma[msk], yarr[msk],\
        ysigma[msk]

    gd_c = complex(0, grid_dens)
    # Define grid of points in x,y where the KDE will be evaluated.
    ext = [xmin, xmax, ymin, ymax]
    x, y = np.mgrid[ext[0]:ext[1]:gd_c, ext[2]:ext[3]:gd_c]
    pos = np.vstack([x.ravel(), y.ravel()])
    values = np.vstack([xarr, yarr])

    # Replace 0 error with very LARGE value.
    np.place(xsigma, xsigma <= 0., 1000.)
    np.place(ysigma, ysigma <= 0., 1000.)
    # Inverse errors as weights
    w = 1. / (xsigma * ysigma)

    d, n = values.shape
    sf = .5 * n**(-1. / (d + 4))
    kernel = stats.gaussian_kde(values, bw_method=sf, weights=w)
    z = np.reshape(kernel(pos).T, x.shape)

    # Max value
    zmax_x, zmax_y = np.unravel_index(z.argmax(), z.shape)

    return x, y, z, zmax_x, zmax_y, [xmin, xmax, ymin, ymax]

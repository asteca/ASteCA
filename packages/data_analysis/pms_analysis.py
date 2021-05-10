
import numpy as np
from scipy import stats
from astropy.stats import sigma_clipped_stats


def main(clp, cld_i, project, flag_make_plot, **kwargs):
    """
    Assume that the pmRA is already corrected by the cosine of DE.
    """
    clreg_PMs, fregs_PMs, allfr_PMs = {}, {}, {}
    cr_KDE_PMs, fr_KDE_PMs, allr_KDE_PMs = {}, {}, {}
    # Check PMs data.
    PM_flag = checkPMs(clp)

    if PM_flag and 'C3' in flag_make_plot:
        print("Processing proper motions")

        # Extract PMs data.
        clreg_PMs, fregs_PMs, allfr_PMs = PMsData(cld_i, clp)
        cr_KDE_PMs, fr_KDE_PMs, allr_KDE_PMs = PMsKDE(
            clreg_PMs, fregs_PMs, allfr_PMs)

        # DEPRECATED 04/2021
        # # Cosine(DE) correction on the RA (if possible and needed)
        # clreg_PMs, fregs_PMs, allfr_PMs = pmRAcosDE(
        #    coords, project, clp['y_offset'], clreg_PMs, fregs_PMs, allfr_PMs)

        # DEPRECATED 04/2021
        # if 'C2' in flag_make_plot:
        #     # Plx data filtered by PMs nan values
        #     plx_pm_flag, pm_Plx_cl, pm_Plx_fr = plx_PMs_data(
        #         clp, clreg_PMs['msk'], fregs_PMs['msk'])

        # DEPRECATED 04/2021
        # # PM distances to the KDE's center.
        # xmax, ymax = cr_KDE_PMs['zmax_x'], cr_KDE_PMs['zmax_y']
        # pm_dist_max = cdist(
        #     np.array([[
        #         cr_KDE_PMs['x'][xmax][ymax],
        #         cr_KDE_PMs['y'][xmax][ymax]]]),
        #     np.array([clreg_PMs['pmRA'], clreg_PMs['pmDE']]).T)[0]

    clp.update({
        'PM_flag': PM_flag, 'clreg_PMs': clreg_PMs, 'fregs_PMs': fregs_PMs,
        'allfr_PMs': allfr_PMs, 'cr_KDE_PMs': cr_KDE_PMs,
        'fr_KDE_PMs': fr_KDE_PMs, 'allr_KDE_PMs': allr_KDE_PMs})
    return clp


def checkPMs(clp):
    """
    Check that non nan PMs are defined within the cluster region.

    id=0 --> plx
    id=1 --> PMs (RA)
    """
    # Extract cluster region Plx / RA PMs data.
    data = np.array(list(zip(*list(zip(*clp['cl_reg_fit']))[7]))[1])
    # Array with no nan values
    clrg = data[~np.isnan(data)]
    if clrg.any() and np.min(clrg) < np.max(clrg):
        return True
    else:
        return False


def PMsData(cld_i, clp):
    """
    Extract data and remove nan PM values from cluster, field regions, and the
    entire frame.
    """

    # Cluster region data.
    pmMP, pmRA, e_pmRA, pmDE, e_pmDE =\
        np.array(list(zip(*clp['cl_reg_fit']))[9]),\
        np.array(list(zip(*list(zip(*clp['cl_reg_fit']))[7]))[1]),\
        np.array(list(zip(*list(zip(*clp['cl_reg_fit']))[8]))[1]),\
        np.array(list(zip(*list(zip(*clp['cl_reg_fit']))[7]))[2]),\
        np.array(list(zip(*list(zip(*clp['cl_reg_fit']))[8]))[2])
    DE_cr = np.array(list(zip(*clp['cl_reg_fit']))[2])
    mmag_cr = np.array(list(zip(*list(zip(*clp['cl_reg_fit']))[3]))[0])
    # Remove nan values from cluster region
    msk_cr = ~np.isnan(pmRA) & ~np.isnan(e_pmRA) & ~np.isnan(pmDE) &\
        ~np.isnan(e_pmDE)
    clreg_PMs = {
        'pmRA': pmRA[msk_cr], 'epmRA': e_pmRA[msk_cr],
        'pmDE': pmDE[msk_cr], 'epmDE': e_pmDE[msk_cr],
        'DE': DE_cr[msk_cr], 'mmag': mmag_cr[msk_cr],
        'MP': pmMP[msk_cr], 'msk': msk_cr}

    # Field regions data
    fregs_PMs = {
        'pmRA': np.array([]), 'epmRA': np.array([]),
        'pmDE': np.array([]), 'epmDE': np.array([]),
        'DE': np.array([]), 'mmag': np.array([]), 'msk': np.array([])}
    if not clp['flag_no_fl_regs_i']:
        msk_fr, mmag_fr, pmRA_fr, e_pmRA_fr, pmDE_fr, e_pmDE_fr, DE_fr =\
            [[] for _ in range(7)]
        for fl_rg in clp['field_regions_i']:
            mmag_fr += list(zip(*list(zip(*fl_rg))[3]))[0]
            pmRA_fr += list(zip(*list(zip(*fl_rg))[7]))[1]
            e_pmRA_fr += list(zip(*list(zip(*fl_rg))[8]))[1]
            pmDE_fr += list(zip(*list(zip(*fl_rg))[7]))[2]
            e_pmDE_fr += list(zip(*list(zip(*fl_rg))[8]))[2]
            DE_fr += list(list(zip(*fl_rg)))[2]
        # To arrays
        mmag_fr, pmRA_fr, e_pmRA_fr, pmDE_fr, e_pmDE_fr, DE_fr,\
            = [np.asarray(_) for _ in (
                mmag_fr, pmRA_fr, e_pmRA_fr, pmDE_fr, e_pmDE_fr, DE_fr)]
        # Mask nan values in field region(s)
        msk_fr = ~np.isnan(pmRA_fr) & ~np.isnan(e_pmRA_fr) &\
            ~np.isnan(pmDE_fr) & ~np.isnan(e_pmDE_fr)
        fregs_PMs = {
            'pmRA': pmRA_fr[msk_fr], 'epmRA': e_pmRA_fr[msk_fr],
            'pmDE': pmDE_fr[msk_fr], 'epmDE': e_pmDE_fr[msk_fr],
            'DE': DE_fr[msk_fr], 'mmag': mmag_fr[msk_fr], 'msk': msk_fr}

    # Entire frame data
    pmRA_all, e_pmRA_all = cld_i['kine'][1], cld_i['ek'][1]
    pmDE_all, e_pmDE_all = cld_i['kine'][2], cld_i['ek'][2]
    # Remove nans
    msk_frame = (~np.isnan(pmRA_all)) & (~np.isnan(pmDE_all))
    pmRA_all, e_pmRA_all = pmRA_all[msk_frame], e_pmRA_all[msk_frame]
    pmDE_all, e_pmDE_all = pmDE_all[msk_frame], e_pmDE_all[msk_frame]
    xRA_all, yDE_all, mmag_all = cld_i['x'][msk_frame],\
        cld_i['y'][msk_frame], cld_i['mags'][0][msk_frame]
    allfr_PMs = {
        'pmRA': pmRA_all, 'epmRA': e_pmRA_all,
        'pmDE': pmDE_all, 'epmDE': e_pmDE_all,
        'DE': yDE_all, 'mmag': mmag_all, 'RA': xRA_all}

    return clreg_PMs, fregs_PMs, allfr_PMs


# DEPRECATED 04/2021
# def pmRAcosDE(
#     cosDE_flag, coords, project, y_offset, clreg_PMs, fregs_PMs,
#         allfr_PMs):
#     """
#     Correct pmRA by the cosine of DE.
#     """
#     def cosCorr(_dict):
#         # Shift declination data if necessary.
#         if project:
#             DE = _dict['DE'] + y_offset
#         else:
#             DE = _dict['DE']
#         # Propagate error in RA*cos(delta) first.
#         _dict['epmRA'] = np.sqrt(
#             (_dict['epmRA'] * _dict['pmRA'] * np.sin(np.deg2rad(DE)))**2
#             + (_dict['epmRA'] * np.cos(np.deg2rad(DE)))**2)
#         # Now apply cosine factor.
#         _dict['pmRA'] = _dict['pmRA'] * np.cos(np.deg2rad(DE))

#         return _dict

#     if cosDE_flag:
#         # Correction is already applied, nothing to do here.
#         pass
#     else:
#         if coords == 'deg':
#             arrs = []
#             for _ in (clreg_PMs, fregs_PMs, allfr_PMs):
#                 arrs.append(cosCorr(_))
#             clreg_PMs, fregs_PMs, allfr_PMs = arrs
#         else:
#             # pmRA data is not corrected by the cosine of DE and pixel
#             # coordinates are being used. The cos factor can not be applied.
#             print("  WARNING: can not apply cos(DE) factor to pmRA")

#     return clreg_PMs, fregs_PMs, allfr_PMs


# DEPRECATED 04/2021
# def plx_PMs_data(clp, msk_cr, msk_fr):
#     """
#     Extract Plx data associated to the values filtered by PMs.
#     """
#     # Cluster region
#     pm_Plx_cl = np.array(
#         list(zip(*list(zip(*clp['cl_reg_fit']))[7]))[0])[msk_cr]
#     if pm_Plx_cl.any():
#         plx_pm_flag = True

#     # Field regions
#     pm_Plx_fr = []
#     if not clp['flag_no_fl_regs_i']:
#         for fl_rg in clp['field_regions_i']:
#             pm_Plx_fr += list(zip(*list(zip(*fl_rg))[7]))[0]
#         pm_Plx_fr = np.array(pm_Plx_fr)[msk_fr]

#     return plx_pm_flag, pm_Plx_cl, pm_Plx_fr


def PMsKDE(clreg_PMs, fregs_PMs, allfr_PMs):
    """
    Error weighted 2D KDE for cluster, field regions, and all frame.
    """

    # Cluster region
    cr_KDE_PMs = kde_2d(
        clreg_PMs['pmRA'], clreg_PMs['epmRA'], clreg_PMs['pmDE'],
        clreg_PMs['epmDE'])

    # Field regions
    fr_KDE_PMs = {}
    if fregs_PMs['pmRA'].any():
        fr_KDE_PMs = kde_2d(
            fregs_PMs['pmRA'], fregs_PMs['epmRA'], fregs_PMs['pmDE'],
            fregs_PMs['epmDE'])

    # All frame.
    allr_KDE_PMs = {}
    if allfr_PMs['pmRA'].any():
        allr_KDE_PMs = kde_2d(
            allfr_PMs['pmRA'], allfr_PMs['epmRA'], allfr_PMs['pmDE'],
            allfr_PMs['epmDE'], 100)

    return cr_KDE_PMs, fr_KDE_PMs, allr_KDE_PMs


def kde_2d(xarr, xsigma, yarr, ysigma, grid_dens=50, Nstd=3, N_max=10000):
    """
    Take an array of x,y data with their errors, create a grid of points in x,y
    and return the 2D KDE density map.

    For better performance, use max N_max random stars when estimating
    the KDE.
    """
    if xarr.size > N_max:
        # print(("  WARNING: used {} stars to estimate the PMs KDE\n"
        #        "  instead of the {} total").format(N_max, xarr.size))
        xarr, xsigma, yarr, ysigma =\
            np.random.choice(xarr, N_max, replace=False),\
            np.random.choice(xsigma, N_max, replace=False),\
            np.random.choice(yarr, N_max, replace=False),\
            np.random.choice(ysigma, N_max, replace=False)

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
    x, y = np.mgrid[xmin:xmax:gd_c, ymin:ymax:gd_c]
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

    return {'x': x, 'y': y, 'z': z, 'zmax_x': zmax_x, 'zmax_y': zmax_y}


# DEPRECATED 18/01/20
# def PMsNNAll(
#     PM_flag_all, flag_PM_coord, xRA_all, yDE_all, pmRA_all, pmDE_all,
#         nnmax):
#     """
#     Nearest-neighbor analysis that combines coordinates and PMs data.
#     """
#     if PM_flag_all:

#         # Normalize data
#         def norm(data):
#             data_norm = []
#             for arr in data:
#                 dmin = np.min(arr)
#                 arr_n = arr - dmin
#                 dmax = np.max(arr_n)
#                 arr_n /= dmax
#                 data_norm.append(arr_n)
#             return np.array(data_norm)

#         if flag_PM_coord is True:
#             data = norm([xRA_all, yDE_all, pmRA_all, pmDE_all]).T
#         else:
#             data = norm([pmRA_all, pmDE_all]).T

#         # Create the tree
#         tree = spatial.cKDTree(data)
#         # Find the closest nnmax-1 neighbors (1st entry is the point itself)
#         dists = tree.query(data, nnmax)
#         # Extract distances
#         nn_dist = dists[0][:, 1:]
#         # Median values
#         PMs_d_median = np.median(nn_dist, 1)

#         return PMs_d_median

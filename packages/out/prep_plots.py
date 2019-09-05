
from ..math_f import exp_function
from ..best_fit.obs_clust_prepare import dataProcess
from ..decont_algors.local_cell_clean import bin_edges_f
import numpy as np
import warnings
from scipy import stats
from astropy.visualization import ZScaleInterval
from astropy.stats import sigma_clipped_stats


def frame_max_min(x_data, y_data):
    '''
    Get max and min values in x,y coordinates.
    '''
    x_min, x_max = min(x_data), max(x_data)
    y_min, y_max = min(y_data), max(y_data)
    return x_min, x_max, y_min, y_max


def aspect_ratio(x_min, x_max, y_min, y_max):
    '''
    Define an optimal aspect ratio for full frame plots.
    '''
    x_range, y_range = abs(x_max - x_min), abs(y_max - y_min)
    asp_ratio = float(min(x_range, y_range) / max(x_range, y_range))

    # If the aspect ratio is smaller than 1:3.
    if asp_ratio < 0.33:
        print("  WARNING: frame's aspect ratio ({:.2f}) is < 1:3.\n"
              "  Cluster's plot will be stretched to 1:1.".format(asp_ratio))
        asp_ratio = 'auto'
    else:
        asp_ratio = 'equal'

    return asp_ratio


def coord_syst(coords):
    '''
    Define system of coordinates used.
    '''
    coord_lst = ['px', 'x', 'y'] if coords == 'px' else ['deg', 'ra', 'dec']
    return coord_lst


def frame_zoomed(x_min, x_max, y_min, y_max, kde_cent, clust_rad):
    '''
    If possible, define zoomed frame.
    '''
    x_zmin, x_zmax = max(x_min, (kde_cent[0] - 1.5 * clust_rad)), \
        min(x_max, (kde_cent[0] + 1.5 * clust_rad))
    y_zmin, y_zmax = max(y_min, (kde_cent[1] - 1.5 * clust_rad)), \
        min(y_max, (kde_cent[1] + 1.5 * clust_rad))
    # Prevent axis stretching.
    if (x_zmax - x_zmin) != (y_zmax - y_zmin):
        lst = [(x_zmax - x_zmin), (y_zmax - y_zmin)]
        val, idx = min((val, idx) for (idx, val) in enumerate(lst))
        if idx == 0:
            x_zmax = x_zmin + lst[1]
        else:
            y_zmax = y_zmin + lst[0]

    return x_zmin, x_zmax, y_zmin, y_zmax


def ax_names(x, y, yaxis):
    '''
    Define names for photometric diagram axes.
    '''
    col_n = []
    c_filts = x[1].split(',')
    for f in c_filts:
        if '_' in f:
            xs = f.split('_')
            col_n.append(xs[0] + '_{' + xs[1] + '}')
        else:
            col_n.append(f)
    x_ax = '(' + '-'.join(col_n) + ')'

    # yaxis indicates if the y axis is a magnitude or a color.
    if yaxis == 'mag':
        y_ax = y[1]
    else:
        y_ax = '(' + y[1].replace(',', '-') + ')'
    if '_' in y_ax:
        ys = y_ax.split('_')
        y_ax = ys[0] + '_{' + ys[1] + '}'

    return x_ax, y_ax


def diag_limits(yaxis, phot_x, phot_y):
    '''
    Define plot limits for *all* photometric diagrams.
    '''
    x_median, x_std = np.nanmedian(phot_x), np.nanstd(phot_x)
    x_min_cmd, x_max_cmd = x_median - 4.5 * x_std, x_median + 4.5 * x_std

    # Use stars within the x limits defined. This prevents stars far away
    # from the x median from affecting the limit in y.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        xmsk = (phot_x > x_min_cmd) & (phot_x < x_max_cmd)

    phot_y_msk = np.array(phot_y)[xmsk]
    y_median, y_std = np.nanmedian(phot_y_msk), np.nanstd(phot_y_msk)

    # y limits.
    if yaxis == 'mag':
        y_min_cmd = y_median + 1.25 * y_std + .75
        # If photometric axis y is a magnitude, make sure the brightest star
        # is always plotted.
        y_max_cmd = np.nanmin(phot_y_msk) - 1.
    else:
        y_max_cmd, y_min_cmd = y_median - 4.5 * y_std, y_median + 4.5 * y_std

    return x_max_cmd, x_min_cmd, y_min_cmd, y_max_cmd


def star_size(mag, N=None, zmin=None, zmax=None):
    '''
    Convert magnitudes into intensities and define sizes of stars in
    finding chart.
    '''
    mag = np.array(mag)
    if N is None:
        N = mag.size
    if zmin is None and zmax is None:
        interval = ZScaleInterval()
        zmin, zmax = interval.get_limits(mag)

    mag = mag.clip(zmin, zmax)
    factor = 500. * (1 - 1 / (1 + 150 / N ** 0.85))
    sizes = .1 + factor * (10 ** ((mag - zmin) / -2.5))
    return sizes


def phot_diag_st_size(x):
    '''
    Calculate optimal size for stars in photometric diagram.
    '''
    a, b, c, d = 2.99, -2.81, 563.36, 15.02
    if x != 0:
        return ((a - d) / (1 + ((x / c) ** b))) + d
    else:
        # If no field regions were defined.
        return 0.


def zoomed_frame(x, y, mags, x_zmin, x_zmax, y_zmin, y_zmax):
    '''
    Separate stars for zoomed frame. Use main magnitude.
    '''
    x_data_z, y_data_z, mag_data_z = [], [], []
    for st_x, st_y, st_mag in list(zip(x, y, mags[0])):
        if x_zmin <= st_x <= x_zmax and y_zmin <= st_y <= y_zmax:
            x_data_z.append(st_x)
            y_data_z.append(st_y)
            mag_data_z.append(st_mag)

    return x_data_z, y_data_z, mag_data_z


def da_colorbar_range(cl_reg_fit, cl_reg_no_fit):
    """
    Extreme values for colorbar.
    """
    lst_comb = cl_reg_fit + cl_reg_no_fit
    v_min_mp, v_max_mp = round(min(list(zip(*lst_comb))[-1]), 2), \
        round(max(list(zip(*lst_comb))[-1]), 2)

    return v_min_mp, v_max_mp


def da_find_chart(
    kde_cent, clust_rad, stars_out, x_zmin, x_zmax, y_zmin, y_zmax,
        cl_reg_fit, cl_reg_no_fit):
    '''
    Finding chart with MPs assigned by the DA.
    '''
    # Arrange stars used in the best fit process.
    cl_reg_fit = list(zip(*cl_reg_fit))
    # Finding chart data. Invert values so higher prob stars are on top.
    chart_fit_inv = [i[::-1] for i in
                     [cl_reg_fit[1], cl_reg_fit[2], cl_reg_fit[9]]]

    # Arrange stars *not* used in the best fit process.
    if cl_reg_no_fit:
        cl_reg_no_fit = list(zip(*cl_reg_no_fit))
        # Finding chart data.
        chart_no_fit_inv = [
            i[::-1] for i in [cl_reg_no_fit[1], cl_reg_no_fit[2],
                              cl_reg_no_fit[9]]]
    else:
        chart_no_fit_inv = [[], [], []]

    # Separate stars outside the cluster's radius.
    out_clust_rad = [[], []]
    for star in stars_out:
        if x_zmin <= star[1] <= x_zmax and y_zmin <= star[2] <= y_zmax:
            dist = np.sqrt((kde_cent[0] - star[1]) ** 2 +
                           (kde_cent[1] - star[2]) ** 2)
            # Only plot stars outside the cluster's radius.
            if dist >= clust_rad:
                out_clust_rad[0].append(star[1])
                out_clust_rad[1].append(star[2])

    return chart_fit_inv, chart_no_fit_inv, out_clust_rad


def da_phot_diag(cl_reg_fit, cl_reg_no_fit):
    '''
    Generate parameters for the photometric diagram plotted with the MPs
    assigned by the DA. The stars are inverted according to their MPs, so that
    those with larger probabilities are plotted last.
    '''
    # Arrange stars used in the best fit process.
    cl_reg_fit = list(zip(*cl_reg_fit))
    # Magnitudes.
    diag_fit_inv = [[i[::-1] for i in list(zip(*cl_reg_fit[3]))]]
    # Colors.
    diag_fit_inv += [[i[::-1] for i in list(zip(*cl_reg_fit[5]))]]
    # membership probabilities.
    diag_fit_inv += [cl_reg_fit[9][::-1]]

    # Arrange stars *not* used in the best fit process.
    if cl_reg_no_fit:
        cl_reg_no_fit = list(zip(*cl_reg_no_fit))
        # Magnitudes.
        diag_no_fit_inv = [[i[::-1] for i in list(zip(*cl_reg_no_fit[3]))]]
        # Colors.
        diag_no_fit_inv += [[i[::-1] for i in list(zip(*cl_reg_no_fit[5]))]]
        # membership probabilities.
        diag_no_fit_inv += [cl_reg_no_fit[9][::-1]]
    else:
        diag_no_fit_inv = [[[]], [[]], []]

    return diag_fit_inv, diag_no_fit_inv


def error_bars(stars_phot, x_min_cmd, err_lst, all_flag=None):
    """
    Calculate error bars for plotting in photometric diagram.
    """
    # Use main magnitude.
    if all_flag == 'all':
        mmag = np.array(stars_phot)
    else:
        if stars_phot:
            mmag = np.array(list(zip(*list(zip(*stars_phot))[3]))[0])
        else:
            mmag = np.array([])

    x_val, mag_y, xy_err = [], [], []
    if mmag.any():
        # List of y values where error bars are plotted.
        mag_y = np.arange(
            int(min(mmag) + 0.5), int(max(mmag) + 0.5) + 0.1)
        # List of x values where error bars are plotted.
        x_val = [x_min_cmd + 0.15] * len(mag_y)
        # Read average fitted values for exponential error fit.
        # Magnitude values are positioned first and colors after in the list
        # 'err_lst'.
        # popt_mag = err_lst[0]
        for popt in err_lst:
            xy_err.append(exp_function.exp_3p(mag_y, *popt))

    return [x_val, mag_y, xy_err]


def param_ranges(best_fit_algor, fundam_params, varIdxs=None, trace=None):
    '''
    Parameter ranges used by several plots.
    '''
    min_max_p = []
    if best_fit_algor in ['brute', 'boot+GA']:

        if varIdxs is not None and trace is not None:
            for cp, param in enumerate(fundam_params):
                if cp in varIdxs:
                    c_model = varIdxs.index(cp)
                    # Use the last 10% of the trace.
                    N = int(trace[c_model].shape[-1] * .1)
                    std = np.std(trace[c_model][-N:])
                    pmin, pmax = np.min(trace[c_model][-N:]),\
                        np.max(trace[c_model][-N:])
                    min_max_p.append([
                        max(param[0], pmin - std),
                        min(param[-1], pmax + std)])
                else:
                    min_max_p.append([min(param) - .001, max(param) + .001])
        else:
            for param in fundam_params:
                # Set the delta for the parameter range. If only one value was
                # used, set a very small delta value.
                delta_p = (max(param) - min(param)) * 0.05 \
                    if max(param) != min(param) else 0.001
                # Store parameter range.
                min_max_p.append([min(param) - delta_p, max(param) + delta_p])

    elif best_fit_algor in ('ptemcee', 'emcee'):
        # Select the ranges given by the limits of the space explored by all
        # the chains, for each parameter.
        for cp, param in enumerate(fundam_params):
            if cp in varIdxs:
                c_model = varIdxs.index(cp)
                # Use the last 10% of the chains.
                N = int(trace[c_model].shape[-1] * .1)
                std = np.std(trace[c_model][:, -N:])
                pmin, pmax = np.min(trace[c_model][:, -N:]),\
                    np.max(trace[c_model][:, -N:])
                min_max_p.append([
                    max(param[0], pmin - std),
                    min(param[-1], pmax + std)])
            else:
                min_max_p.append([min(param) - .001, max(param) + .001])

    elif best_fit_algor == 'abc':
        # Select the ranges given by the limits of the space explored by all
        # the chains, for each parameter.
        for cp, param in enumerate(fundam_params):
            if cp in varIdxs:
                c_model = varIdxs.index(cp)
                mean = np.mean(post_bi[c_model])
                std3 = 3 * np.std(post_bi[c_model])
                min_max_p.append([
                    max(param[0], mean - std3),
                    min(param[-1], mean + std3)])
            else:
                min_max_p.append([min(param) - .001, max(param) + .001])

    return min_max_p


def p2_ranges(p2, min_max_p):
    '''
    Parameter ranges used by the MCMC 2-param density plots.
    '''
    par_idx = {
        'metal': 0, 'age': 1, 'ext': 2, 'dist': 3, 'mass': 4, 'binar': 5}
    par = p2.split('-')

    min_max_p2 = min_max_p[par_idx[par[0]]] + min_max_p[par_idx[par[1]]]

    return min_max_p2


def likl_y_range(opt_method, lkl_best, lkl_mean):
    '''
    Obtain y axis range for the likelihood axis.
    '''
    # if opt_method == 'emcee':
    #     l_min_max = [
    #         max(0., min(lkl_old) - .2 * min(lkl_old)),
    #         np.median(lkl_old[:int(.1 * len(lkl_old))]) * 1.5]
    # elif opt_method == 'boot+GA':
    # Take limits from L_min curve.
    lkl_range = max(lkl_mean) - min(lkl_best)
    l_min_max = [max(0., min(lkl_best) - 0.1 * lkl_range),
                 max(lkl_mean) + 0.1 * lkl_range]

    return l_min_max


def packData(
    lkl_method, lkl_binning, cl_max_mag, synth_clst, shift_isoch,
        colors, filters, col_0_comb, mag_0_comb, col_1_comb):
    """
    Properly select and pack data for CMD/CCD of observed and synthetic
    clusters, and their Hess diagram.
    """
    bin_method = 'auto' if lkl_method == 'tolstoy' else lkl_binning
    mags_cols_cl, dummy = dataProcess(cl_max_mag)
    # Obtain bin edges for each dimension, defining a grid.
    bin_edges = bin_edges_f(bin_method, mags_cols_cl)

    N_mags, N_cols = len(filters), len(colors)

    # CMD of main magnitude and first color defined.
    # Used to defined limits.
    x_phot_all, y_phot_all = col_0_comb, mag_0_comb
    frst_obs_mag, frst_obs_col = list(zip(*list(zip(*cl_max_mag))[3]))[0],\
        list(zip(*list(zip(*cl_max_mag))[5]))[0]
    frst_synth_col, frst_synth_mag = synth_clst[0][0][1],\
        synth_clst[0][0][0]
    # Indexes of binary systems.
    binar_idx = synth_clst[1][0]
    frst_col_edgs, frst_mag_edgs = bin_edges[1], bin_edges[0]
    # Filters and colors are appended continuously in 'shift_isoch'. If
    # there are 3 defined filters, then the first color starts at the
    # index 3. This is why 'N_mags' is used as the 'first color' index.
    frst_col_isoch, frst_mag_isoch = shift_isoch[N_mags], shift_isoch[0]
    # gs plot coords.
    gs_y1, gs_y2 = 0, 2
    # Index of observed filter/color to scatter plot.
    i_obs_x, i_obs_y = 0, 0
    hr_diags = [
        [x_phot_all, y_phot_all, frst_obs_col, frst_obs_mag, frst_synth_col,
         frst_synth_mag, binar_idx, frst_col_edgs, frst_mag_edgs,
         frst_col_isoch, frst_mag_isoch, colors[0], filters[0], 'mag',
         i_obs_x, i_obs_y, gs_y1, gs_y2]]

    # If more than one color was defined, plot an extra CMD (main magnitude
    # versus first color), and an extra CCD (first color versus second color)
    if N_cols > 1:
        scnd_obs_col = list(zip(*list(zip(*cl_max_mag))[5]))[1]
        scnd_synth_col = synth_clst[0][0][2]
        scnd_col_edgs = bin_edges[2]
        scnd_col_isoch = shift_isoch[N_mags + 1]
        # CMD of main magnitude and second color defined.
        x_phot_all, y_phot_all = col_1_comb, mag_0_comb
        gs_y1, gs_y2 = 2, 4
        i_obs_x, i_obs_y = 1, 0
        hr_diags.append(
            [x_phot_all, y_phot_all, scnd_obs_col, frst_obs_mag,
             scnd_synth_col, frst_synth_mag, binar_idx, scnd_col_edgs,
             frst_mag_edgs, shift_isoch[2], frst_mag_isoch, colors[1],
             filters[0], 'mag', i_obs_x, i_obs_y, gs_y1, gs_y2])
        # CCD of first and second color defined.
        x_phot_all, y_phot_all = col_0_comb, col_1_comb
        gs_y1, gs_y2 = 4, 6
        i_obs_x, i_obs_y = 0, 1
        hr_diags.append(
            [x_phot_all, y_phot_all, frst_obs_col, scnd_obs_col,
             frst_synth_col, scnd_synth_col, binar_idx, frst_col_edgs,
             scnd_col_edgs, frst_col_isoch, scnd_col_isoch, colors[0],
             colors[1], 'col', i_obs_x, i_obs_y, gs_y1, gs_y2])

    return hr_diags


def get_hess(obs_mags_cols, synth_phot, hess_xedges, hess_yedges):
    """
    Hess diagram of observed minus best match synthetic cluster.
    """
    # 2D histogram of the observed cluster.
    cl_histo = np.histogram2d(
        *obs_mags_cols, bins=[hess_xedges, hess_yedges])[0]
    # 2D histogram of the synthetic cluster.
    syn_histo = np.histogram2d(*synth_phot, bins=[hess_xedges, hess_yedges])[0]

    # Grid for pcolormesh.
    hess_x, hess_y = np.meshgrid(hess_xedges, hess_yedges)

    # Hess diagram: observed minus synthetic.
    hess_diag = np.array([])
    if syn_histo.size:
        hess_diag = cl_histo - syn_histo
        if hess_diag.size:
            HD = np.rot90(hess_diag)
            HD = np.flipud(HD)
        else:
            HD = np.array([])

    if not HD.size:
        print("  WARNING: the Hess diagram could no be obtained")

    return hess_x, hess_y, HD


def plxPlot(
    mmag_clp, mp_clp, plx_clp, e_plx_clp, flag_no_fl_regs_i,
        field_regions_i):
    """
    Parameters for the parallax plot.
    """
    # Put large MP stars in cluster region on top.
    mp_i = mp_clp.argsort()
    mmag_clp, mp_clp, plx_clp, e_plx_clp = mmag_clp[mp_i],\
        mp_clp[mp_i], plx_clp[mp_i], e_plx_clp[mp_i]

    if not flag_no_fl_regs_i:
        plx_flrg, mag_flrg = [], []
        # Extract parallax data.
        for fl_rg in field_regions_i:
            plx_flrg += list(zip(*list(zip(*fl_rg))[7]))[0]
            mag_flrg += list(zip(*list(zip(*fl_rg))[3]))[0]
        plx_flrg, mag_flrg = np.asarray(plx_flrg), np.asarray(mag_flrg)
        # Mask 'nan' and set range.
        msk0 = ~np.isnan(plx_flrg)
        plx_flrg, mag_flrg = plx_flrg[msk0], mag_flrg[msk0]
        msk = (plx_flrg > -5.) & (plx_flrg < 10.)
        plx_flrg, mag_flrg = plx_flrg[msk], mag_flrg[msk]
    else:
        plx_flrg, mag_flrg = np.array([]), []

    return plx_flrg, mag_flrg, mmag_clp, mp_clp, plx_clp, e_plx_clp


def PMsPlot(pmMP, pmRA_DE, e_pmRA_DE, pmDE, e_pmDE, mmag_pm, pm_dist_max):
    """
    Parameters for the proper motions plot.
    """

    # Re-arrange so stars with larger MPs are on top.
    mp_i = pmMP.argsort()
    arr_lst = []
    for arr in (pmMP, pmRA_DE, e_pmRA_DE, pmDE, e_pmDE, mmag_pm, pm_dist_max):
        if arr is not None:
            arr_lst.append(arr[mp_i])
        else:
            arr_lst.append([])

    return arr_lst


def CIEllipse(points, prob=.95):
    """
    Generate a 'prob' confidence interval ellipse based on the mean and
    covariance of a point "cloud".

    Source: https://stackoverflow.com/q/12301071/1391441
    Definition: https://stats.stackexchange.com/a/217377/10416

    Definition (Wikipedia): "Were this procedure to be repeated on numerous
    samples, the fraction of calculated confidence intervals (which would
    differ for each sample) that encompass the true population parameter would
    tend toward 90%."


    Parameters
    ----------
        points : An Nx2 array of the data points.
        prob : probability value for the CI region.
    """
    def eigsorted(cov):
        '''
        Eigenvalues and eigenvectors of the covariance matrix.
        '''
        vals, vecs = np.linalg.eigh(cov)
        order = vals.argsort()[::-1]
        return vals[order], vecs[:, order]

    # Location of the center of the ellipse.
    mean_pos = points.mean(axis=0)

    # The 2x2 covariance matrix to base the ellipse on.
    cov = np.cov(points, rowvar=False)

    vals, vecs = eigsorted(cov)
    theta = np.degrees(np.arctan2(*vecs[:, 0][::-1]))

    k = np.sqrt(stats.chi2.ppf(prob, 2))
    # Width and height are "full" widths, not radius
    width, height = 2 * np.sqrt(vals) * k

    return mean_pos, width, height, theta


def complSeparate(cl_region_i, membs_i, memb_prob_avrg_sort):
    """
    Separate stars into complete and incomplete arrays.
    """
    mags_c, mags_i, cols_c, cols_i, colors_c, colors_i = [[] for _ in range(6)]
    ids_c = list(zip(*memb_prob_avrg_sort))[0]
    for i, star in enumerate(cl_region_i):
        id_i = star[0]
        if id_i in ids_c:
            j = ids_c.index(id_i)
            mags_c.append(memb_prob_avrg_sort[j][3][0])
            cols_c.append(memb_prob_avrg_sort[j][5])
            colors_c.append(memb_prob_avrg_sort[j][9])
        else:
            mags_i.append(star[3][0])
            cols_i.append(star[5])
            colors_i.append(membs_i[i])

    cols_c = np.array(cols_c).T
    if cols_i:
        cols_i = np.array(cols_i).T
        idx_i = np.argsort(colors_i)
        mags_i = np.array(mags_i)[idx_i].tolist()
        cols_i = np.array([_[idx_i] for _ in cols_i]).tolist()
        colors_i = np.array(colors_i)[idx_i].tolist()
    else:
        cols_i = [[] for _ in range(cols_c.shape[0])]

    idx_c = np.argsort(colors_c)
    mags_c = np.array(mags_c)[idx_c].tolist()
    cols_c = np.array([_[idx_c] for _ in cols_c]).tolist()
    colors_c = np.array(colors_c)[idx_c].tolist()

    return mags_c, mags_i, cols_c, cols_i, colors_c, colors_i


def PMsrange(pmRA_DE, pmDE):
    """
    """
    ra_mean, ra_median, ra_std = sigma_clipped_stats(pmRA_DE)
    de_mean, de_median, de_std = sigma_clipped_stats(pmDE)
    x_range = min(6. * ra_std, np.ptp(pmRA_DE))
    y_range = min(6. * de_std, np.ptp(pmDE))
    xyrange = max(x_range, y_range)

    raPMrng = ra_median - .5 * xyrange, ra_median + .5 * xyrange
    dePMrng = de_median - .5 * xyrange, de_median + .5 * xyrange

    return raPMrng, dePMrng

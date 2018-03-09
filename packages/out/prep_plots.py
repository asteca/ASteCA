
from ..math_f import exp_function
from ..best_fit.obs_clust_prepare import dataProcess
from ..decont_algors.local_cell_clean import bin_edges_f
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt


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
        print ("  WARNING: frame's aspect ratio ({:.2f}) is < 1:3.\n"
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
    # Create photometric axis names.
    x_ax = '(' + x[1].replace(',', '-') + ')'
    # yaxis indicates if the y axis is a magnitude or a color.
    if yaxis == 'mag':
        y_ax = y[1]
    else:
        y_ax = '(' + y[1].replace(',', '-') + ')'
    return x_ax, y_ax


def kde_limits(phot_x, phot_y):
    '''
    Return photometric diagram limits taken from a 2D KDE.
    '''

    xmin, xmax = min(phot_x), max(phot_x)
    ymin, ymax = min(phot_y), max(phot_y)
    # Stack photometric data.
    values = np.vstack([phot_x, phot_y])
    # Obtain Gaussian KDE.
    kernel = stats.gaussian_kde(values)
    # Grid density (number of points).
    gd = 25
    gd_c = complex(0, gd)
    # Define x,y grid.
    x, y = np.mgrid[xmin:xmax:gd_c, ymin:ymax:gd_c]
    positions = np.vstack([x.ravel(), y.ravel()])
    # Evaluate kernel in grid positions.
    k_pos = kernel(positions)

    # Generate 30 contour lines.
    plt.figure()
    cs = plt.contour(x, y, np.reshape(k_pos, x.shape), 30)
    plt.close()
    # Extract (x,y) points delimiting each line.
    x_v, y_v = np.asarray([]), np.asarray([])
    # Only use the outer curve.
    col = cs.collections[0]
    # If more than one region is defined by this curve (ie: the main sequence
    # region plus a RC region or some other detached region), obtain x,y from
    # all of them.
    for lin in col.get_paths():
        x_v = np.append(x_v, lin.vertices[:, 0])
        y_v = np.append(y_v, lin.vertices[:, 1])

    return x_v, y_v


def diag_limits(yaxis, phot_x, phot_y):
    '''
    Define plot limits for *all* photometric diagrams.
    '''
    x_v, y_v = kde_limits(phot_x, phot_y)

    # Define diagram limits.
    x_min_cmd, x_max_cmd = min(x_v) - 1.25, max(x_v) + 1.25
    y_min_cmd = max(y_v) + 1.25
    # If photometric axis y is a magnitude, make sure the brightest star
    # is always plotted.
    if yaxis == 'mag':
        y_max_cmd = min(phot_y) - 1.
    else:
        y_max_cmd = min(y_v) - 1.

    return x_max_cmd, x_min_cmd, y_min_cmd, y_max_cmd


def star_size(mag, N=None, min_m=None):
    '''
    Convert magnitudes into intensities and define sizes of stars in
    finding chart.
    '''
    # Scale factor.
    if N is None:
        N = len(mag)
    if min_m is None:
        min_m = min(mag)
    factor = 500. * (1 - 1 / (1 + 150 / N ** 0.85))
    return 0.1 + factor * 10 ** ((np.array(mag) - min_m) / -2.5)


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
    for st_x, st_y, st_mag in zip(x, y, mags[0]):
        if x_zmin <= st_x <= x_zmax and y_zmin <= st_y <= y_zmax:
            x_data_z.append(st_x)
            y_data_z.append(st_y)
            mag_data_z.append(st_mag)

    return x_data_z, y_data_z, mag_data_z


def field_region_stars(stars_out_rjct, field_regions):
    """
    Generate list with *all* rejected stars outside of the cluster region, and
    all stars within a defined field region.
    """
    stars_f_rjct = [[], []]
    for star in stars_out_rjct:
        stars_f_rjct[0].append(star[5][0])
        stars_f_rjct[1].append(star[3][0])

    stars_f_acpt = [[], []]
    if field_regions:
        for fr in field_regions:
            for star in fr:
                stars_f_acpt[0].append(star[5][0])
                stars_f_acpt[1].append(star[3][0])

    return stars_f_rjct, stars_f_acpt


def da_colorbar_range(cl_reg_fit, cl_reg_no_fit):
    """
    Extreme values for colorbar.
    """
    lst_comb = cl_reg_fit + cl_reg_no_fit
    v_min_mp, v_max_mp = round(min(zip(*lst_comb)[-1]), 2), \
        round(max(zip(*lst_comb)[-1]), 2)

    return v_min_mp, v_max_mp


def da_find_chart(
    kde_cent, clust_rad, stars_out, x_zmin, x_zmax, y_zmin, y_zmax,
        cl_reg_fit, cl_reg_no_fit):
    '''
    Finding chart with MPs assigned by the DA.
    '''
    # Arrange stars used in the best fit process.
    cl_reg_fit = zip(*cl_reg_fit)
    # Finding chart data. Invert values so higher prob stars are on top.
    chart_fit_inv = [i[::-1] for i in
                     [cl_reg_fit[1], cl_reg_fit[2], cl_reg_fit[7]]]

    # Arrange stars *not* used in the best fit process.
    if cl_reg_no_fit:
        cl_reg_no_fit = zip(*cl_reg_no_fit)
        # Finding chart data.
        chart_no_fit_inv = [
            i[::-1] for i in [cl_reg_no_fit[1], cl_reg_no_fit[2],
                              cl_reg_no_fit[7]]]
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


def da_phot_diag(cl_reg_fit, cl_reg_no_fit, v_min_mp, v_max_mp):
    '''
    Generate parameters for the photometric diagram plotted with the MPs
    assigned by the DA. The stars are inverted according to their MPs, so that
    those with larger probabilities are plotted last.
    '''
    # Decide if colorbar should be plotted.
    plot_colorbar = True if v_min_mp != v_max_mp else False

    # Arrange stars used in the best fit process.
    cl_reg_fit = zip(*cl_reg_fit)
    # Magnitudes.
    diag_fit_inv = [[i[::-1] for i in zip(*cl_reg_fit[3])]]
    # Colors.
    diag_fit_inv += [[i[::-1] for i in zip(*cl_reg_fit[5])]]
    # membership probabilities.
    diag_fit_inv += [cl_reg_fit[7][::-1]]

    # Arrange stars *not* used in the best fit process.
    if cl_reg_no_fit:
        cl_reg_no_fit = zip(*cl_reg_no_fit)
        # Magnitudes.
        diag_no_fit_inv = [[i[::-1] for i in zip(*cl_reg_no_fit[3])]]
        # Colors.
        diag_no_fit_inv += [[i[::-1] for i in zip(*cl_reg_no_fit[5])]]
        # membership probabilities.
        diag_no_fit_inv += [cl_reg_no_fit[7][::-1]]
    else:
        diag_no_fit_inv = [[[]], [[]], []]

    return plot_colorbar, diag_fit_inv, diag_no_fit_inv


def error_bars(stars_phot, x_min_cmd, err_lst, all_flag=None):
    """
    Calculate error bars for plotting in photometric diagram.
    """
    # Use main magnitude.
    if all_flag == 'all':
        mmag = stars_phot.tolist()
    else:
        mmag = zip(*zip(*stars_phot)[3])[0]
    x_val, mag_y, x_err, y_err = [], [], [], []
    if mmag:
        # List of y values where error bars are plotted.
        mag_y = np.arange(
            int(min(mmag) + 0.5), int(max(mmag) + 0.5) + 0.1)
        # List of x values where error bars are plotted.
        x_val = [x_min_cmd + 0.4] * len(mag_y)
        # Read average fitted values for exponential error fit.
        # Magnitude values are positioned first and colors after in the list
        # 'err_lst'.
        # TODO generalize to N dimensions
        popt_mag, popt_col1 = err_lst[0], err_lst[1]
        x_err = exp_function.exp_3p(mag_y, *popt_col1)
        y_err = exp_function.exp_3p(mag_y, *popt_mag)
    err_bar = [x_val, mag_y, x_err, y_err]

    return err_bar


def param_ranges(fundam_params):
    '''
    Parameter ranges used by several plots.
    '''
    min_max_p = []
    for param in fundam_params:
        # Set the delta for the parameter range. If only one value was
        # used, set a very small delta value.
        delta_p = (max(param) - min(param)) * 0.05 \
            if max(param) != min(param) else 0.001
        # Store parameter range.
        min_max_p.append([min(param) - delta_p, max(param) + delta_p])

    return min_max_p


def p2_ranges(min_max_p, varIdxs, model_done, nwalkers, nsteps):
    '''
    Parameter ranges used by 'emcee' 2-param density plots.
    '''
    min_max_p2 = []
    for vi in range(6):  # TODO hard-coded to 6 parameters
        if vi in varIdxs:
            model = varIdxs.index(vi)
            hx, edge = np.histogram(
                model_done[model], range=min_max_p[vi], bins=50)
            # x% of the total number of synthetic clusters explored.
            non_z = np.nonzero(hx > int(.005 * nwalkers * nsteps))
            min_max_p2.append([edge[non_z[0][0]], edge[non_z[0][-1]]])
        else:
            min_max_p2.append(min_max_p[vi])

    return min_max_p2


def likl_y_range(opt_method, lkl_old):
    '''
    Obtain y axis range for the likelihood axis.
    '''
    if opt_method == 'emcee':
        l_min_max = [
            max(0., min(lkl_old) - .2 * min(lkl_old)),
            np.median(lkl_old[:int(.1 * len(lkl_old))]) * 1.5]
    elif opt_method == 'genet':
        # Take limits from L_min curve.
        lkl_range = max(lkl_old[1]) - min(lkl_old[0])
        l_min_max = [max(0., min(lkl_old[0]) - 0.1 * lkl_range),
                     max(lkl_old[1]) + 0.1 * lkl_range]

    return l_min_max


# def BestTick(minv, maxv, max_char):
#     '''
#     Find optimal number and length of ticks for a given fixed maximum
#     number of characters in the axis.
#     '''

#     st, diff_chars, st_indx = [], 1000, 0
#     # Check these 4 possible sizes for the ticks and keep the best one.
#     for i in range(4):
#         mostticks = i + 4

#         minimum = (maxv - minv) / mostticks
#         magnitude = 10 ** math.floor(math.log(minimum) / math.log(10))
#         residual = minimum / magnitude
#         if residual > 5:
#             tick = 10 * magnitude
#         elif residual > 2:
#             tick = 5 * magnitude
#         elif residual > 1:
#             tick = 2 * magnitude
#         else:
#             tick = magnitude

#         st.append(tick)
#         # Count the number of chars used by this step.
#         ms = (i + 4) * (len(str(tick)) - 1)
#         # Only use if it is less than the fixed max value of chars.
#         if ms <= max_char:
#             if (max_char - ms) < diff_chars:
#                 # Store the closest value to max_chars.
#                 diff_chars = (max_char - ms)
#                 st_indx = i

#     # Set min tick value according to the best step length selected above.
#     if minv <= 0.:
#         xmin = 0.
#     elif minv <= st[st_indx]:
#         xmin = st[st_indx]
#     else:
#         xmin = int(round(minv / st[st_indx])) * st[st_indx]

#     return xmin, st[st_indx]


def packData(lkl_method, lkl_binning, cl_max_mag, synth_clst, shift_isoch,
             colors, filters, cld):
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
    x_phot_all, y_phot_all = cld['cols'][0], cld['mags'][0]
    frst_obs_mag, frst_obs_col = list(zip(*zip(*cl_max_mag)[3])[0]),\
        list(zip(*zip(*cl_max_mag)[5])[0])
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
        scnd_obs_col = list(zip(*zip(*cl_max_mag)[5])[1])
        scnd_synth_col = synth_clst[0][0][2]
        scnd_col_edgs = bin_edges[2]
        scnd_col_isoch = shift_isoch[N_mags + 1]
        # CMD of main magnitude and second color defined.
        x_phot_all, y_phot_all = cld['cols'][1], cld['mags'][0]
        gs_y1, gs_y2 = 2, 4
        i_obs_x, i_obs_y = 1, 0
        hr_diags.append(
            [x_phot_all, y_phot_all, scnd_obs_col, frst_obs_mag,
             scnd_synth_col, frst_synth_mag, binar_idx, scnd_col_edgs,
             frst_mag_edgs, shift_isoch[2], frst_mag_isoch, colors[1],
             filters[0], 'mag', i_obs_x, i_obs_y, gs_y1, gs_y2])
        # CCD of first and second color defined.
        x_phot_all, y_phot_all = cld['cols'][0], cld['cols'][1]
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
        print("  WARNING: the Hess diagram could no be obtained.")

    return hess_x, hess_y, HD

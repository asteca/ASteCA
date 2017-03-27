
from ..math_f import exp_function
from ..best_fit import obs_clust_prepare
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import math


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


def frame_zoomed(x_min, x_max, y_min, y_max, clust_cent, clust_rad):
    '''
    If possible, define zoomed frame.
    '''
    x_zmin, x_zmax = max(x_min, (clust_cent[0] - 1.5 * clust_rad)), \
        min(x_max, (clust_cent[0] + 1.5 * clust_rad))
    y_zmin, y_zmax = max(y_min, (clust_cent[1] - 1.5 * clust_rad)), \
        min(y_max, (clust_cent[1] + 1.5 * clust_rad))
    # Prevent axis stretching.
    if (x_zmax - x_zmin) != (y_zmax - y_zmin):
        lst = [(x_zmax - x_zmin), (y_zmax - y_zmin)]
        val, idx = min((val, idx) for (idx, val) in enumerate(lst))
        if idx == 0:
            x_zmax = x_zmin + lst[1]
        else:
            y_zmax = y_zmin + lst[0]

    return x_zmin, x_zmax, y_zmin, y_zmax


def ax_names(filters, colors):
    '''
    Define names for photometric diagram axes.
    '''
    # y_axis == 0 indicates that the y axis is a magnitude.
    y_axis = 0
    # Create photometric axis names.
    x_ax = '(' + colors[0][1].replace(',', '-') + ')'
    y_ax = filters[0][1]
    return x_ax, y_ax, y_axis


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


def diag_limits(y_axis, phot_x, phot_y):
    '''
    Define plot limits for *all* photometric diagrams.
    '''
    x_v, y_v = kde_limits(phot_x[0], phot_y[0])

    # Define diagram limits.
    x_min_cmd, x_max_cmd = min(x_v) - 1.25, max(x_v) + 1.25
    y_min_cmd = max(y_v) + 1.25
    # If photometric axis y is a magnitude, make sure the brightest star
    # is always plotted.
    if y_axis == 0:
        y_max_cmd = min(phot_y[0]) - 1.
    else:
        y_max_cmd = min(y_v) - 1.

    return x_max_cmd, x_min_cmd, y_min_cmd, y_max_cmd


def star_size(mag):
    '''
    Convert magnitudes into intensities and define sizes of stars in
    finding chart.
    '''
    # Scale factor.
    factor = 500. * (1 - 1 / (1 + 150 / len(mag) ** 0.85))
    return 0.1 + factor * 10 ** ((np.array(mag) - min(mag)) / -2.5)


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
    clust_cent, clust_rad, stars_out, x_zmin, x_zmax, y_zmin, y_zmax,
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
            dist = np.sqrt((clust_cent[0] - star[1]) ** 2 +
                           (clust_cent[1] - star[2]) ** 2)
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
    # Decides if colorbar should be plotted.
    plot_colorbar = True if v_min_mp != v_max_mp else False

    # Arrange stars used in the best fit process.
    cl_reg_fit = zip(*cl_reg_fit)
    # Photometric diagram.
    diag_fit_inv = [
        i[::-1] for i in [zip(*cl_reg_fit[5])[0], zip(*cl_reg_fit[3])[0],
                          cl_reg_fit[7]]]

    # Arrange stars *not* used in the best fit process.
    if cl_reg_no_fit:
        cl_reg_no_fit = zip(*cl_reg_no_fit)
        # Photometric diagram.
        diag_no_fit_inv = [
            i[::-1] for i in [
                zip(*cl_reg_no_fit[5])[0], zip(*cl_reg_no_fit[3])[0],
                cl_reg_no_fit[7]]]
    else:
        diag_no_fit_inv = [[], [], []]

    return plot_colorbar, diag_fit_inv, diag_no_fit_inv


def error_bars(stars_phot, x_min_cmd, err_lst):
    """
    Calculate error bars for plotting in photometric diagram.
    """
    # Use main magnitude.
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
    Parameter ranges used by GA plots.
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


def likl_y_range(lkl_old):
    '''
    Obtain y axis range for the likelihood axis.
    '''
    # Take limits from L_min curve.
    lkl_range = max(lkl_old[1]) - min(lkl_old[0])
    l_min_max = [max(0., min(lkl_old[0]) - 0.1 * lkl_range),
                 max(lkl_old[1]) + 0.1 * lkl_range]

    return l_min_max


def BestTick(minv, maxv, max_char):
    '''
    Find optimal number and length of ticks for a given fixed maximum
    number of characters in the axis.
    '''

    st, diff_chars, st_indx = [], 1000, 0
    # Check these 4 possible sizes for the ticks and keep the best one.
    for i in range(4):
        mostticks = i + 4

        minimum = (maxv - minv) / mostticks
        magnitude = 10 ** math.floor(math.log(minimum) / math.log(10))
        residual = minimum / magnitude
        if residual > 5:
            tick = 10 * magnitude
        elif residual > 2:
            tick = 5 * magnitude
        elif residual > 1:
            tick = 2 * magnitude
        else:
            tick = magnitude

        st.append(tick)
        # Count the number of chars used by this step.
        ms = (i + 4) * (len(str(tick)) - 1)
        # Only use if it is less than the fixed max value of chars.
        if ms <= max_char:
            if (max_char - ms) < diff_chars:
                # Store the closest value to max_chars.
                diff_chars = (max_char - ms)
                st_indx = i

    # Set min tick value according to the best step length selected above.
    if minv <= 0.:
        xmin = 0.
    elif minv <= st[st_indx]:
        xmin = st[st_indx]
    else:
        xmin = int(round(minv / st[st_indx])) * st[st_indx]

    return xmin, st[st_indx]


def get_hess(lkl_method, bin_method, cl_reg_fit, synth_clust):
    """
    Hess diagram of observed minus best match synthetic cluster.
    """
    if lkl_method == 'tolstoy':
        lkl_method, bin_method = 'dolphin', 'bb'
    # Observed cluster's histogram and bin edges for each dimension.
    bin_edges, cl_histo = obs_clust_prepare.main(
        cl_reg_fit, lkl_method, bin_method)[:2]

    # Histogram of the synthetic cluster, using the bin edges calculated
    # with the observed cluster.
    hess_diag = np.array([])
    if synth_clust:
        synth_phot = synth_clust[0][0]
        if synth_phot:
            syn_histo = np.histogramdd(synth_phot, bins=bin_edges)[0]
            hess_nd = cl_histo - syn_histo
            # TODO this uses the first two defined photometric dimensions.
            hess_diag = hess_nd.reshape(hess_nd.shape[:2] + (-1,)).sum(axis=-1)

    if not hess_diag.size:
        print("  WARNING: the synthetic cluster is empty.")

    hess_data = {'hess_diag': hess_diag, 'hess_edges': bin_edges}

    return hess_data

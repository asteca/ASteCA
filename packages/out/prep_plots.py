
from ..inp import input_params as g
from ..math_f import exp_function
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


def coord_syst():
    '''
    Define system of coordinates used.
    '''
    px_deg = g.gd_params[-1]
    coord_lst = ['px', 'x', 'y'] if px_deg == 'px' else ['deg', 'ra', 'dec']
    return coord_lst


def frame_zoomed(x_min, x_max, y_min, y_max, center_cl, clust_rad):
    '''
    If possible, define zoomed frame.
    '''
    x_zmin, x_zmax = max(x_min, (center_cl[0] - 1.5 * clust_rad)), \
        min(x_max, (center_cl[0] + 1.5 * clust_rad))
    y_zmin, y_zmax = max(y_min, (center_cl[1] - 1.5 * clust_rad)), \
        min(y_max, (center_cl[1] + 1.5 * clust_rad))
    # Prevent axis stretching.
    if (x_zmax - x_zmin) != (y_zmax - y_zmin):
        lst = [(x_zmax - x_zmin), (y_zmax - y_zmin)]
        val, idx = min((val, idx) for (idx, val) in enumerate(lst))
        if idx == 0:
            x_zmax = x_zmin + lst[1]
        else:
            y_zmax = y_zmin + lst[0]

    return x_zmin, x_zmax, y_zmin, y_zmax


def ax_names():
    '''
    Define names for photometric diagram axes.
    '''
    # y_axis == 0 indicates that the y axis is a magnitude.
    y_axis = 0
    # Create photometric axis names.
    y_ax, x_ax0, m_ord = g.axes_params[0:3]
    if m_ord == 21:
        x_ax = '(' + x_ax0 + '-' + y_ax + ')'
    elif m_ord == 12:
        x_ax = '(' + y_ax + '-' + x_ax0 + ')'

    return x_ax, y_ax, x_ax0, y_axis


def ax_data(mag_data, col_data):
    '''
    Unpack coordinates and photometric data.
    '''
    # x_data, y_data = id_coords[1:]
    phot_x = col_data
    phot_y = mag_data
    return phot_x, phot_y


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
    cs = plt.contour(x, y, np.reshape(k_pos, x.shape), 30)
    # Extract (x,y) points delimitating each line.
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
    x_v, y_v = kde_limits(phot_x, phot_y)

    # Define diagram limits.
    x_min_cmd, x_max_cmd = min(x_v) - 1.25, max(x_v) + 1.25
    y_min_cmd = max(y_v) + 1.25
    # If photometric axis y is a magnitude, make sure the brightest star
    # is always plotted.
    if y_axis == 0:
        y_max_cmd = min(phot_y) - 1.
    else:
        y_max_cmd = min(y_v) - 1.

    return x_max_cmd, x_min_cmd, y_min_cmd, y_max_cmd


def star_size(mag_data):
    '''
    Convert magnitudes into intensities and define sizes of stars in
    finding chart.
    '''
    # Scale factor.
    factor = 500. * (1 - 1 / (1 + 150 / len(mag_data) ** 0.85))
    return 0.1 + factor * 10 ** ((np.array(mag_data) - min(mag_data)) / -2.5)


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


def separate_stars(x_data, y_data, mag_data, x_zmin, x_zmax, y_zmin, y_zmax,
                   stars_out_rjct, field_regions):
    '''
    Separate stars in lists.
    '''
    # Separate stars in zoomed frame.
    x_data_z, y_data_z, mag_data_z = [], [], []
    for i, st_x in enumerate(x_data):
        st_y, st_mag = y_data[i], mag_data[i]
        if x_zmin <= st_x <= x_zmax and y_zmin <= st_y <= y_zmax:
            x_data_z.append(st_x)
            y_data_z.append(st_y)
            mag_data_z.append(st_mag)

    # Generate list with *all* rejected stars outside of the cluster region.
    stars_f_rjct = [[], []]
    for star in stars_out_rjct:
        stars_f_rjct[0].append(star[5])
        stars_f_rjct[1].append(star[3])

    # Generate list with stars within a defined field region.
    stars_f_acpt = [[], []]
    if field_regions:
        for fr in field_regions:
            for star in fr:
                stars_f_acpt[0].append(star[5])
                stars_f_acpt[1].append(star[3])

    return x_data_z, y_data_z, mag_data_z, stars_f_rjct, stars_f_acpt


def da_plots(center_cl, clust_rad, stars_out, x_zmin, x_zmax, y_zmin, y_zmax,
             x_max_cmd, col_data, err_lst, red_return):
    '''
    Generate parameters for the finding chart and the photometric diagram
    plotted with the MPs assigned by the DA.
    '''

    red_memb_fit, red_memb_no_fit = red_return[:2]

    # Get extreme values for colorbar.
    lst_comb = red_memb_fit + red_memb_no_fit
    v_min_mp, v_max_mp = round(min(zip(*lst_comb)[-1]), 2), \
        round(max(zip(*lst_comb)[-1]), 2)

    # Decides if colorbar should be plotted.
    plot_colorbar = False
    if v_min_mp != v_max_mp:
        plot_colorbar = True

    # Finding chart.
    # Arrange stars used in the best fit process.
    m_p_m_temp = [[], [], []]
    for star in red_memb_fit:
        m_p_m_temp[0].append(star[1])
        m_p_m_temp[1].append(star[2])
        m_p_m_temp[2].append(star[7])
    # Create new list with inverted values so higher prob stars are on top.
    chart_fit_inv = [i[::-1] for i in m_p_m_temp]
    # Arrange stars *not* used in the best fit process.
    m_p_m_temp = [[], [], []]
    for star in red_memb_no_fit:
        m_p_m_temp[0].append(star[1])
        m_p_m_temp[1].append(star[2])
        m_p_m_temp[2].append(star[7])
    # Create new list with inverted values so higher prob stars are on top.
    chart_no_fit_inv = [i[::-1] for i in m_p_m_temp]

    # Separate stars outside the cluster's radius.
    out_clust_rad = [[], []]
    for star in stars_out:
        if x_zmin <= star[1] <= x_zmax and y_zmin <= star[2] <= y_zmax:
            dist = np.sqrt((center_cl[0] - star[1]) ** 2 +
                           (center_cl[1] - star[2]) ** 2)
            # Only plot stars outside the cluster's radius.
            if dist >= clust_rad:
                out_clust_rad[0].append(star[1])
                out_clust_rad[1].append(star[2])

    # Photometric diagram.
    # Arrange stars used in the best fit process.
    m_p_m_temp = [[], [], []]
    for star in red_memb_fit:
        m_p_m_temp[0].append(star[5])
        m_p_m_temp[1].append(star[3])
        m_p_m_temp[2].append(star[7])
    # Create new list with inverted values so higher prob stars are on top.
    diag_fit_inv = [i[::-1] for i in m_p_m_temp]
    # Arrange stars *not* used in the best fit process.
    m_p_m_temp = [[], [], []]
    for star in red_memb_no_fit:
        m_p_m_temp[0].append(star[5])
        m_p_m_temp[1].append(star[3])
        m_p_m_temp[2].append(star[7])
    # Create new list with inverted values so higher prob stars are on top.
    diag_no_fit_inv = [i[::-1] for i in m_p_m_temp]

    # For plotting error bars in photom diag.
    x_val, mag_y, x_err, y_err = [], [], [], []
    if zip(*lst_comb)[3]:
        mag_y = np.arange(int(min(zip(*lst_comb)[3]) + 0.5),
                          int(max(zip(*lst_comb)[3]) + 0.5) + 0.1)
        x_val = [min(x_max_cmd, max(col_data) + 0.2) - 0.4] * len(mag_y)
        # Read average fitted values for exponential error fit.
        popt_mag, popt_col1 = err_lst[:2]
        x_err = exp_function.exp_3p(mag_y, *popt_col1)
        y_err = exp_function.exp_3p(mag_y, *popt_mag)
    err_bar = [x_val, mag_y, x_err, y_err]

    return v_min_mp, v_max_mp, plot_colorbar, chart_fit_inv, \
        chart_no_fit_inv, out_clust_rad, diag_fit_inv, diag_no_fit_inv, err_bar


def param_ranges(ip_list):
    '''
    Set parameter ranges used by GA plots.
    '''
    min_max_p = []
    for param in ip_list[1]:
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

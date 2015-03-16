# -*- coding: utf-8 -*-
"""
Created on Thu Dic 18 12:00:00 2014

@author: gabriel
"""

from .._in import get_in_params as g
import numpy as np
import math


def reject_outliers(data, m=6.5):
    '''
    Reject outliers from array.
    http://stackoverflow.com/a/16562028/1391441
    '''
    d = np.abs(data - np.median(data))
    mdev = np.median(d)
    s = d / mdev if mdev else np.asarray([0. for _ in d])
    return data[s < m]


def frame_max_min(x_data, y_data):
    '''
    Get max and min values in x,y coordinates.
    '''
    x_min, x_max = min(x_data), max(x_data)
    y_min, y_max = min(y_data), max(y_data)
    return x_min, x_max, y_min, y_max


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
    y_axis = 0
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
    #x_data, y_data = id_coords[1:]
    phot_x = col_data
    phot_y = mag_data
    return phot_x, phot_y


def diag_limits(y_axis, phot_x, phot_y):
    '''
    Define plot limits for *all* photometric diagrams.
    '''
    phot_x_s, phot_y_s = reject_outliers(phot_x), reject_outliers(phot_y)
    x_max_cmd, x_min_cmd = max(phot_x_s) + 0.5, min(phot_x_s) - 0.5
    y_min_cmd, y_max_cmd = max(phot_y_s) + 0.5, min(phot_y_s) - 0.5
    # If photometric axis y is a magnitude, make sure the brighter stars
    # are plotted.
    y_max_cmd = (min(phot_y) - 1.) if y_axis == 0 else y_max_cmd

    return x_max_cmd, x_min_cmd, y_min_cmd, y_max_cmd


def star_size(mag_data):
    '''
    Convert magnitudes into intensities and define sizes of stars in
    finding chart.
    '''
    # Scale factor.
    factor = 500. * (1 - 1 / (1 + 150 / len(mag_data) ** 0.85))
    return 0.1 + factor * 10 ** ((np.array(mag_data) - min(mag_data)) / -2.5)


def error_vals(isoch_fit_errors):
    '''
    Set errors to zero if bootstrap was not used, for plotting purposes.
    '''
    N_b = g.bf_params[-1]
    if N_b >= 2:
        p_errs = isoch_fit_errors
    else:
        p_errs = [0.] * len(isoch_fit_errors)
    return p_errs


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


def BestTick(minv, maxv):
    '''
    Find optimal number and length of ticks for a given fixed maximum
    number of characters in the axis.
    '''

    # Define maximum number of characters in axis.
    max_char = 30
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

# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 20:55:57 2013

@author: gabriel
"""

import numpy as np


def calc_integ_mag(st_reg):
    '''
    Calculate integrated magnitude up to a certain maximum magnitude value.
    '''

    # Define magnitude range.
    mag_range = np.linspace(min(st_reg), max(st_reg), 100)

    # Define final output lists.
    reg_mag = [[], []]

    # Sort magnitude list.
    sort_lis = sorted(st_reg)

    int_mag_val = 0.
    for mag_limit in mag_range:

        # Calculate integrated magnitude up to this mag value.
        trim_lis = [i for i in sort_lis if i <= mag_limit]

        if len(trim_lis) == 0:
            int_mag_val = mag_limit
        else:
            energ_sum = 0.
            for mag_i in trim_lis:
                energ_sum = energ_sum + 10 ** (mag_i / -2.5)

            int_mag_val = -2.5 * np.log10(energ_sum)

        reg_mag[0].append(mag_limit)
        reg_mag[1].append(int_mag_val)

    return reg_mag


def integ_mag(center_cl, clust_rad, cluster_region, field_region,
    flag_area_stronger):
    '''
    Obtain integrated magnitude using all stars inside the cluster's radius for
    several limits in magnitude.
    '''

    # Only use stars inside cluster's radius.
    cl_region_r = [[], []]
    for star in cluster_region:
        dist = np.sqrt((center_cl[0] - star[1]) ** 2 +
        (center_cl[1] - star[2]) ** 2)
        if dist <= clust_rad:
            # Append magnitude.
            cl_region_r[0].append(star[3])
            # Append color.
            cl_region_r[1].append(star[5] + star[3])

    cl_reg_mag = calc_integ_mag(cl_region_r[0])
    cl_reg_col = calc_integ_mag(cl_region_r[1])

    if flag_area_stronger is not True:
        # Run for every field region defined.
        fl_reg_m_1, fl_reg_c_1 = [[], []], [[], []]
        for f_reg in field_region:
            # For magnitude values.
            fl_reg_m_0 = calc_integ_mag(zip(*f_reg)[3])
            fl_reg_m_1[0].extend(fl_reg_m_0[0])
            fl_reg_m_1[1].extend(fl_reg_m_0[1])
            # For color values.
            fl_reg_c_0 = calc_integ_mag(np.asarray(zip(*f_reg)[5]) +
                np.asarray(zip(*f_reg)[3]))
            fl_reg_c_1[0].extend(fl_reg_c_0[0])
            fl_reg_c_1[1].extend(fl_reg_c_0[1])

        # Sort arrays reversing the integrated magnitudes/colors.
        fl_reg_m_2, fl_reg_c_2 = [[], []], [[], []]
        fl_reg_m_2[0] = np.sort(fl_reg_m_1[0])
        fl_reg_m_2[1] = np.sort(fl_reg_m_1[1])[::-1]
        fl_reg_c_2[0] = np.sort(fl_reg_c_1[0])
        fl_reg_c_2[1] = np.sort(fl_reg_c_1[1])[::-1]

        # Interpolate all curves to obtain final field integrated and color
        # magnitude.
        fl_reg_mag, fl_reg_col = [[], []], [[], []]
        fl_reg_mag[0] = np.linspace(min(fl_reg_m_2[0]), max(fl_reg_m_2[0]), 200)
        fl_reg_mag[1] = np.interp(fl_reg_mag[0], fl_reg_m_2[0], fl_reg_m_2[1])
        fl_reg_col[0] = np.linspace(min(fl_reg_c_2[0]), max(fl_reg_c_2[0]), 200)
        fl_reg_col[1] = np.interp(fl_reg_col[0], fl_reg_c_2[0], fl_reg_c_2[1])

        # Obtain integrated magnitude of clean cluster region, ie: substracting
        # the field contribution.
        if min(fl_reg_mag[1]) >= min(cl_reg_mag[1]):
            integ_mag = -2.5 * np.log10(1 - 10 ** ((min(fl_reg_mag[1]) -
            min(cl_reg_mag[1])) / -2.5)) + min(cl_reg_mag[1])
        else:
            # If the field is brighter than the cluster.
            integ_mag = min(cl_reg_mag[1])

        # Obtain integrated second magnitude of clean cluster region.
        if min(fl_reg_col[1]) >= min(cl_reg_col[1]):
            integ_col = -2.5 * np.log10(1 - 10 ** ((min(fl_reg_col[1]) -
            min(cl_reg_col[1])) / -2.5)) + min(cl_reg_col[1])
        else:
            # If the field is brighter than the cluster.
            integ_col = min(cl_reg_col[1])

    else:
        print '  WARNING: no field regions defined. Integrated magnitude'
        print '  is not cleaned from field star contamination.'
        # Pass dummy lists.
        fl_reg_mag, fl_reg_col = [np.array([]), np.array([])], \
        [np.array([]), np.array([])]
        integ_mag = min(cl_reg_mag[1])
        integ_col = min(cl_reg_col[1])

    integr_return = [cl_reg_mag, fl_reg_mag, integ_mag, cl_reg_col,
        fl_reg_col, integ_col]

    return integr_return

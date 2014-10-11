# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 20:55:57 2013

@author: gabriel
"""

import numpy as np
from scipy.interpolate import spline
import get_in_params as g


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


def field_reg_integ_mag_curve(fl_reg_m):
    '''
    Return smooth averaged curve for integrated magnitude.
    '''

    # Average all field regions for the integrated magnitude.
    fl_reg_m1 = np.average(fl_reg_m, axis=0)

    # Smooth curve for the averaged integrated magnitude.
    fl_reg_mag = []
    # Magnitude values.
    fl_reg_mag.append(np.linspace(min(fl_reg_m1[0]), max(fl_reg_m1[0]), 100))
    # Integ magnitude values.
    fl_reg_mag.append(spline(fl_reg_m1[0], fl_reg_m1[1], fl_reg_mag[0]))

    return fl_reg_mag


def integ_mag(cl_region, field_region, flag_area_stronger):
    '''
    Obtain integrated magnitude using all stars inside the cluster's radius for
    several limits in magnitude.
    '''

    if g.im_flag:

        # This variable tells me how the color is created, if the first
        # magnitude is substracted from the second one or the other way around.
        m_ord = g.axes_params[2]
        # Check how the second magnitude whould be formed.
        sig = 1. if m_ord == 21 else -1.

        # Only use stars inside cluster's radius.
        cl_region_r = [[], []]
        for star in cl_region:
            # Append first magnitude.
            cl_region_r[0].append(star[3])
            # Append second magnitude.
            cl_region_r[1].append(star[3] + sig * star[5])

        cl_reg_mag1 = calc_integ_mag(cl_region_r[0])
        cl_reg_mag2 = calc_integ_mag(cl_region_r[1])

        if flag_area_stronger is not True:

            # Run for every field region defined.
            fl_reg_m = [[], []]
            for f_reg in field_region:
                # First magnitude values.
                fl_reg_m[0].append(calc_integ_mag(zip(*f_reg)[3]))
                # Second magnitude values.
                fl_reg_m[1].append(calc_integ_mag(np.asarray(zip(*f_reg)[3]) +
                    sig * np.asarray(zip(*f_reg)[5])))

            fl_reg_mag1 = field_reg_integ_mag_curve(fl_reg_m[0])
            fl_reg_mag2 = field_reg_integ_mag_curve(fl_reg_m[1])

            # Obtain integrated magnitude of clean cluster region, ie:
            # substracting the field contribution.
            if min(fl_reg_mag1[1]) >= min(cl_reg_mag1[1]):
                integ_mag1 = -2.5 * np.log10(1 - 10 ** ((min(fl_reg_mag1[1]) -
                min(cl_reg_mag1[1])) / -2.5)) + min(cl_reg_mag1[1])
            else:
                # If the field is brighter than the cluster.
                integ_mag1 = min(cl_reg_mag1[1])

            # Obtain integrated second magnitude of clean cluster region.
            if min(fl_reg_mag2[1]) >= min(cl_reg_mag2[1]):
                integ_mag2 = -2.5 * np.log10(1 - 10 ** ((min(fl_reg_mag2[1]) -
                min(cl_reg_mag2[1])) / -2.5)) + min(cl_reg_mag2[1])
            else:
                # If the field is brighter than the cluster.
                integ_mag2 = min(cl_reg_mag2[1])

        else:
            print '  WARNING: no field regions defined. Integrated magnitude'
            print '  is not cleaned from field star contamination.'
            # Pass dummy lists.
            fl_reg_mag1, fl_reg_mag2 = [np.array([]), np.array([])], \
            [np.array([]), np.array([])]
            integ_mag1 = min(cl_reg_mag1[1])
            integ_mag2 = min(cl_reg_mag2[1])

        integr_return = [cl_reg_mag1, fl_reg_mag1, integ_mag1, cl_reg_mag2,
            fl_reg_mag2, integ_mag2]

        int_col = sig * (integ_mag2 - integ_mag1)
        print 'Integrated color magnitude distribution obtained (%0.2f).' % \
            int_col

    else:
        print 'Skipping integrated magnitudes function.'
        integr_return = []

    return integr_return

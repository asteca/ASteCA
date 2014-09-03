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


def movingaverage(interval, window_size):
    window = np.ones(int(window_size)) / float(window_size)
    return np.convolve(interval, window, 'same')


def integ_mag(center_cl, clust_rad, cluster_region, field_region, axes_params,
    flag_area_stronger):
    '''
    Obtain integrated magnitude using all stars inside the cluster's radius for
    several limits in magnitude.
    '''

    # This variable tells me how the color is created, if the first magnitude
    # is substracted from the second one or the other way around.
    m_ord = axes_params[2]
    # Check how the second magnitude whould be formed.
    if m_ord == 21:
        sig = 1.
    elif m_ord == 12:
        sig = -1.

    # Only use stars inside cluster's radius.
    cl_region_r = [[], []]
    for star in cluster_region:
        dist = np.sqrt((center_cl[0] - star[1]) ** 2 +
        (center_cl[1] - star[2]) ** 2)
        if dist <= clust_rad:
            # Append first magnitude.
            cl_region_r[0].append(star[3])
            # Append second magnitude.
            cl_region_r[1].append(star[3] + sig * star[5])

    cl_reg_mag1 = calc_integ_mag(cl_region_r[0])
    cl_reg_mag2 = calc_integ_mag(cl_region_r[1])

    if flag_area_stronger is not True:
        # Run for every field region defined.
        fl_reg_m11, fl_reg_m21 = [[], []], [[], []]
        for f_reg in field_region:

            # First magnitude values.
            fl_reg_m10 = calc_integ_mag(zip(*f_reg)[3])
            fl_reg_m11[0].extend(fl_reg_m10[0])
            fl_reg_m11[1].extend(fl_reg_m10[1])
            # Second magnitude values.
            fl_reg_m20 = calc_integ_mag(np.asarray(zip(*f_reg)[3]) +
                sig * np.asarray(zip(*f_reg)[5]))
            fl_reg_m21[0].extend(fl_reg_m20[0])
            fl_reg_m21[1].extend(fl_reg_m20[1])

        # Sort arrays reversing the integrated magnitudes.
        fl_reg_m12, fl_reg_m22 = [[], []], [[], []]
        fl_reg_m12[0] = np.sort(fl_reg_m11[0])
        fl_reg_m12[1] = np.sort(fl_reg_m11[1])[::-1]
        fl_reg_m22[0] = np.sort(fl_reg_m21[0])
        fl_reg_m22[1] = np.sort(fl_reg_m21[1])[::-1]

        # Interpolate all curves to obtain final field integrated and color
        # magnitude.
        fl_reg_mag1, fl_reg_mag2 = [[], []], [[], []]
        fl_reg_mag1[0] = np.linspace(min(fl_reg_m12[0]), max(fl_reg_m12[0]),
            20)
        fl_reg_mag1[1] = np.interp(fl_reg_mag1[0], fl_reg_m12[0], fl_reg_m12[1])
        fl_reg_mag2[0] = np.linspace(min(fl_reg_m22[0]), max(fl_reg_m22[0]),
            20)
        fl_reg_mag2[1] = np.interp(fl_reg_mag2[0], fl_reg_m22[0], fl_reg_m22[1])

        smoothed = movingaverage(fl_reg_m21[0], 50)
        import matplotlib.pyplot as plt
        print fl_reg_m21
        #for i, f_reg in enumerate(field_region):
            #print i, len(zip(*f_reg)[3])
            #fl_reg_m20 = calc_integ_mag(np.asarray(zip(*f_reg)[3]) +
                #sig * np.asarray(zip(*f_reg)[5]))
            #plt.scatter(fl_reg_m20[0], fl_reg_m20[1], label=i)
        #plt.plot(fl_reg_mag2[0], fl_reg_mag2[1], label='interp', c='r')
        plt.scatter(fl_reg_m21[0], fl_reg_m21[1], label='mov_avr', c='b')
        #plt.gca().invert_yaxis()
        #plt.legend()
        plt.show()
        raw_input()

        # Obtain integrated magnitude of clean cluster region, ie: substracting
        # the field contribution.
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

    return integr_return

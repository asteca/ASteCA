"""
@author: gabriel
"""

import numpy as np
import matplotlib.pyplot as plt
from display_rad import disp_rad as d_r
from .._in import get_in_params as g


def main_rad_algor(rdp_params, field_dens, bin_width, coord):
    '''
    This function holds the main algorithm that returns a radius value.
    '''
    # Unpack values.
    radii, rdp_points = rdp_params[:2]

    # Difference between max RDP density value and the field density value.
    delta_total = (max(rdp_points) - field_dens)
    # If the difference between the max density value and the field density is
    # less than 3 times the value of the field density, raise a flag.
    flag_delta_total = False
    if delta_total < 3 * field_dens:
        flag_delta_total = True

    # Find maximum density value and assume this is the central density.
    # Do not use previous values.
    #max_dens_ind = np.argmax(rdp_points)
    #rdp_points_c = rdp_points[max_dens_ind:]
    #radii_c = radii[max_dens_ind:]

    # Interpolate extra RDP points between those calculated.
    N = 100
    t, xp = np.linspace(0, 1, N), np.linspace(0, 1, len(rdp_points))
    # Store RDP interpolated values.
    rdp_points_i = np.interp(t, xp, rdp_points)
    radii_i = np.interp(t, xp, radii)
    p_width = (max(radii) - min(radii)) / N

    #
    # Core algorithm for radius finding.
    #
    # Store each point in the RDP as a bar of given area, given by:
    # (height of the point - field density value) x (bin width).
    rad_found = []
    for smooth_param in np.arange(0.01, 0.091, 0.01):

        rdp_perc = []
        for dens_rad in rdp_points_i:
            if dens_rad < (field_dens + smooth_param * delta_total):
                rdp_perc.append(0.)
            else:
                rdp_perc.append(max((dens_rad - field_dens) * p_width, 0.))

        # Sum all areas to obtain total RDP area.
        tot_rdp_dens = sum(rdp_perc)

        # Convert each area to the given percentage it represents in the total
        # RDP area, and add each previous percentage value. This way we generate
        # a list containing the *cummulative* RDP percentages.
        perc_dens = [rdp_perc[0] / tot_rdp_dens]
        for i, perc in enumerate(rdp_perc[1:]):
            d_perc = perc / tot_rdp_dens
            perc_dens.append(perc_dens[i] + d_perc)

        from scipy.optimize import curve_fit
        def fit_func(x, a, b):
            return 1 - 1 / (a * (x ** b) + 1)
        ab_p, ab_e = curve_fit(fit_func, np.asarray(radii_i) / max(radii_i),
            np.asarray(perc_dens))

        #with open('/media/rest/github/asteca/output/cc.dat', "a") as f_out:
            #f_out.write(''' {:>15.2f}   {:>15.2f}'''.format(*ab_p))
            #f_out.write('\n')

        import bisect
        cc_lim = [10, 50, 100, 200, 400, float("inf")]
        perc_lim = [0.9, 0.92, 0.94, 0.96, 0.97, 1.]
        indx_lim = bisect.bisect(cc_lim, ab_p[0])
        print ab_p[0], indx_lim, perc_lim[indx_lim]

        #rad_found = []
        #for p in np.arange(0.89, perc_lim[indx_lim], 0.01):
        #for p in np.arange(0.9, 0.99, 0.01):
        #p = 1.0 - smooth_param
        p = perc_lim[indx_lim]
        i_9x = min(range(len(perc_dens)), key=lambda i: abs(perc_dens[i] - p))
        slope = []
        for i, r in enumerate(radii_i[:i_9x]):
            slope.append((perc_dens[i] - p) / (r - radii_i[i_9x]))
        rad_found.append(radii_i[slope.index(min(slope))])
            #print p, radii_i[slope.index(min(slope))]
            #plt.scatter(radii_i[:i_9x], slope)
        #plt.show()

    # Check mode.
    mode_r = g.cr_params[0]
    if mode_r not in {'auto', 'manual'}:
        print "  WARNING: CR mode is not valid. Default to 'auto'."
        mode_r = 'auto'

    # Raise a flag if less than two radius values were found under all the
    # n_deltas values.
    flag_delta, flag_not_stable = False, False
    #if len(index_rad) < 3:
    if len(rad_found) < 3:
        flag_delta = True

    #rad_found = []
    #for ind in index_rad:
        ## Use the stored indexes to obtain the actual radii values.
        #rad_found.append(radii_i[ind])

    # If at least one radius value was found.
    if rad_found:

        print rad_found
        raw_input()
        from scipy import stats
        import warnings
        clust_rad, _std = np.mean(rad_found), np.std(rad_found)
        # Catch warning if stats fails to obtain confidence interval.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            conf_int = stats.norm.interval(0.68, loc=clust_rad, scale=_std /
                np.sqrt(len(rad_found)))
        # If stats returns a confidence interval with a NaN, discard it.
        if np.any(np.isnan(conf_int)):
            conf_int = [0., 0.]
        conf_dif = abs(clust_rad - max(conf_int))
        e_rad = max(conf_dif, bin_width)
        # Prevent too small radius by fixing the minimum value to the second
        # RDP point.
        if clust_rad < radii[1]:
            clust_rad = radii[1]
    else:
        flag_not_stable = True
        # No radius value found. Assign radius value as the middle element
        # in the full radii list.
        clust_rad, e_rad = radii_i[int(len(radii_i) / 2.)], 0.
        print '  WARNING: no radius found, setting value to: {:g} {}'.format(
            clust_rad, coord)

    return clust_rad, e_rad, flag_delta_total, flag_not_stable, flag_delta


def get_clust_rad(phot_data, field_dens, center_params, rdp_params,
    semi_return, bin_width, coord_lst):
    """
    Obtain the value for the cluster's radius by counting the number of points
    that fall within a given interval of the field density or lower. If this
    number is equal to a minimum fixed number of points 'n_left', then assign
    the radius as the point closest to the field density value among those
    first n_left points counting from the first one that fell below the
    (field dens + delta limit) range.
    Iterate increasing the interval around the field density and finally
    average all the radius values found for each interval.
    """

    coord = coord_lst[0]
    # Call function that holds the radius finding algorithm.
    clust_rad, e_rad, flag_delta_total, flag_not_stable, flag_delta = \
    main_rad_algor(rdp_params, field_dens, bin_width, coord)

    # Check if semi or manual mode are set.
    flag_radius_manual = False
    if g.mode == 'auto':
        print 'Auto radius found: {:g} {}.'.format(clust_rad, coord)

    elif g.mode == 'semi':
        # Unpack semi values.
        cent_cl_semi, cl_rad_semi, cent_flag_semi, rad_flag_semi, \
        err_flag_semi = semi_return

        if rad_flag_semi == 1:
            # Update values.
            clust_rad, e_rad = cl_rad_semi, 0.
            print 'Semi radius set: {:g} {}.'.format(clust_rad, coord)
        else:
            print 'Auto radius found: {:g} {}.'.format(clust_rad, coord)

    # If Manual mode is set, display radius and ask the user to accept it or
    # input new one.
    elif g.mode == 'manual':

        print 'Radius found: {:g} {}.'.format(clust_rad, coord)
        d_r(phot_data, bin_width, center_params, clust_rad, field_dens,
            rdp_params)
        plt.show()

        wrong_answer = True
        while wrong_answer:
            answer_rad = raw_input('Accept radius (otherwise input new one \
manually)? (y/n) ')
            if answer_rad == 'y':
                print 'Value accepted.'
                wrong_answer = False
            elif answer_rad == 'n':
                clust_rad_m = float(raw_input('Input new radius value (in \
px): '))
                # Update radius value.
                clust_rad = clust_rad_m
                wrong_answer = False
                flag_radius_manual = True
            else:
                print 'Wrong input. Try again.\n'

    radius_params = [clust_rad, e_rad, flag_delta_total, flag_not_stable,
        flag_delta, flag_radius_manual]

    return radius_params
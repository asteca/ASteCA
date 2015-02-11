"""
@author: gabriel
"""

import numpy as np
import matplotlib.pyplot as plt
from display_rad import disp_rad as d_r
from scipy.optimize import curve_fit
from scipy import stats
from scipy import integrate
#from scipy.optimize import fsolve
import warnings
from .._in import get_in_params as g


def fit_func(x, a, b):
    '''
    Function to fit the cummulative RDP.
    '''
    return 1 - 1 / (a * (x ** b) + 1)


def fit_func_inv(y, a, b):
    '''
    Inverted fit_func.
    '''
    return ((1 / a) * ((1 / (1 - y)) - 1)) ** (1 / b)


def fit_func_diff(x, a, b):
    '''
    Derivative of the cummulative RDP.
    '''
    return (a * b * x ** (b - 1)) / (a * x ** b + 1) ** 2


def main_rad_algor(rdp_params, field_dens, bin_width, coord, semi_return):
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
    N = 1000
    t, xp = np.linspace(0, 1, N), np.linspace(0, 1, len(rdp_points))
    # Store interpolated RDP values.
    rdp_points_i = np.interp(t, xp, rdp_points)
    # Store interpolated normalized radii values.
    radii_i = np.interp(t, xp, radii) / max(radii)
    p_width = (max(radii) - min(radii)) / N

    #
    # Core algorithm for radius finding.
    #
    #
    rad_found = []

    try:
        # Store each point in the RDP as a bar of given area, given by:
        # (height of the point - field density value) x delta_radius / N.
        rdp_perc = []
        for dens_rad in rdp_points_i:
            #if dens_rad < (field_dens + 0.1 * delta_total):
                #rdp_perc.append(0.)
            #else:
                rdp_perc.append(max((dens_rad - field_dens) * p_width, 0.))

        # Sum all areas to obtain total RDP area.
        tot_rdp_dens = sum(rdp_perc)

        # Convert each area to the given percentage it represents in the
        # total RDP area, and add each previous percentage value. This way
        # we generate a list containing the *cummulative* RDP percentages.
        perc_dens = [rdp_perc[0] / tot_rdp_dens]
        for i, perc in enumerate(rdp_perc[1:]):
            d_perc = perc / tot_rdp_dens
            perc_dens.append(perc_dens[i] + d_perc)

        # Obtain best fit params for the RDP cumulative function.
        ab, ab_e = curve_fit(fit_func, np.asarray(radii_i),
            np.asarray(perc_dens))
        a, b = ab

        # Integrate the fit_func using the fitted a & b values.
        integ_tot = integrate.quad(fit_func, 0., 1., (a, b))[0]

        # Obtain radii values for the following cummulative RDP values.
        rdp_i = [0.1, 0.25, 0.5, 0.75]
        r_i = []
        for i in rdp_i:
            r_i.append(fit_func_inv(i, a, b))

        # Obtain average slope in the range between these radii.
        diff_avrg = []
        for i in range(len(rdp_i) - 1):
            diff_ij = []
            for j in np.linspace(r_i[i], r_i[i + 1], 10):
                diff_ij.append(fit_func_diff(j, a, b))
            # Get average.
            diff_avrg.append(np.mean(diff_ij))
        #slope_1_2 = (rdp_1 - rdp_2) / (r_1 - r_2)
        print r_i
        print diff_avrg
        y = fit_func(radii_i, a, b)
        plt.plot(radii_i, y)
        plt.scatter(radii_i, perc_dens)
        plt.show()

        # The concentration of the cluster is defined loosely by the 'a'
        # parameter in the RDP cummulative function.
        # A low value of 'a' means a loose RDP and a high value means the
        # cummulative RDP is concentrated.
        diff_0205 = fit_func_diff(0.2, a, b) / fit_func_diff(0.5, a, b)
        if diff_0205 > 20.:  # Very concentrated cluster.
            dens_lim, diff_lim, int_lim = 0.10, 0.25, 0.99
        elif 15. < diff_0205 <= 20.:  # Somewhat concentrated cluster.
            dens_lim, diff_lim, int_lim = 0.10, 0.50, 0.95
        elif 10. < diff_0205 <= 15.:  # Somewhat concentrated cluster.
            dens_lim, diff_lim, int_lim = 0.15, 0.75, 0.9
        elif 5. < diff_0205 <= 10.:  # Not so concentrated cluster.
            dens_lim, diff_lim, int_lim = 0.20, 1.00, 0.8
        elif diff_0205 < 5.:  # Loose cluster.
            dens_lim, diff_lim, int_lim = 0.25, 1.25, 0.7

        # Iterate throught the RDP values, starting from the closest one
        # to the cluster's center.
        rad_lim = 50
        for i, r in enumerate(radii_i):

            # Evaluate the derivative of the fitted function in 'r'.
            r_diff = fit_func_diff(r, a, b)

            # Evaluate the integral of the fitted function in 'r'.
            r_inte = integrate.quad(fit_func, 0., r, (a, b))[0]

            # Evaluate the RDP in 'r'.
            r_dens = rdp_points_i[i]

            #print r, r * max(radii), abs(r_dens - field_dens), r_diff, r-r_inte
            if len(rad_found) > rad_lim:
                break
            # All conditions met.
            elif abs(r_dens - field_dens) <= dens_lim * delta_total and\
                r_diff < diff_lim and \
                (r - r_inte) >= int_lim * (1 - integ_tot):
                    rad_found.append(r * max(radii))
                    break
            # Conditions density and diff.
            elif abs(r_dens - field_dens) <= dens_lim * delta_total and \
                r_diff < diff_lim:
                rad_found.append(r * max(radii))
            # Conditions density and integral.
            elif abs(r_dens - field_dens) <= dens_lim * delta_total and \
                (r - r_inte) >= int_lim * (1 - integ_tot):
                rad_found.append(r * max(radii))
            # Conditions diff and integral.
            elif r_diff < diff_lim and \
                (r - r_inte) >= int_lim * (1 - integ_tot):
                rad_found.append(r * max(radii))

        #for c in np.arange(sl_rang[i_s][0], sl_rang[i_s][1], 0.01):
            ## Estimated aprox solution (initial guess for fsolve)
            #est = (ab[0] * c / ab[1]) ** (-1 / (ab[1] + 1))
            #with warnings.catch_warnings():
                #warnings.simplefilter("ignore")
                #rad_p = fsolve(fit_func_diff, est, (ab[0], ab[1], c))
            ##print rad_p, est

            ## Store radius value found.
            #rad_found.append(rad_p * max(radii))

        ##
        ##
        ## Integrate the fit_func using the fitted a & b values.
        #integ_tot = integrate.quad(fit_func, 0., 1., (ab[0], ab[1]))[0]

        #int_res = [[], []]
        #for lim in np.arange(0.05, 0.991, 0.01):
            ## Integrate the fit_func using the fitted a & b values.
            #integ_result = integrate.quad(fit_func, 0., lim,
                #(ab[0], ab[1]))[0] / integ_tot
            #int_res[0].append(lim)
            #int_res[1].append(fit_func(lim, ab[0], ab[1]) / integ_result)

            ## Find index of radius value closest to the integral value.
            #idx_integ = min(range(len(perc_dens)),
                #key=lambda i: abs(perc_dens[i] - integ_result))

            ## Store radius value found.
            #rad_found.append(radii_i[idx_integ] * max(radii))

        #rad_semi = semi_return[1] / max(radii)
        #per_in_rad_semi = fit_func(rad_semi, ab[0], ab[1])
        #rad_di = []
        #integ_tot = integrate.quad(fit_func, 0., 1., (ab[0], ab[1]))[0]
        #for p in [0.5, 0.9]:
            #r_p = ((1 / ab[0]) * ((1 / (1 - p)) - 1)) ** (1 / ab[1])
            ## Derivative and integral evaluated in several values.
            #rad_di.append(fit_func_diff(r_p, ab[0], ab[1]))
            ## Integral
            #rad_di.append(integrate.quad(fit_func, 0., r_p,
                #(ab[0], ab[1]))[0] / integ_tot)

        #rad_di = []
        #for p in [0.5, 0.9]:
            #r_p = ((1 / ab[0]) * ((1 / (1 - p)) - 1)) ** (1 / ab[1])
            ## Derivative and integral evaluated in several values.
            #rad_di.append(fit_func_diff(r_p, ab[0], ab[1]))

        #x = ((ab[0] / rad_di[0]) ** (1 / rad_di[1])) * ab[1]
        #slope = (0.97 - 0.9) / (0. - 400.)
        #for x_i in np.linspace((x - x * 0.1), (x + x * 0.1), 10):
            #r_p = max(slope * x_i + 0.97, 0.2)
            #rad = ((1 / ab[0]) * ((1 / (1 - r_p)) - 1)) ** (1 / ab[1]) * max(radii)
            #print x, x_i, r_p, rad

            #rad_found.append(rad)

    except:
        pass

        #import bisect
        #cc_lim = [10, 50, 100, 200, 400, float("inf")]
        #perc_lim = [0.9, 0.92, 0.94, 0.96, 0.97, 1.]
        #indx_lim = bisect.bisect(cc_lim, ab_p[0])
        #print ab_p[0], indx_lim, perc_lim[indx_lim]

        ##rad_found = []
        ##for p in np.arange(0.89, perc_lim[indx_lim], 0.01):
        ##for p in np.arange(0.9, 0.99, 0.01):
        ##p = 1.0 - smooth_param
        #p = perc_lim[indx_lim]
        #i_9x = min(range(len(perc_dens)), key=lambda i: abs(perc_dens[i] - p))
        #slope = []
        #for i, r in enumerate(radii_i[:i_9x]):
            #slope.append((perc_dens[i] - p) / (r - radii_i[i_9x]))
        #rad_found.append(radii_i[slope.index(min(slope))])
            ##print p, radii_i[slope.index(min(slope))]
            ##plt.scatter(radii_i[:i_9x], slope)
        ##plt.show()

    ## Check mode.
    #mode_r = g.cr_params[0]
    #if mode_r not in {'auto', 'manual'}:
        #print "  WARNING: CR mode is not valid. Default to 'auto'."
        #mode_r = 'auto'

    # Raise a flag if less than two radius values were found.
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

        p_r = []
        for r in [0.2, 0.5]:
            p_r.append(fit_func_diff(r, a, b))
        with open('/media/rest/github/asteca/output/cc.dat', "a") as f_out:
            f_out.write(''' {:>10.2f} {:>10.2f} {:>10.2f} {:>10.2f} {:>10.2f} {:>10.2f}'''.format(
                a, b, fit_func(0.1, a, b), p_r[0] / p_r[1], semi_return[1], clust_rad))
            f_out.write('\n')

    else:
        flag_not_stable = True
        # No radius value found. Assign radius value as the middle element
        # in the full radii list.
        clust_rad, e_rad = radii_i[int(len(radii_i) / 2.)] * max(radii), 0.
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
    main_rad_algor(rdp_params, field_dens, bin_width, coord, semi_return)

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
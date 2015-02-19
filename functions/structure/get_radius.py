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
    rad_found = []
    #try:
    # Store each point in the RDP as a bar of given area, given by:
    # (height of the point - field density value) x delta_radius / N.
    rdp_perc = []
    for dens_rad in rdp_points_i:
        #if dens_rad < (field_dens + 0.1 * delta_total):
            #rdp_perc.append(0.)
        #else:
            rdp_perc.append(max((dens_rad - field_dens) * p_width, 0.))

    max_r = np.argmax(rdp_points_i)
    ##f_lim = field_dens / max(rdp_points)
    #area = 0
    #for i, r in enumerate(radii_i[max_r:]):
        #area_prev = area
        #area = area + rdp_perc[max_r + i]
        #if abs(area_prev / area - 0.9) < 0.01:
            ## Evaluate the RDP in 'r'.
            #r_dens = rdp_points_i[max_r + i]
            #r_lim = (r_dens - field_dens) / delta_total
            #print area_prev / area, r_dens, r_lim
            #break
    #raw_input()

    # Sum all areas to obtain total RDP area.
    tot_rdp_dens = sum(rdp_perc)

    for perc in np.arange(0.1, 0.2, 0.01):
        area, points_added = 0., 0
        for i, r in enumerate(radii_i[max_r:]):
            #area_prev = area
            area = area + rdp_perc[max_r + i]
            # Evaluate the RDP in 'r'.
            r_dens = rdp_points_i[max_r + i]
            if area != 0:
                if (r_dens - field_dens) <= perc * delta_total and \
                area > (0.5 + perc) * tot_rdp_dens:
                    print '-->', points_added, r * max(radii), perc
                    rad_found.append(r * max(radii))
                    points_added += 1
            if points_added >= 10:
                # Once a minimum of 10 points have been found that satisfy
                # the condition, break out of this inner for loop.
                break

    clust_rad = np.mean(rad_found)
    with open('/media/rest/github/asteca/output/cc.dat', "a") as f_out:
        f_out.write(''' {:>10.2f}'''.format(clust_rad))
    #with open('/media/rest/github/asteca/output/cc.dat', "a") as f_out:
        f_out.write('\n')
    #print rad_found, np.mean(rad_found)
    #raw_input()

        ## Sum all areas to obtain total RDP area.
        #tot_rdp_dens = sum(rdp_perc)

        ## Convert each area to the given percentage it represents in the
        ## total RDP area, and add each previous percentage value. This way
        ## we generate a list containing the *cummulative* RDP percentages.
        #perc_dens = [rdp_perc[0] / tot_rdp_dens]
        #for i, perc in enumerate(rdp_perc[1:]):
            #d_perc = perc / tot_rdp_dens
            #perc_dens.append(perc_dens[i] + d_perc)

        ## Obtain best fit params for the RDP cumulative function.
        #ab, ab_e = curve_fit(fit_func, np.asarray(radii_i),
            #np.asarray(perc_dens))
        #a, b = ab

        ## Integrate the fit_func using the fitted a & b values.
        #integ_tot = integrate.quad(fit_func, 0., 1., (a, b))[0]

        ## Obtain average slope in the range between these radii.
        #diff_ij = []
        #for j in np.linspace(0.1, 0.2, 10):
            #diff_ij.append(fit_func_diff(j, a, b))
        ## Get average.
        #diff_avrg = np.mean(diff_ij)

        #r_08 = min(fit_func_inv(0.8, a, b), 0.9)
        #integ_08_1 = 100 * ((1 - r_08) -
            #integrate.quad(fit_func, r_08, 1., (a, b))[0])

        ## The concentration of the cluster is defined loosely by the 'diff_avrg'
        ## parameter. A low value means a loose RDP and a high value means the
        ## cummulative RDP is concentrated.

        #if diff_avrg > 20.:  # Very concentrated.
            #if integ_08_1 > 2:  # Very variable RDP tail.
                #dens_lim, diff_lim, int_lim = 0.1, 1.5, 0.75
            #elif integ_08_1 <= 2:  # More stable RDP tail.
                #dens_lim, diff_lim, int_lim = 0.1, 0.75, 0.95

        #elif 10. < diff_avrg <= 20.:  # Somewhat concentrated.
            #dens_lim, diff_lim, int_lim = 0.1, 1., 0.95

        #elif 5. < diff_avrg <= 10.:  # Not so concentrated.
            #if integ_08_1 > 4:  # Very variable RDP tail.
                #dens_lim, diff_lim, int_lim = 0.1, 1.5, 0.95
            #elif integ_08_1 <= 4:  # More stable RDP tail.
                #dens_lim, diff_lim, int_lim = 0.1, 0.5, 0.95

        #elif diff_avrg <= 5.:  # Loose.
            #dens_lim, diff_lim, int_lim = 0.1, 0.75, 0.95

        #par_lst = [[1.5, 0.75], [1.5, 0.9], [1.5, 0.95],
            #[1., 0.75], [1., 0.9], [1., 0.95],
            #[0.75, 0.75], [0.75, 0.9], [0.75, 0.95],
            #[0.5, 0.75], [0.5, 0.9], [0.5, 0.95],
            #[0.25, 0.75], [0.25, 0.9], [0.25, 0.95]]
        #for params in par_lst:
            #dens_lim = 0.1
            #diff_lim, int_lim = params
            #rad_found = []

            ## Iterate throught the RDP values, starting from the closest one
            ## to the cluster's center.
            #rad_lim = 10
            #for i, r in enumerate(radii_i):

                ## Evaluate the derivative of the fitted function in 'r'.
                #r_diff = fit_func_diff(r, a, b)

                ## Evaluate the integral of the fitted function in 'r'.
                #r_inte = integrate.quad(fit_func, 0., r, (a, b))[0]

                ## Evaluate the RDP in 'r'.
                #r_dens = rdp_points_i[i]

                #if len(rad_found) > rad_lim:
                    #break
                ## All conditions met.
                #elif abs(r_dens - field_dens) <= dens_lim * delta_total and\
                    #r_diff < diff_lim and \
                    #(r - r_inte) >= int_lim * (1 - integ_tot):
                        #rad_found.append(r * max(radii))
                        #break
                ## Conditions density and diff.
                #elif abs(r_dens - field_dens) <= dens_lim * delta_total and \
                    #r_diff < diff_lim:
                    #rad_found.append(r * max(radii))
                ## Conditions density and integral.
                #elif abs(r_dens - field_dens) <= dens_lim * delta_total and \
                    #(r - r_inte) >= int_lim * (1 - integ_tot):
                    #rad_found.append(r * max(radii))
                ### Conditions diff and integral.
                ##elif r_diff < diff_lim and \
                    ##(r - r_inte) >= int_lim * (1 - integ_tot):
                    ##rad_found.append(r * max(radii))


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

    #except:
        #pass

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
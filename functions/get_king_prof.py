# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 10:54:00 2013

@author: gabriel
"""

import numpy as np
from scipy.optimize import curve_fit
from scipy.ndimage.filters import gaussian_filter
from get_dens_prof import get_dens_prof as gdp


def new_hist(x_data, y_data, d_b):
    '''
    Generate new 2D histogram with increased bin width. Obtain also new center
    value.
    '''
    xmin, xmax = min(x_data), max(x_data)
    ymin, ymax = min(y_data), max(y_data)
    rang = [[xmin, xmax], [ymin, ymax]]
    binsxy = [int((xmax - xmin) / d_b), int((ymax - ymin) / d_b)]
    # hist is the 2D histogran, xedges & yedges store the edges of the bins
    hist, xedges, yedges = np.histogram2d(x_data, y_data, range=rang,
                                          bins=binsxy)
    # H_g is the 2D histogram with a gaussian filter applied
    h_g = gaussian_filter(hist, 2, mode='constant')
    # x_c_b, y_c_b store the x,y coordinates of the bin with the
    # maximum value in the 2D filtered histogram.
    x_c_b, y_c_b = np.unravel_index(h_g.argmax(), h_g.shape)

    return hist, x_c_b, y_c_b


def three_params(x, a, b, c):
    '''
    Three parameters King profile fit.
    '''
    d = backg_dens

    return c * (1 / np.sqrt(1 + (np.asarray(x) / a) ** 2) -
        1 / np.sqrt(1 + (b / a) ** 2)) ** 2 + d


def get_king_profile(clust_rad, backg_value, radii, ring_density, delta_xy,
    x_data, y_data, width_bin):
    '''
    Function to fit the 3-params King profile to a given radial density.
    The background density value is fixed and the core radius, tidal radius and
    maximum density are fitted.
    '''
    global backg_dens

    # Define values used by King profiles.
    max_dens = max(ring_density)
    backg_dens = backg_value
    r_t = clust_rad

    # Initial guesses for fit: r_c, r_t, max_dens.
    guess = (r_t / 2., r_t, max_dens)

    # Skip first radius value if it is smaller than the second value.
    if ring_density[0] > ring_density[1]:
        radii_k, ring_dens_k = radii, ring_density
    else:
        radii_k, ring_dens_k = radii[1:], ring_density[1:]

    flag_king_no_conver = True  # Flag that indicates no convergence.
    i, d_b_k = 1, width_bin
    # Iterate increasing the bin width until either King profile converges or
    # a value of twice the original bin width is reached.
    while flag_king_no_conver is True and i < width_bin:

        try:
            k_prof, k_pr_err = curve_fit(three_params, radii_k, ring_dens_k,
                guess)

            # If fit converged to tidal radius that extends beyond the maximum
            # range of the frame, discard it.
            if k_prof[1] > delta_xy:
                k_prof[1] = -1.
                k_pr_err[1][1] = 1.
                # Raise flag.
                flag_king_no_conver = True
            else:
                flag_king_no_conver = False

        except RuntimeError:
            flag_king_no_conver = True

        if flag_king_no_conver is True:

            # Increase bin width by 1 px.
            d_b_k = i + width_bin
            # Obtain 2D histo with new bin width.
            hist, x_c_b, y_c_b = new_hist(x_data, y_data, d_b_k)
            bin_center = [x_c_b, y_c_b]
            # Get new density profile.
            radii_k, ring_dens_k, poisson_error = gdp(hist, bin_center, d_b_k)

        i += 1

    # If 3-P King profile converged, calculate approximate number of cluster
    # members with eq (3) from Froebrich et al. (2007); 374, 399-408
    if flag_king_no_conver is False:
        print '3-P King profile obtained: %0.2f px.' % k_prof[1]
        x = 1 + (k_prof[1] / k_prof[0]) ** 2
        n_c_k = int(round((np.pi * max_dens * k_prof[0] ** 2) * (np.log(x) -
        4 + (4 * np.sqrt(x) + (x - 1)) / x)))
    else:
        print '  WARNING: King profile fitting did not converge.'
        # If 3-P King profile did not converge, pass dummy values
        n_c_k, k_prof, k_pr_err = -1, [-1, -1, -1], [[-1, -1, -1],
            [-1, -1, -1], [-1, -1, -1]]

    return k_prof, k_pr_err, d_b_k, n_c_k, flag_king_no_conver
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 10:54:00 2013

@author: gabriel
"""

import numpy as np
from scipy.optimize import curve_fit


def three_params(x, a, b, c):
    '''
    Three parameters King profile fit.
    '''
    d = backg_dens
    return c * (1 / np.sqrt(1 + (x / a) ** 2) -
        1 / np.sqrt(1 + (b / a) ** 2)) ** 2 + d


def get_king_profile(clust_rad, backg_value, radii, ring_density):
    '''
    Function to fit the 3-params King profile to a given radial density.
    The background density value is fixed and the core radius, tidal radius and
    maximum density are fitted.
    '''
    global backg_dens

    # Define values used by King profiles.
    max_dens = max(ring_density[0])
    backg_dens = backg_value
    r_t = clust_rad

    # Initial guesses for fit: r_c, r_t, max_dens.
    guess = (r_t / 2., r_t, max_dens)

    # Skip first radius value if it is smaller than the second value.
    flag_king_no_conver = False  # Flag that indicates no convergence.
    if ring_density[0][0] > ring_density[0][1]:
        try:
            k_prof, k_pr_err = curve_fit(three_params, radii[0],
                ring_density[0], guess)
        except RuntimeError:
            flag_king_no_conver = True
            k_prof, k_pr_err = [-1, -1, -1], [[-1, -1, -1], [-1, -1, -1],
                                [-1, -1, -1]]
    else:
        try:
            k_prof, k_pr_err = curve_fit(three_params, radii[0][1:],
                                         ring_density[0][1:], guess)
        except RuntimeError:
            flag_king_no_conver = True
            k_prof, k_pr_err = [-1, -1, -1], [[-1, -1, -1], [-1, -1, -1],
                                [-1, -1, -1]]

    # If 3-P King profile converged, calculate approximate number of cluster
    # members with eq (3) from Froebrich et al. (2007); 374, 399-408
    if flag_king_no_conver is False:
        x = 1 + (k_prof[1] / k_prof[0]) ** 2
        n_c_k = int(round((np.pi * max_dens * k_prof[0] ** 2) * (np.log(x) -
        4 + (4 * np.sqrt(x) + (x - 1)) / x)))
    else:
        n_c_k = -1

    # If fit converged to tidal radius too large.
    if k_prof[1] > 1000.:
        k_prof[1] = 999.99
        k_pr_err[1][1] = 998001.

    return k_prof, k_pr_err, n_c_k, flag_king_no_conver
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

    return c * (1 / np.sqrt(1 + (np.asarray(x) / a) ** 2) -
        1 / np.sqrt(1 + (b / a) ** 2)) ** 2 + d


def get_king_profile(clust_rad, backg_value, radii, ring_density):
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

    flag_king_no_conver = False  # Flag that indicates no convergence.
    try:
        k_prof, k_pr_err = curve_fit(three_params, radii_k, ring_dens_k, guess)
        print '3-P King profile obtained.'
    except RuntimeError:
        flag_king_no_conver = True
        print '  WARNING: King profile fitting did not converge.'

    # If 3-P King profile converged, calculate approximate number of cluster
    # members with eq (3) from Froebrich et al. (2007); 374, 399-408
    if flag_king_no_conver is False:

        # If fit converged to tidal radius too large, trim it.
        if k_prof[1] > 1000.:
            k_prof[1] = 999.99
            k_pr_err[1][1] = 998001.

        x = 1 + (k_prof[1] / k_prof[0]) ** 2
        n_c_k = int(round((np.pi * max_dens * k_prof[0] ** 2) * (np.log(x) -
        4 + (4 * np.sqrt(x) + (x - 1)) / x)))
    else:
        # If 3-P King profile did not converge, pass dummy values
        n_c_k, k_prof, k_pr_err = -1, [-1, -1, -1], [[-1, -1, -1], [-1, -1, -1],
                            [-1, -1, -1]]

    return k_prof, k_pr_err, n_c_k, flag_king_no_conver
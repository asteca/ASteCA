# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 10:54:00 2013

@author: gabriel
"""

import numpy as np
from scipy.optimize import curve_fit


def two_params(x, cd, rc, bg):
    '''
    Two parameters King profile fit.
    '''
    return bg + cd / (1 + (np.asarray(x) / rc) ** 2)


def three_params(x, rt, cd, rc, bg):
    '''
    Three parameters King profile fit.
    '''
    return cd * (1 / np.sqrt(1 + (np.asarray(x) / rc) ** 2) -
        1 / np.sqrt(1 + (rt / rc) ** 2)) ** 2 + bg


def get_king_profile(kp_flag, clust_rad, backg_value, radii, ring_density):
    '''
    Function to fit the 3-params King profile to a given radial density.
    The background density value is fixed and the core radius, tidal radius and
    maximum central density are fitted.
    '''

    # Flags that indicate either no convergence or that the fits were not
    # attempted.
    flag_2pk_conver, flag_3pk_conver = False, False

    # Check flag to run or skip.
    if kp_flag:

        # Background value is fixed.
        bg = backg_value
        # Initial guesses for fit: max_dens, rt, rc
        max_dens, rt, rc = max(ring_density), clust_rad, clust_rad / 2.
        guess2 = (max_dens, rc)
        guess3 = (max_dens, rc, rt)

        # Skip first radius value if it is smaller than the second value. This
        # makes it easier for the KP to converge.
        if ring_density[0] > ring_density[1]:
            radii_k, ring_dens_k = radii, ring_density
        else:
            radii_k, ring_dens_k = radii[1:], ring_density[1:]

        # Attempt to fit a 3-P King profile with the background value fixed.
        try:
            popt, pcov = curve_fit(lambda x, cd, rc,
                rt: three_params(x, rt, cd, rc, bg), radii_k, ring_dens_k,
                guess3)

            # Unpack tidal radius and its error.
            cd, rc, rt = popt
            e_rc = np.sqrt(pcov[1][1]) if pcov[1][1] > 0 else -1.
            e_rt = np.sqrt(pcov[2][2]) if pcov[2][2] > 0 else -1.
            flag_3pk_conver = True

            # If fit converged to tidal radius that extends beyond 100 times
            # the core radius, discard it.
            if rt > rc * 100. or rt <= 0. or rc <= 0.:
                # Raise flag to reject fit.
                flag_3pk_conver = False
        except RuntimeError:
            flag_3pk_conver = False

        # If 3-P King profile converged, ie: the tidal radius was found,
        # calculate approximate number of cluster members with eq (3) from
        # Froebrich et al. (2007); 374, 399-408
        if flag_3pk_conver is True:
            print 'Tidal radius obtained: %0.2f px.' % rt
            x = 1 + (rt / rc) ** 2
            n_c_k = int(round((np.pi * cd * rc ** 2) * (np.log(x) -
            4 + (4 * np.sqrt(x) + (x - 1)) / x)))
        else:
            print '  WARNING: tidal radius could not be obtained.'
            # If 3-P King profile did not converge, pass dummy values
            rt, e_rt, n_c_k = -1., -1., -1.

        if flag_3pk_conver is False:
            # Fit a 2P King profile first to obtain the maximum central
            # density and core radius.
            try:
                popt, pcov = curve_fit(lambda x, cd,
                    rc: two_params(x, cd, rc, bg), radii_k, ring_dens_k, guess2)
                # Unpaxk max density and core radius.
                cd, rc = popt
                # Obtain error in core radius.
                e_rc = np.sqrt(pcov[1][1]) if pcov[1][1] > 0 else -1.
                flag_2pk_conver = True
            except RuntimeError:
                flag_2pk_conver = False
                # Pass dummy values
                rc, e_rc, rt, e_rt, n_c_k, cd = -1., -1., -1., -1., -1., -1.
    else:
        # Pass dummy values
        rc, e_rc, rt, e_rt, n_c_k, cd = -1., -1., -1., -1., -1., -1.

    return rc, e_rc, rt, e_rt, n_c_k, cd, flag_2pk_conver, flag_3pk_conver
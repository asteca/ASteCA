# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 15:08:39 2014

@author: gabriel
"""

import numpy as np
from girardi_isochs_format import isoch_format as i_format
from get_isochs import get_isochs as gi
from get_met_ages_values import get_m_a_vls as gmav


def interp_isoch(isochrone):
    '''
    Interpolate extra color, magnitude and masses into the isochrone.
    '''
    N = 1500
    t, xp = np.linspace(0, 1, N), np.linspace(0, 1, len(isochrone[0]))
    # Store isochrone's interpolated values.
    isoch_inter = np.asarray([np.interp(t, xp, _) for _ in isochrone])

    return isoch_inter


def ip(ps_params, bf_flag):
    '''
    Read isochrones and parameters if best fit function is set to run.
    '''

    ip_list = []
    # Only read files of best fit method is set to run.
    if bf_flag:

        # Unpack.
        iso_path, cmd_select, iso_select, par_ranges = ps_params

        # Read Girardi metallicity files format.
        isoch_format = i_format(iso_select, cmd_select)

        # Obtain allowed metallicities and ages. Use the first photometric
        # system defined.
        # *WE ASUME ALL PHOTOMETRIC SYSTEMS CONTAIN THE SAME NUMBER OF
        # METALLICITY FILES*
        param_ranges, param_rs, met_f_filter, met_values, age_values = \
        gmav(iso_path, isoch_format, par_ranges)

        # Get isochrones and their parameter values.
        isoch_list = gi(cmd_select, met_f_filter, age_values, isoch_format)

        # Interpolate extra points into all isochrones.
        isochs_interp = [[] for _ in isoch_list]
        for i, _ in enumerate(isoch_list):
            for isoch in _:
                isochs_interp[i].append(interp_isoch(isoch))

        # Pack params.
        param_values = [met_values, age_values] + param_ranges[2:]
        ip_list = [isochs_interp, param_values, param_rs]

        iso_ver = {'10': '1.0', '11': '1.1', '12': '1.2S'}
        print ("PARSEC v{} theoretical isochrones read,".format(
            iso_ver[iso_select[-2:]]))
        lens = [len(_) for _ in param_values]
        total = reduce(lambda x, y: x * y, lens, 1)
        print ("interpolated and stored:\n"
        "  {} metallicity values (z),\n"
        "  {} age values (per z),\n"
        "  {} reddening values,\n"
        "  {} distance values,\n"
        "  {} mass values,\n"
        "  {} binary fraction values.".format(*lens))
        print "  = {:.1e} approx total models.".format(total)

    return ip_list
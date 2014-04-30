# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 16:38:10 2014

@author: gabriel
"""

import numpy as np
from isoch_likelihood import isoch_likelihood as i_l


def brute_force(err_lst, obs_clust, completeness, ip_list, sc_params,
                ga_params, cmd_sel):
    '''
    Brute force algorithm that computes the likelihoods for *all* the defined
    isochrones.
    '''

    isoch_list, isoch_ma, isoch_ed, ranges_steps = ip_list

    # Create lists with all possible
    e_lst = isoch_ed[0]
    d_lst = isoch_ed[1]

    # Initiate list that will hold the likelihood values telling us how well
    # each isochrone (syhtnetic cluster) fits the observed data.
    score = []

    tot_sols, i = len(isoch_list) * len(isoch_list[0]) * len(e_lst) *\
    len(d_lst), 0

    milestones = [25, 50, 75, 100]
    # Iterate through all metallicities.
    for m, m_isochs in enumerate(isoch_list):

        # Iterate through all ages.
        for a, isochs in enumerate(m_isochs):

            # Iterate through all extinction values.
            for e in e_lst:

                # Iterate through all distance modulus.
                for d in d_lst:

                    # Pass metallicity and age values for plotting purposes.
                    params = [isoch_ma[m][a][0], isoch_ma[m][a][1], e, d]

                    # Call likelihood function with m,a,e,d values.
                    likel_val = i_l(err_lst, obs_clust, completeness, sc_params,
                                    isoch_list[m][a], params, cmd_sel)

                    # Store the likelihood for each synthetic cluster.
                    score.append([likel_val, isoch_ma[m][a][0],
                                  isoch_ma[m][a][1], e, d])
                    i += 1

                    percentage_complete = (100.0 * (i + 1) / tot_sols)
                    while len(milestones) > 0 and percentage_complete >= \
                    milestones[0]:
                        print "  {}% done".format(milestones[0])
                        # Remove that milestone from the list.
                        milestones = milestones[1:]

    # Find index of function with smallest likelihood value.
    # This index thus points to the isochrone that best fits the observed
    # group of stars.
    best_fit_indx = np.argmin(zip(*score)[0])

    isoch_fit_params = [score[best_fit_indx][1:], score[best_fit_indx][0]]

    return isoch_fit_params
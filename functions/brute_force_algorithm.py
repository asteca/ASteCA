# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 16:38:10 2014

@author: gabriel
"""

import numpy as np
from get_likelihood import isoch_likelihood as i_l


def brute_force(err_lst, obs_clust, completeness, ip_list, st_d_bin_mr,
                ga_params, cmd_sel):
    '''
    Brute force algorithm that computes the likelihoods for *all* the defined
    isochrones.
    '''

    isoch_list, param_values = ip_list[:-1]

    # Unpack parameters values.
    m_lst, a_lst, e_lst, d_lst, mass_lst, bin_lst = param_values

    # Initiate list that will hold the likelihood values telling us how well
    # each isochrone (syhtnetic cluster) fits the observed data.
    score = [[], []]

    tot_sols, i = reduce(lambda x, y: x * y,
        [len(_) for _ in param_values], 1), 0

    milestones = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    # Iterate through all metallicities.
    for m_i, m in enumerate(m_lst):

        # Iterate through all ages.
        for a_i, a in enumerate(a_lst):

            # Iterate through all extinction values.
            for e in e_lst:

                # Iterate through all distance modulus.
                for d in d_lst:

                    # Iterate through all masses.
                    for mass in mass_lst:

                        # Iterate through all binary fractions.
                        for bin_frac in bin_lst:

                            params = [m, a, e, d, mass, bin_frac]

                            # Call likelihood function with m,a,e,d values.
                            isochrone = isoch_list[m_i][a_i]
                            # Call likelihood function with m,a,e,d values.
                            likel_val = i_l(err_lst, obs_clust, completeness,
                                st_d_bin_mr, isochrone, params, cmd_sel)
                            # Store the likelihood for each synthetic cluster.
                            score[0].append(likel_val)
                            score[1].append(params)

                            # Print percentage done.
                            i += 1
                            percentage_complete = (100.0 * (i + 1) / tot_sols)
                            while len(milestones) > 0 and \
                            percentage_complete >= milestones[0]:
                                print "  {}% done".format(milestones[0])
                                # Remove that milestone from the list.
                                milestones = milestones[1:]

    # Find index of function with smallest likelihood value.
    # This index thus points to the isochrone that best fits the observed
    # group of stars.
    best_fit_indx = np.argmin(score[0])

    isoch_fit_params = [score[1][best_fit_indx], score]

    return isoch_fit_params
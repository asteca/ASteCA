
import numpy as np
import likelihood


def main(lkl_method, e_max, bin_mr, err_lst, completeness, max_mag_syn,
         fundam_params, obs_clust, theor_tracks, R_V, ext_coefs, st_dist_mass,
         N_fc):
    '''
    Brute force algorithm that computes the likelihoods for *all* the defined
    isochrones.
    '''
    # Unpack parameters values.
    m_lst, a_lst, e_lst, d_lst, mass_lst, bin_lst = fundam_params

    # Initiate list that will hold the likelihood values telling us how well
    # each isochrone (synthetic cluster) fits the observed data.
    model_done = [[], []]

    # Count the number of total models/solutions to explore.
    tot_sols, i = reduce(lambda x, y: x * y,
                         [len(_) for _ in fundam_params], 1), 0
    milestones = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

    # Iterate through all metallicities.
    for m_i, m in enumerate(m_lst):

        # Iterate through all ages.
        for a_i, a in enumerate(a_lst):

            # Iterate through all extinction values.
            for e in e_lst:

                # Iterate through all distance moduli.
                for d in d_lst:

                    # Iterate through all masses.
                    for mass in mass_lst:

                        # Iterate through all binary fractions.
                        for bin_frac in bin_lst:
                            model = [m, a, e, d, mass, bin_frac]

                            # Call likelihood function with m,a,e,d values.
                            isochrone = theor_tracks[m_i][a_i]
                            # Call likelihood function with m,a,e,d values.
                            lkl = likelihood.main(
                                lkl_method, e_max, bin_mr, err_lst, obs_clust,
                                completeness, max_mag_syn, st_dist_mass,
                                isochrone, R_V, ext_coefs, N_fc, model)
                            # Store the likelihood for each synthetic
                            # cluster.
                            model_done[0].append(model)
                            model_done[1].append(lkl)

                            # Print percentage done.
                            i += 1
                            percentage_complete = (100.0 * (i + 1) /
                                                   tot_sols)
                            while len(milestones) > 0 and \
                                    percentage_complete >= milestones[0]:
                                best_fit_indx = np.argmin(model_done[1])
                                print (" {:>3}%  L={:.1f} ({:g}, {:g}, {:g},"
                                       " {:g}, {:g}, {:g})".format(
                                           milestones[0],
                                           model_done[1][best_fit_indx],
                                           *model_done[0][best_fit_indx]))
                                # Remove that milestone from the list.
                                milestones = milestones[1:]

    # Find index of function with smallest likelihood value.
    # This index thus points to the isochrone that best fits the observed
    # group of stars.
    best_fit_indx = np.argmin(model_done[1])

    isoch_fit_params = [model_done[0][best_fit_indx], model_done]

    return isoch_fit_params

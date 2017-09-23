
import numpy as np
import likelihood


def main(lkl_method, e_max, bin_mr, err_lst, completeness, max_mag_syn,
         fundam_params, obs_clust, theor_tracks, R_V, ext_coefs, st_dist_mass,
         N_fc, N_bf=1):
    """
    Brute force algorithm that computes the likelihoods for *all* the defined
    isochrones.

    # TODO what's stated below is addressed by issue #347

    It is possible that this algorithm returns a *larger* likelihood than the
    GA. This is counter-intuitive, but it is because the GA samples the IMF
    *each time* a mass value is checked, and the same mass value can be checked
    several times. The BF algorithm on the other hand does this only *once*
    per mass value.
    """

    m_lst, a_lst, e_lst, d_lst, mass_lst, bin_lst = fundam_params

    # Initiate list that will hold the likelihood values telling us how well
    # each isochrone (synthetic cluster) fits the observed data.
    model_done = [[], []]

    # Count the number of total models/solutions to explore.
    tot_sols, i = np.prod([len(_) for _ in fundam_params]), 0
    percs = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

    # Iterate through all metallicities.
    for m_i, met in enumerate(m_lst):

        # Iterate through all ages.
        for a_i, age in enumerate(a_lst):

            # Call likelihood function with m,a,e,d values.
            isochrone = theor_tracks[m_i][a_i]

            # Iterate through all extinction values.
            for e in e_lst:

                # Iterate through all distance moduli.
                for d in d_lst:

                    # Iterate through all masses.
                    for mass in mass_lst:

                        # Iterate through all binary fractions.
                        for bf in bin_lst:

                            # In place for #347
                            for _ in range(N_bf):
                                model = [met, age, e, d, mass, bf]
                                # Call likelihood function for this model.
                                lkl = likelihood.main(
                                    lkl_method, e_max, bin_mr, err_lst,
                                    obs_clust, completeness, max_mag_syn,
                                    st_dist_mass, isochrone, R_V, ext_coefs,
                                    N_fc, model)
                                # Store likelihood and synthetic cluster.
                                model_done[0].append(model)
                                model_done[1].append(lkl)

                            # Print percentage done.
                            i += 1
                            p_comp = (100.0 * (i + 1) / tot_sols)
                            while len(percs) > 0 and p_comp >= percs[0]:
                                best_fit_indx = np.argmin(model_done[1])
                                print (" {:>3}%  L={:.1f} ({:g}, {:g}, {:g},"
                                       " {:g}, {:g}, {:g})".format(
                                           percs[0],
                                           model_done[1][best_fit_indx],
                                           *model_done[0][best_fit_indx]))
                                # Remove that percentage value from the list.
                                percs = percs[1:]

    # Index of function with smallest likelihood value.
    # This index points to the model that best fits the observed cluster.
    best_fit_indx = np.argmin(model_done[1])

    isoch_fit_params = [model_done[0][best_fit_indx], model_done]

    return isoch_fit_params


import numpy as np
from ..synth_clust import synth_cluster

#############################################################
# # Timer function: http://stackoverflow.com/a/21860100/1391441
# from contextlib import contextmanager
# import time


# @contextmanager
# def timeblock(label):
#     start = time.clock()
#     try:
#         yield
#     finally:
#         end = time.clock()
#         print ('{} elapsed: {}'.format(label, end - start))
#############################################################


def tolstoy(synth_clust, obs_clust):
    '''
    Takes a synthetic cluster, compares it to the observed cluster and
    returns the weighted (log) likelihood value.
    This function follows the recipe given in Tolstoy & Saha (1996),
    Hernandez & Valls-Gabaud (2008) and Monteiro, Dias & Caetano (2010).
    '''

    # If synthetic cluster is empty, assign high likelihood value.
    if not synth_clust:
        tlst_lkl = 1.e09
    else:
        # synthetic cluster's photometry and errors.
        synth_phot, synth_errors = synth_clust[0]
        # Observed cluster's photometry and membership probabilities.
        obs_st, mem_probs = obs_clust

        # Square synthetic photometric errors.
        synth_errors = np.square(synth_errors)
        # Create list with the proper format for the synthetic cluster.
        syn_st = []
        for st_phot, st_e_phot in zip(zip(*synth_phot), zip(*synth_errors)):
            syn_st.append(zip(*[st_phot, st_e_phot]))

        # Probability for each observed star.
        cl_stars_probs = []
        # For each observed star.
        for o_st in obs_st:
            # o_st = [phot_1, phot_2, phot_3, ...]
            # o_phot_j = [phot_val, error]

            # Summatory over all synthetic stars.
            syn_sum = 0.
            # For each synthetic star.
            for s_st in syn_st:
                # s_st = [phot_1, phot_2, phot_3, ...]
                # s_phot_j = [phot_val, error]

                # Summatory for each photometric dimension.
                exp_sum, err_mult = 0., 1.
                # For each photometric dimension stored.
                for k, o_p in enumerate(o_st):
                    # Photometric difference squared, divided by the sum of
                    # the squared photometric errors, for this observed star,
                    # compared to this synthetic star, in this photometric
                    # dimension.
                    exp_sum += np.square(o_p[0] - s_st[k][0]) /\
                        (o_p[1] + s_st[k][1])
                    # Multiplication of square root of the sum of the squared
                    # photometric errors.
                    err_mult = err_mult * (o_p[1] + s_st[k][1])

                # Probability that the observed star 'o_st' is a star from
                # this synthetic model.
                syn_sum += np.exp(-0.5 * exp_sum) / np.sqrt(err_mult)

            # The final prob for this observed star is the sum over all
            # synthetic stars.
            cl_stars_probs.append(syn_sum)

        # Weight the probability of each observed star by its MP.
        # Use 1e-10 to avoid inf values in the np.log that follows.
        cl_stars_probs_w = cl_stars_probs * mem_probs
        cl_stars_probs_w[cl_stars_probs_w == 0.] = 1.e-10
        # Normalize by the number of synthetic stars, and sum the (negative)
        # logarithms to obtain the final likelihood.
        tlst_lkl = -sum(np.log(cl_stars_probs_w / len(syn_st)))

    return tlst_lkl


def dolphin(synth_clust, obs_clust):
    '''
    Takes a synthetic cluster, compares it to the observed cluster and
    returns the weighted (log) likelihood value.
    This function follows the recipe given in Dolphin (2002).
    '''

    # If synthetic cluster is empty, assign high likelihood value.
    if not synth_clust:
        # If synthetic cluster is empty, assign a large likelihood value.
        dolph_lkl = 1e09
    else:
        # with timeblock("  Histodd"):
        synth_phot = synth_clust[0][0]
        # Observed cluster's histogram and bin edges for each dimension.
        cl_histo, bin_edges = obs_clust

        # with timeblock("  log(lkl)"):
        # Histogram of the synthetic cluster, using the bin edges calculated
        # with the observed cluster.
        syn_histo = np.histogramdd(synth_phot, bins=bin_edges)[0]

        # Flatten N-dimensional histograms.
        cl_histo_f = np.array(cl_histo).flatten()
        syn_histo_f = np.array(syn_histo).flatten()
        # Assign small value to the zero elements in 'syn_histo_f'.
        syn_histo_f[syn_histo_f == 0] = 1.e-9

        # Obtain inverse logarithmic 'Poisson likelihood ratio'.
        dolph_lkl = 2. * (
            len(synth_phot[0]) - sum(cl_histo_f * np.log(syn_histo_f)))

    return dolph_lkl


def mighell(Q, P):
    '''
    '''

    if not Q.any():
        chi = 10000.
    else:

        # Observed cluster's histogram.
        cl_histo = P[0]
        # Bin edges for each dimension.
        b_rx, b_ry = P[1]

        # Magnitude and color for the synthetic cluster.
        syn_mags_cols = np.array(zip(*[Q[0], Q[2]]))
        # Histogram of the synthetic cluster, using the bin edges calculated
        # with the observed cluster.
        syn_histo = np.histogramdd(syn_mags_cols, bins=[b_rx, b_ry])[0]

        chi = 0.
        for el1 in zip(*(cl_histo, syn_histo)):
            for el2 in zip(*(el1[0], el1[1])):
                c = np.square(el2[0] + min(el2[0], 1) - el2[1]) / (el2[0] + 1)
                chi += c

        # print len(Q[0]), chi
        # import matplotlib.pyplot as plt
        # fig = plt.figure()
        # ax1 = fig.add_subplot(1, 3, 1)
        # ax2 = fig.add_subplot(1, 3, 2)
        # ax3 = fig.add_subplot(1, 3, 3)
        # ax1.imshow(d_1.transpose(), origin='lower', aspect='auto')
        # ax2.imshow(d_2.transpose(), origin='lower', aspect='auto')
        # ax3.scatter(P[4], P[2], c='r')
        # ax3.scatter(Q[0], Q[2], c='b')
        # text1 = 'chi = {:.2f}\n'.format(chi)
        # text2 = 'N = {}'.format(len(Q[0]))
        # text = text1 + text2
        # ax3.text(0.05, 0.9, text, transform=ax3.transAxes)
        # fig.subplots_adjust(hspace=1)
        # plt.show()

    return chi


def main(lkl_method, e_max, bin_mass_ratio, err_lst, obs_clust,
         completeness, st_dist_mass, isochrone, ext_coefs, N_fc,
         synth_cl_params):
    '''
    Generate synthetic cluster with an isochrone of given values for
    metallicity and age. Match the synthetic cluster to the observed cluster.
    '''

    # Generate synthetic cluster..
    # with timeblock("synth_cl"):
    synth_clust = synth_cluster.main(
        e_max, bin_mass_ratio, err_lst, completeness, st_dist_mass, isochrone,
        ext_coefs, N_fc, synth_cl_params)
    # synth_clust = [photom_dims, photom_errs, binar_idx, extra_info]
    # photom_dims = [mag1, mag2, ..., magN, col1, col2, ..., colM]
    # photom_errs = [e_m1, e_m2, ..., e_mN, e_c1, e_c2, ..., e_CM]
    # binar_idx = [ids of binary systems for stars in photom_dims]
    # extra info = lists with additional information, from the theoretical
    #              isochrones.

    # Obtain the likelihood matching the synthetic and observed clusters.
    if lkl_method == 'tolstoy':
        likelihood = tolstoy(synth_clust, obs_clust)
    else:
        likelihood = dolphin(synth_clust, obs_clust)
        # likelihood = mighell(synth_clust, obs_clust)

    return likelihood

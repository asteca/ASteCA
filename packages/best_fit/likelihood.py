
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
    synth_phot, synth_errors = synth_clust[0]

    # If synthetic cluster is empty, assign high likelihood value.
    if not synth_phot[0].any():
        likelihood = 1.e09
    else:

        # Unpack observed cluster's photometry, its squared errors, and its
        # membership probabilities.
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
            # synthetic stars. Use 1e-10 to avoid nan and inf values in the
            # calculations that follow.
            cl_stars_probs.append(max(syn_sum, 1.e-10))

        # Final score: sum log likelihoods for each star in cluster, weighted
        # by its probability.
        likelihood = -sum(np.log(cl_stars_probs * mem_probs / len(syn_st)))

    return likelihood


def dolphin_plot(Q, P, b_rx, b_ry, cl_histo, syn_histo, likel):
    '''
    Plots Dolphin histograms.
    '''

    print len(Q[0]), likel
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MultipleLocator

    # Define grid for 2D histograms.
    X, Y = np.meshgrid(b_rx, b_ry)
    # Extent of plots.
    x_extend = [min(P[-1][0]) - 1., max(P[-1][0]) + 1.]
    y_extend = [max(P[-1][1]) + 1., min(P[-1][1]) - 6.]

    # Define subplots.
    fig = plt.figure()
    ax0 = fig.add_subplot(1, 4, 1)
    ax0.set_xlim(x_extend)
    ax0.set_ylim(y_extend)
    ax1 = fig.add_subplot(1, 4, 2)
    ax2 = fig.add_subplot(1, 4, 3)
    ax2.minorticks_on()
    ax2.yaxis.set_major_locator(MultipleLocator(1.0))
    ax3 = fig.add_subplot(1, 4, 4)
    ax3.minorticks_on()
    ax3.yaxis.set_major_locator(MultipleLocator(1.0))

    # Scatter plot for cluster region + synth cluster.
    ax0.scatter(P[-1][0], P[-1][1], c='r', label='Obs clust')
    ax0.scatter(Q[0], Q[2], c='b', label='Synth clust')
    for x_ed in b_rx:
        # vertical lines
        ax0.axvline(x_ed, linestyle=':', color='k', zorder=1)
    for y_ed in b_ry:
        # horizontal lines
        ax0.axhline(y_ed, linestyle=':', color='k', zorder=1)
    ax0.legend(fancybox=True, loc='lower left', scatterpoints=1)

    # Scatter plot for synth clust.
    ax1.scatter(Q[0], Q[2], c='b', label='Synth clust')
    ax1.invert_yaxis()
    ax1.legend(fancybox=True, loc='lower left', scatterpoints=1)

    # Observed cluster 2D histo.
    # Rotate and flip H.
    H = np.rot90(cl_histo)
    H = np.flipud(H)
    ax2.pcolormesh(X, Y, H, cmap=plt.cm.Reds)
    ax2.text(0.05, 0.95, 'Obs clust 2D histo', transform=ax2.transAxes,
             fontsize=15)

    # Synthetic cluster 2D histo.
    H = np.rot90(syn_histo)
    H = np.flipud(H)
    ax3.pcolormesh(X, Y, H, cmap=plt.cm.Reds)
    text = r'Synth clust 2D histo ($-\ln(PLR) \approx {:.2f}$)'.format(
        likel)
    ax3.text(0.05, 0.95, text, transform=ax3.transAxes, fontsize=15)

    # Invert histograms axis.
    ax2.invert_yaxis()
    ax3.invert_yaxis()

    # Set limits for 2D histos.
    ax2.set_xlim(x_extend)
    ax2.set_ylim(y_extend)
    ax3.set_xlim(x_extend)
    ax3.set_ylim(y_extend)

    fig.subplots_adjust(hspace=1)
    plt.show()


def dolphin(Q, obs_clust):
    '''
    Takes a synthetic cluster, compares it to the observed cluster and
    returns the weighted (log) likelihood value.
    This function follows the recipe given in Dolphin (2002).
    '''

    if not Q.any():
        # If synthetic cluster is empty, assign a large likelihood value.
        poiss_lkl = 1e09
    else:
        # with timeblock("  Histodd"):
        # Observed cluster's histogram.
        cl_histo = obs_clust[0]
        # Bin edges for each dimension.
        b_rx, b_ry = obs_clust[1]

        # Magnitude sand colors for the synthetic cluster.
        syn_mags_cols = np.array(zip(*[Q[0], Q[2]]))
        # Histogram of the synthetic cluster, using the bin edges calculated
        # with the observed cluster.
        syn_histo = np.histogramdd(syn_mags_cols, bins=[b_rx, b_ry])[0]

        # with timeblock("  log(lkl)"):
        # Small value used to replace zeros.
        epsilon = 1e-10
        # Obtain inverse logarithmic 'Poisson likelihood ratio'.
        poiss_lkl = len(Q[0])
        for el1 in zip(*(cl_histo, syn_histo)):
            for el2 in zip(*(el1[0], el1[1])):
                c = -1. * el2[0] * np.log(max(el2[1], epsilon))
                poiss_lkl += c

        # Call this function to see histograms produced.
        # *IMPORTANT*: The list passed in obs_clust_prepare must be modified
        # for this block to work.
        # dolphin_plot(Q, P, b_rx, b_ry, cl_histo, syn_histo, poiss_lkl)

    return poiss_lkl


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
    Call with an isochrone of given values for metallicity and age and supply
    the extinction and distance modulus values to move that isochrone. Use
    that isochrone to generate a synthetic cluster with those parameters and
    finally compare it with the observed cluster.

    synth_clust = [photom_dims, photom_errs, binar_idx, extra_info]
    photom_dims = [mag1, mag2, ..., magN, col1, col2, ..., colM]
    photom_errs = [e_m1, e_m2, ..., e_mN, e_c1, e_c2, ..., e_CM]
    binar_idx = [ids of binary systems for stars in photom_dims]
    extra info = lists with additional information, from the theoretical
                 isochrones.
    '''

    # Generate synthetic cluster using this "moved" isochrone and a mass
    # distribution.
    # with timeblock("synth_cl"):
    synth_clust = synth_cluster.main(
        e_max, bin_mass_ratio, err_lst, completeness, st_dist_mass, isochrone,
        ext_coefs, N_fc, synth_cl_params)

    # Call function to obtain the likelihood by comparing the synthetic cluster
    # with the observed cluster.
    if lkl_method == 'tolstoy':
        likelihood = tolstoy(synth_clust, obs_clust)
    else:
        likelihood = dolphin(synth_clust, obs_clust)
        # likelihood = mighell(synth_clust, obs_clust)

    return likelihood

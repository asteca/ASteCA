
import numpy as np
from ..synth_clust import synth_cluster

# ############################################################
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
# ############################################################


def tolstoy(synth_clust, obs_clust):
    '''
    Compares the observed cluster with a synthetic cluster,
    returns the weighted (log) likelihood value.
    This function follows the recipe given in Tolstoy & Saha (1996),
    Hernandez & Valls-Gabaud (2008) and Monteiro, Dias & Caetano (2010).

    Likelihood form:

    \begin{equation}
    -\ln L = N\ln M -
             \sum\limits_{i=1}^{N} \ln
             \left\{
                P_i \sum\limits_{j=1}^M \frac{\exp
                    \left[
                        -\frac{1}{2}
                            \left(
                                \sum_{k=1}^D \frac{(f_{ik}-g_{jk})^2}
                                {\sigma_{ik}^2 + \rho_{jk}^2}
                            \right)
                    \right]}
                {\prod_{k=1}^D \sqrt{\sigma_{ik}^2 + \rho_{jk}^2}}
            \right\}
    \end{equation}
    '''

    # If synthetic cluster is empty, assign large likelihood value.
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

        # Photometric difference (observed - synthetic), for all dimensions.
        phot_dif = np.array(obs_st)[:, None, :, 0] -\
            np.array(syn_st)[None, :, :, 0]
        # Sum of squared photometric errors, for all dimensions.
        sigma_sum = np.array(obs_st)[:, None, :, 1] +\
            np.array(syn_st)[None, :, :, 1]

        # Sum for all photometric dimensions.
        Dsum = (np.square(phot_dif) / sigma_sum).sum(axis=-1)
        # Product of summed squared sigmas.
        sigma_prod = np.prod(sigma_sum, axis=-1)
        # All elements inside synthetic stars summatory.
        sum_M_j = np.exp(-0.5 * Dsum) / np.sqrt(sigma_prod)
        # Sum for all synthetic stars.
        sum_M = np.sum(sum_M_j, axis=-1)

        # Multiply by membership probabilities.
        sum_M_MP = mem_probs * sum_M
        # Remove 0. elements before applying the logarithm below.
        sum_N = sum_M_MP[sum_M_MP != 0.]

        # Final negative logarithmic likelihood
        tlst_lkl = len(obs_st) * np.log(len(syn_st)) - sum(np.log(sum_N))

    return tlst_lkl


def dolphin(synth_clust, obs_clust):
    '''
    Takes a synthetic cluster, compares it to the observed cluster and
    returns the weighted (log) likelihood value.
    This function follows the recipe given in Dolphin (2002).
    '''

    # If synthetic cluster is empty, assign large likelihood value.
    if not synth_clust:
        dolph_lkl = 1.e09
    else:
        synth_phot = synth_clust[0][0]
        # Observed cluster's histogram and bin edges for each dimension.
        bin_edges = obs_clust[0]
        # Indexes of n_i=0 elements in flattened observed cluster array,
        # and the array with no n_i=0 elements.
        cl_z_idx, cl_histo_f_z, dolphin_cst = obs_clust[-3:]

        # Histogram of the synthetic cluster, using the bin edges calculated
        # with the observed cluster.
        syn_histo = np.histogramdd(synth_phot, bins=bin_edges)[0]
        # Flatten N-dimensional histogram.
        syn_histo_f = np.array(syn_histo).ravel()
        # Remove all bins where n_i = 0 (no observed stars).
        syn_histo_f_z = syn_histo_f[cl_z_idx]

        # Assign small value to the m_i = 0 elements in 'syn_histo_f_z'.
        # The value equals 1 star divided among all empty bins.
        syn_histo_f_z[syn_histo_f_z == 0] =\
            1. / max(np.count_nonzero(syn_histo_f_z == 0), 1.)

        # Obtain inverse logarithmic 'Poisson likelihood ratio'.
        # 2 * (M - sum(n_i * ln(m_i)) + sum(n_i * ln(n_i)) - N)
        dolph_lkl = 2. * (
            synth_phot[0].size - np.sum(cl_histo_f_z * np.log(syn_histo_f_z)) +
            dolphin_cst)

    return dolph_lkl


def mighell(synth_clust, obs_clust):
    '''
    '''
    if not synth_clust:
        mig_chi = 1.e09
    else:
        # Observed cluster's histogram and bin edges for each dimension.
        bin_edges = obs_clust[0]
        # Observed cluster's flattened histogram and indexes of n_i=0 elements.
        cl_histo_f, cl_z_idx = obs_clust[2:4]

        # Synthetic cluster.
        synth_phot = synth_clust[0][0]
        # Histogram of the synthetic cluster, using the bin edges calculated
        # with the observed cluster.
        syn_histo = np.histogramdd(synth_phot, bins=bin_edges)[0]

        # Flatten N-dimensional histogram.
        syn_histo_f = np.array(syn_histo).ravel()
        # Indexes of bins that are empty in both arrays.
        z = cl_z_idx[0] | (syn_histo_f != 0)
        # Remove those bins.
        cl_histo_f_z, syn_histo_f_z = cl_histo_f[z], syn_histo_f[z]

        # Final chi.
        mig_chi = np.sum(np.square(
            cl_histo_f_z + np.clip(cl_histo_f_z, 0, 1) - syn_histo_f_z) /
            (cl_histo_f_z + 1.))

    return mig_chi


def main(lkl_method, e_max, bin_mass_ratio, err_lst, obs_clust, completeness,
         max_mag_syn, st_dist_mass, isochrone, R_V, ext_coefs, N_fc,
         synth_cl_params):
    '''
    Generate synthetic cluster with an isochrone of given values for
    metallicity and age. Match the synthetic cluster to the observed cluster.
    '''

    # Generate synthetic cluster.
    # with timeblock("synth_cl"):
    synth_clust = synth_cluster.main(
        e_max, bin_mass_ratio, err_lst, completeness, max_mag_syn,
        st_dist_mass, isochrone, R_V, ext_coefs, N_fc, synth_cl_params)
    # synth_clust = [photom_dims, photom_errs, binar_idx, extra_info]
    # photom_dims = [mag1, mag2, ..., magN, col1, col2, ..., colM]
    # photom_errs = [e_m1, e_m2, ..., e_mN, e_c1, e_c2, ..., e_CM]
    # binar_idx = [ids of binary systems for stars in photom_dims]
    # extra info = lists with additional information, from the theoretical
    #              isochrones.

    # Obtain the likelihood matching the synthetic and observed clusters.
    if lkl_method == 'tolstoy':
        likelihood = tolstoy(synth_clust, obs_clust)
    elif lkl_method == 'dolphin':
        likelihood = dolphin(synth_clust, obs_clust)
    elif lkl_method == 'mighell':
        likelihood = mighell(synth_clust, obs_clust)

    return likelihood


import numpy as np
from ..synth_clust import tracksPrep, extin_coefs, imf
from ..best_fit import completenessPercEstim


def main(pd, clp, td):
    """
    Prepare tracks and other required data
    """
    if pd['best_fit_algor'] == 'n':
        return td

    # Keys added: 'ext_coefs', 'N_fc', 'm_ini_idx',
    # 'st_dist_mass', 'theor_tracks', 'mean_bin_mr', 'err_norm_rand',
    # 'binar_probs', 'ext_unif_rand'
    td = tracksPrep.main(td, **pd)

    # Obtain extinction coefficients.
    td['ext_coefs'] = extin_coefs.main(
        pd['cmd_systs'], pd['filters'], pd['colors'])

    # Obtain mass distribution using the selected IMF.
    Nmets = len(td['fundam_params'][0])
    td['st_dist_mass'] = imf.main(pd['IMF_name'], Nmets, pd['Max_mass'])

    td['err_norm_rand'], td['binar_probs'], td['ext_unif_rand'] = randVals(
        td['st_dist_mass'], td['theor_tracks'])

    td['ed_compl_vals'] = []
    if clp['completeness'][-1] is True:
        td['ed_compl_vals'] = completenessPercEstim.main(clp, td)

    return td


def randVals(st_dist_mass, theor_tracks):
    """
    Generate lists of random values used by the synthetic cluster generating
    function.
    """
    # This is the maximum number of stars that will ever be interpolated into
    # an isochrone
    N_mass = 0
    for sdm in st_dist_mass:
        N_mass = max(len(sdm[0]), N_mass)

    # Number of metallicity values defined
    N_mets = theor_tracks.shape[0]

    binar_probs = []
    for _ in range(N_mets):
        # Uniform probabilities, used for the binary assignment (one per
        # metallicity)
        binar_probs.append(np.random.uniform(0., 1., len(st_dist_mass[_][0])))

    err_norm_rand = np.random.normal(0., 1., (N_mets, N_mass))

    # For the move_isoch() function. In place for #174
    # ext_unif_rand.append(np.random.uniform(0., 1., N_mass))

    ext_unif_rand = np.random.uniform(0., 1., (N_mets, N_mass))

    return err_norm_rand, binar_probs, ext_unif_rand

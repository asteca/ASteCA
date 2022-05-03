
import numpy as np
from ..synth_clust import tracksPrep, extin_coefs, imf


def main(pd, clp, td):
    """
    Prepare tracks and other required data
    """
    if pd['best_fit_algor'] == 'n':
        return td

    # Interpolate the isochrones and generate binary systems
    td = tracksPrep.main(td, **pd)

    # Obtain extinction coefficients.
    td['ext_coefs'] = extin_coefs.main(
        pd['cmd_systs'], pd['filters'], pd['colors'])

    # Obtain mass distribution using the selected IMF.
    Nmets = len(td['fundam_params'][0])
    td['st_dist_mass'] = imf.main(pd['IMF_name'], Nmets, pd['Max_mass'])

    td['rand_norm_vals'], td['rand_unif_vals'] = randVals(
        td['theor_tracks'], td['st_dist_mass'])

    return td


def randVals(theor_tracks, st_dist_mass):
    """
    Generate lists of random values used by the synthetic cluster generating
    function.
    """
    # This is the maximum number of stars that will ever be interpolated into
    # an isochrone
    N_isoch, N_mass = theor_tracks.shape[-1], 0
    for sdm in st_dist_mass:
        N_mass = max(len(sdm[0]), N_mass, N_isoch)

    # Used by `move_isochrone()` and `add_errors`
    rand_norm_vals = np.random.normal(0., 1., (2, N_mass))

    # Used by `move_isochrone()`, `binarity()`, `completeness_rm()`
    rand_unif_vals = np.random.uniform(0., 1., (3, N_mass))

    return rand_norm_vals, rand_unif_vals

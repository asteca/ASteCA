
import numpy as np
from . import zaWAverage
from . import move_isochrone
from . import cut_max_mag
from . import mass_distribution
from . import mass_interp
from . import binarity
from . import completeness_rm
from . import add_errors


def main(
    fundam_params, varIdxs, model, theor_tracks, completeness, max_mag_syn,
    st_dist_mass, R_V, ext_coefs, N_fc, err_pars, m_ini_idx, binar_flag,
        extra_pars_flag=False):
    """
    Takes an isochrone and returns a synthetic cluster created according to
    a certain mass distribution.

    The synthetic cluster returned has the shape:

    synth_clust = [photometry, binary_idxs + extra_pars]

    photometry = [photom, errors]
    photom = [f1, f2, ..., fF, c1, c2, ..., cC]
    (where F and C are the total number of filters and colors defined)

    errors = [ef1, ef2, ..., efF, ec1, ec2, ..., ecC]
    (photometric errrors for each photometric dimension defined)

    Correct indexes of binary systems after completeness removal.
    binary_idxs = [i1, i2, ..., iN]

    Lists containing the theoretical tracks extra parameters.
    extra_pars = [l1, l2, ..., l6]
    """

    # If 'extra_pars_flag' is True and synth_clust = np.array([])
    sigma, extra_pars, isoch_moved, mass_dist, isoch_binar, isoch_compl =\
        [[] for _ in range(6)]

    # Return proper values for fixed parameters and parameters required
    # for the (z, log(age)) isochrone averaging.
    model_proper, z_model, a_model, ml, mh, al, ah = properModel(
        fundam_params, model, varIdxs)

    # Generate a weighted average isochrone from the (z, log(age)) values in
    # the 'model'.
    isochrone = zaWAverage.main(
        theor_tracks, m_ini_idx, fundam_params, z_model, a_model, ml, mh, al,
        ah)

    # Extract parameters
    e, d, M_total, bin_frac = model_proper

    # Move theoretical isochrone using the values 'e' and 'd'.
    isoch_moved = move_isochrone.main(
        isochrone, e, d, R_V, ext_coefs, N_fc, binar_flag, m_ini_idx)

    # Get isochrone minus those stars beyond the magnitude cut.
    isoch_cut = cut_max_mag.main(isoch_moved, max_mag_syn)

    # # In place for #358
    # return isoch_cut.T[:, :3]

    # Empty list to pass if at some point no stars are left.
    synth_clust = np.array([])
    # Check for an empty array.
    if isoch_cut.any():
        # Mass distribution to produce a synthetic cluster based on
        # a given IMF and total mass.
        mass_dist = mass_distribution.main(st_dist_mass, M_total)

        # Interpolate masses in mass_dist into the isochrone rejecting those
        # masses that fall outside of the isochrone's mass range.
        # This destroys the order by magnitude.
        isoch_mass = mass_interp.main(isoch_cut, mass_dist, m_ini_idx)

        if isoch_mass.any():
            # Assignment of binarity.
            isoch_binar = binarity.main(isoch_mass, bin_frac, m_ini_idx, N_fc)

            # Completeness limit removal of stars.
            isoch_compl = completeness_rm.main(isoch_binar, completeness)

            if isoch_compl.any():
                # Get errors according to errors distribution.
                synth_clust, sigma, extra_pars = add_errors.main(
                    isoch_compl, err_pars, m_ini_idx, binar_flag,
                    extra_pars_flag)

    if extra_pars_flag is False:
        # Only pass the photometry, used by the likelihood function
        return synth_clust
    return synth_clust, sigma, extra_pars, isoch_moved, mass_dist,\
        isoch_binar, isoch_compl


def properModel(fundam_params, model, varIdxs):
    """
    Define the 'proper' model with values for (z, a) taken from its grid,
    and filled values for those parameters that are fixed.

    Parameters
    ----------
    model : array
      Array of *free* fundamental parameters only (ie: in varIdxs).

    Returns
    -------
    model_proper : list
      Stores (E_BV, dm, Mass, b_fr) including the fixed parameters that are
      missing from 'model'.
    z_model, a_model : floats
      The (z, a) values for this model's isochrone.
    ml, mh, al, ah : ints
      Indexes of the (z, a) values in the grid that define the box that enclose
      the (z_model, a_model) values.

    """

    model_proper, j = [], 0
    for i, par in enumerate(fundam_params):
        # Check if this parameter is one of the 'free' parameters.
        if i in varIdxs:
            # If it is the parameter metallicity.
            if i == 0:
                # Select the closest value in the array of allowed values.
                mh = min(len(par) - 1, np.searchsorted(par, model[i - j]))
                ml = mh - 1
                # Define the model's z value
                z_model = model[i - j]
            # If it is the parameter log(age).
            elif i == 1:
                # Select the closest value in the array of allowed values.
                ah = min(len(par) - 1, np.searchsorted(par, model[i - j]))
                al = ah - 1
                a_model = model[i - j]
            else:
                model_proper.append(model[i - j])
        else:
            if i == 0:
                ml = mh = 0
                z_model = fundam_params[0][0]
            elif i == 1:
                al = ah = 0
                a_model = fundam_params[1][0]
            else:
                model_proper.append(par[0])
            j += 1

    return model_proper, z_model, a_model, ml, mh, al, ah

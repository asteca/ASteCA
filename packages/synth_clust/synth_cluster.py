
import numpy as np
from . import zaWAverage
from . import move_isochrone
from . import cut_max_mag
from . import mass_distribution
from . import mass_interp
from . import binarity
from . import completeness_rm
from . import add_errors

# import matplotlib.pyplot as plt


def main(
    fundam_params, varIdxs, model, completeness, err_lst,
    max_mag_syn, ext_coefs, binar_flag, mean_bin_mr, N_fc, m_ini_idx,
    st_dist_mass, theor_tracks, err_norm_rand, binar_probs, ext_unif_rand,
        transpose_flag=True):
    """
    Takes an isochrone and returns a synthetic cluster created according to
    a certain mass distribution.

    The synthetic cluster returned has the shape:

    * transpose_flag = True

    synth_clust = [mag, c1, (c2)]

    where c1 and c2 colors defined.

    * transpose_flag = False

    synth_clust = [mag, c1, (c2), mag_b, c1_b, (c2_b), m_ini_1, m_ini_1]

    where 'm_ini_1, m_ini_2' are the masses of the single/binary systems. The
    single systems only have 'm_ini_1' defined, and a '-99' stored in 'm_ini_2'

    """

    # Return proper values for fixed parameters and parameters required
    # for the (z, log(age)) isochrone averaging.
    model_proper, z_model, a_model, ml, mh, al, ah = properModel(
        fundam_params, model, varIdxs)

    # Generate a weighted average isochrone from the (z, log(age)) values in
    # the 'model'.
    isochrone = zaWAverage.main(
        theor_tracks, fundam_params, binar_flag, m_ini_idx, z_model, a_model,
        ml, mh, al, ah)

    # Extract parameters
    e, d, M_total, bin_frac, R_V = model_proper

    # Move theoretical isochrone using the values 'e' and 'd'.
    isoch_moved = move_isochrone.main(
        isochrone, e, d, R_V, ext_coefs, N_fc, ext_unif_rand[ml],
        m_ini_idx, binar_flag)

    # Get isochrone minus those stars beyond the magnitude cut.
    isoch_cut = cut_max_mag.main(isoch_moved, max_mag_syn)

    # # In place for #358
    # return isoch_cut.T[:, :3]

    # Empty list to pass if at some point no stars are left.
    synth_clust = np.array([])
    if isoch_cut.any():

        # Return the isochrone with the proper total mass.
        mass_dist = mass_distribution.main(
            st_dist_mass[ml], mean_bin_mr, bin_frac, M_total)

        # Interpolate masses in mass_dist into the isochrone rejecting those
        # masses that fall outside of the isochrone's mass range.
        # This destroys the order by magnitude.
        isoch_mass = mass_interp.main(isoch_cut, m_ini_idx, mass_dist)

        # # In place for testing #508
        # plt.subplot(121)
        # # plt.scatter(isoch_mass[4], isoch_mass[3], c='r')
        # plt.scatter(isoch_mass[1], isoch_mass[0], c='g')
        # # plt.scatter(isoch_cut[4], isoch_cut[3], c='b', marker='x')
        # # plt.scatter(isoch_cut[1], isoch_cut[0], c='cyan', marker='x')
        # plt.gca().invert_yaxis()
        # # plt.hist(isoch_mass[0], 25)

        if isoch_mass.any():
            # Assignment of binarity.
            isoch_binar = binarity.main(
                isoch_mass, bin_frac, m_ini_idx, N_fc, binar_probs[ml])

            # # In place for testing #508
            # plt.subplot(122)
            # plt.scatter(isoch_binar[1], isoch_binar[0], c='g')
            # plt.gca().invert_yaxis()
            # plt.show()

            # Completeness limit removal of stars.
            isoch_compl = completeness_rm.main(isoch_binar, completeness)

            # # In place for testing #508
            # plt.subplot(122)
            # plt.scatter(isoch_compl[1], isoch_compl[0], c='g')
            # plt.gca().invert_yaxis()
            # # plt.hist(isoch_compl[0], 25)
            # plt.show()

            if isoch_compl.any():
                # Get errors according to errors distribution.
                synth_clust = add_errors.main(
                    isoch_compl, err_lst, err_norm_rand[ml])

                # Transposing is necessary for np.histogramdd() in the
                # likelihood
                if transpose_flag:
                    synth_clust = synth_clust[:sum(N_fc)].T

    return synth_clust


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
      Stores (E_BV, dm, Mass, b_fr, R_V) including the fixed parameters that
      are missing from 'model'.
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


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
    model, transpose_flag, DR_dist, alpha, model_proper, varIdxs, completeness,
    err_lst, max_mag_syn, N_obs_stars, fundam_params, ed_compl_vals, ext_coefs,
    N_fc, m_ini_idx, st_dist_mass, theor_tracks, rand_norm_vals,
        rand_unif_vals):
    """
    Takes an isochrone and returns a synthetic cluster created according to
    a certain mass distribution.

    The synthetic cluster returned has the shape:

    * transpose_flag = True

    synth_clust = [mag, c1, (c2)]

    where c1 and c2 colors defined.

    * transpose_flag = False

    synth_clust = [mag, c1, (c2), mag_b, c1_b, (c2_b), m_ini_1, m_ini_2]

    where 'm_ini_1, m_ini_2' are the primary and secondary masses of the
    binary systems. The single systems only have a '0' stored in 'm_ini_2'.
    """

    # Return proper values for fixed parameters and parameters required
    # for the (z, log(age)) isochrone averaging.
    model_proper, ml, mh, al, ah = properModel(
        fundam_params, model_proper, model, varIdxs)

    # Extract parameters
    met, age, e, dr, d, beta, R_V = model_proper

    # Generate a weighted average isochrone from the (z, log(age)) values in
    # the 'model'.
    isochrone = zaWAverage.main(
        theor_tracks, fundam_params, m_ini_idx, met, age, ml, mh, al, ah)

    # Move theoretical isochrone using the values 'e' and 'd'.
    isoch_moved = move_isochrone.main(
        isochrone, e, dr, d, R_V, ext_coefs, N_fc, DR_dist,
        rand_norm_vals[0], rand_unif_vals[0], m_ini_idx)

    # Get isochrone minus those stars beyond the magnitude cut.
    isoch_cut = cut_max_mag.main(isoch_moved, max_mag_syn)
    # In place for MiMO testing
    # return isoch_cut
    # # In place for #358
    # return isoch_cut.T[:, :3]

    # Empty list to pass if at some point no stars are left.
    synth_clust, M_total = np.array([]), 0.

    if not isoch_cut.any():
        return synth_clust, M_total

    mass_ini = isoch_cut[m_ini_idx]

    # Return the mass distribution within the 'isoch_cut' range
    mass_dist, msk_m, M_total_arr = mass_distribution.main(
        mass_ini, st_dist_mass[ml])

    if msk_m.sum() == 0:
        return synth_clust, M_total

    # # The estimation of the total fraction of stars removed by the
    # # completeness function is replaced by an exact removal *after*
    # this function is applied. This does not affect the performance
    # significantly. REMOVE
    #
    # compl_rm_perc = 0
    # if completeness[-1] is True:
    #     # Indexes of elements in ed_compl_vals closest the (e, d) values
    #     idx_a = np.argmin(abs(a_model - ed_compl_vals[0]))
    #     idx_e = np.argmin(abs(e - ed_compl_vals[1]))
    #     idx_d = np.argmin(abs(d - ed_compl_vals[2]))
    #     compl_rm_perc = ed_compl_vals[3][idx_a, idx_e, idx_d]
    # # Total number of stars, corrected by the removal process by the
    # # completeness function below
    # N_stars = int(N_obs_stars / (1 - compl_rm_perc))
    # mass_dist = mass[msk_m][:N_stars]
    # # Total mass estimation for this 'N_stars' value
    # M_total = st_dist_mass[ml][1][msk_m][:N_stars][-1]

    # Interpolate masses in mass_dist into the isochrone rejecting those
    # masses that fall outside of the isochrone's mass range.
    # This destroys the order by magnitude.
    isoch_mass = mass_interp.main(isoch_cut, mass_ini, mass_dist)

    if not isoch_mass.any():
        return synth_clust, M_total

    # Assignment of binarity.
    isoch_binar = binarity.main(
        isoch_mass, alpha, beta, m_ini_idx, N_fc, rand_unif_vals[1])

    # Use a different array of random uniform values to de-couple the
    # process from the one applied in binarity.main()
    # Completeness removal of stars.
    isoch_compl, msk_cr = completeness_rm.main(
        isoch_binar, completeness, rand_unif_vals[2])

    # Keep the same number of synthetic stars as observed stars
    isoch_compl = isoch_compl[:, :N_obs_stars]

    if not isoch_compl.any():
        return synth_clust, M_total

    # Apply 'msk_cr' to correct the mass for the completeness (if applied)
    if msk_cr.sum() > 0:
        M_total = M_total_arr[msk_cr][:N_obs_stars][-1]
    else:
        M_total = M_total_arr[:N_obs_stars][-1]

    # Percentage of mass added by the binaries:
    # binar_mass_perc = binary_sists_mass / single_sists_mass
    binar_mass_perc = isoch_binar[-1].sum() /\
        isoch_binar[m_ini_idx].sum()
    # Total mass corrected by the added mass as binary systems
    M_total = M_total * (1 + binar_mass_perc)

    # Assign errors according to errors distribution.
    synth_clust = add_errors.main(
        isoch_compl, err_lst, rand_norm_vals[1])

    # Transposing is necessary for np.histogramdd() in the
    # likelihood
    if transpose_flag:
        synth_clust = synth_clust[:sum(N_fc)].T

    return synth_clust, M_total


def properModel(fundam_params, model_proper, model, varIdxs):
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
      All fundamental parameters, including the fixed parameters that are
      missing from 'model'.
    ml, mh, al, ah : ints
      Indexes of the (z, a) values in the grid that define the box that enclose
      the proper (z, a) values.

    """
    model_proper[varIdxs] = model

    ml = mh = 0
    if 0 in varIdxs:
        par = fundam_params[0]
        mh = min(len(par) - 1, np.searchsorted(par, model_proper[0]))
        ml = mh - 1

    al = ah = 0
    if 1 in varIdxs:
        par = fundam_params[1]
        mh = min(len(par) - 1, np.searchsorted(par, model_proper[1]))
        al = ah - 1

    return model_proper, ml, mh, al, ah

    # model_proper, j = [], 0
    # for i, par in enumerate(fundam_params):
    #     # Check if this parameter is one of the 'free' parameters.
    #     if i in varIdxs:
    #         # If it is the parameter metallicity.
    #         if i == 0:
    #             # Select the closest value in the array of allowed values.
    #             mh = min(len(par) - 1, np.searchsorted(par, model[i - j]))
    #             ml = mh - 1
    #             # Define the model's z value
    #             z_model = model[i - j]
    #         # If it is the parameter log(age).
    #         elif i == 1:
    #             # Select the closest value in the array of allowed values.
    #             ah = min(len(par) - 1, np.searchsorted(par, model[i - j]))
    #             al = ah - 1
    #             a_model = model[i - j]
    #         else:
    #             model_proper.append(model[i - j])
    #     else:
    #         if i == 0:
    #             ml = mh = 0
    #             z_model = fundam_params[0][0]
    #         elif i == 1:
    #             al = ah = 0
    #             a_model = fundam_params[1][0]
    #         else:
    #             model_proper.append(par[0])
    #         j += 1

    # return model_proper, z_model, a_model, ml, mh, al, ah

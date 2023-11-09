
import numpy as np
from . import zaWAverage
from . import move_isochrone
from . import cut_max_mag
from . import mass_distribution
from . import mass_interp
from . import extinction
from . import binarity
from . import completeness_rm
from . import add_errors


# import matplotlib.pyplot as plt
# def plot(arr, name):
#     plt.title(name + f", N={len(arr[0])}")
#     print(name, arr[2].sum(), arr[-1].sum())
#     plt.scatter(arr[1], arr[0])
#     plt.gca().invert_yaxis()
#     plt.show()


def main(
    model, transpose_flag, DR_dist, DR_percentage, alpha, model_proper,
    varIdxs, completeness, err_lst, max_mag_syn, N_obs_stars, fundam_params,
    ext_coefs, N_fc, m_ini_idx, st_dist_mass, theor_tracks, rand_norm_vals,
        rand_unif_vals):
    """
    Takes an isochrone and returns a synthetic cluster created according to
    a certain mass distribution.

    The synthetic cluster returned has the shape:

    * transpose_flag = True

    synth_clust = [mag, c1, (c2)]

    where c1 and c2 colors defined.

    * transpose_flag = False

    synth_clust = [mag, c1, (c2), m_ini_1, mag_b, c1_b, (c2_b), m_ini_2]

    where 'm_ini_1, m_ini_2' are the primary and secondary masses of the
    binary systems. The single systems only have a '0' stored in 'm_ini_2'.
    """
    # max_mag_syn, N_obs_stars = 30, 1_000_000
    # for age in (8, 8.5, 9, 9.5):
    #     model_proper = np.array([0., age, 0.15, 0., 0., 3.1, 10.])
    #     model = np.array([age, 3.1])

    # Return proper values for fixed parameters and parameters required
    # for the (z, log(age)) isochrone averaging.
    met, age, beta, av, dr, rv, dm, ml, mh, al, ah = properModel(
        fundam_params, model_proper, model, varIdxs)

    # Generate a weighted average isochrone from the (z, log(age)) values in
    # the 'model'.
    isochrone = zaWAverage.main(
        theor_tracks, fundam_params, m_ini_idx, met, age, ml, mh, al, ah)
    # plot(isochrone, 'isoch')

    # Move theoretical isochrone using the distance modulus
    isoch_moved = move_isochrone.main(isochrone, N_fc, dm)
    # plot(isoch_moved, 'isoch_moved')

    # Remove stars beyond the maximum magnitude
    isoch_cut, _ = cut_max_mag.main(isoch_moved, max_mag_syn)
    # plot(isoch_cut, 'isoch_cut')

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
    mass_dist, M_total_arr = mass_distribution.main(mass_ini, st_dist_mass[ml])

    if not mass_dist.any():
        return synth_clust, M_total

    # Interpolate masses in mass_dist into the isochrone rejecting those
    # masses that fall outside of the isochrone's mass range.
    # This destroys the order by magnitude.
    isoch_mass = mass_interp.main(isoch_cut, mass_ini, mass_dist)
    # plot(isoch_mass, 'isoch_mass')

    # Assignment of binarity.
    isoch_binar = binarity.main(
        isoch_mass, alpha, beta, m_ini_idx, N_fc, rand_unif_vals[1])
    # plot(isoch_binar, 'isoch_binar')

    isoch_extin = extinction.main(
        isoch_binar, av, dr, rv, ext_coefs, N_fc, DR_dist, DR_percentage,
        rand_norm_vals[0], rand_unif_vals[0], m_ini_idx)
    # Remove stars moved beyond the maximum magnitude
    isoch_extin, msk_ct = cut_max_mag.main(isoch_extin, max_mag_syn)
    # plot(isoch_extin, 'isoch_extin')

    if not isoch_extin.any():
        return synth_clust, M_total

    # Completeness removal of stars.
    isoch_compl, msk_cr = completeness_rm.main(
        isoch_extin, completeness, rand_unif_vals[2])

    if not isoch_compl.any():
        return synth_clust, M_total

    # Keep the same number of synthetic stars as observed stars
    isoch_compl = isoch_compl[:, :N_obs_stars]
    # plot(isoch_compl, 'isoch_compl')

    # Apply 'msk_ct' and 'msk_cr' to correct the mass for the stars removed
    # by extinction() and completeness_rm() (if applied)
    if msk_cr.any():
        M_total = M_total_arr[msk_ct][msk_cr][:N_obs_stars][-1]
    else:
        M_total = M_total_arr[msk_ct][:N_obs_stars][-1]
    # print(M_total)

    # Correct M_total applying the percentage of mass added by the binaries:
    # binar_mass_perc = secondary_masses / primary_masses
    binar_mass_perc = isoch_binar[-1].sum() / isoch_binar[m_ini_idx].sum()
    M_total = M_total * (1 + binar_mass_perc)

    # import pandas as pd
    # dd = {}
    # for i, coln in enumerate(('Jmag', 'H-Ks', 'm1', 'Gmag_b', 'BP-RP_b', 'm2')):
    #     if i == 3 or i == 4:
    #         continue
    #     dd[coln]=isoch_compl[i]
    # df=pd.DataFrame(dd)
    # fname = str(age) + '_03_2000.csv'
    # df.to_csv("/home/gabriel/Descargas/" + fname,index=False)
    # print(age)
    # breakpoint()

    # Assign errors according to errors distribution.
    synth_clust = add_errors.main(isoch_compl, err_lst, rand_norm_vals[1])

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
    model_proper0  : list
      All fundamental parameters, including the fixed parameters that are
      missing from 'model'.
    ml, mh, al, ah : ints
      Indexes of the (z, a) values in the grid that define the box that enclose
      the proper (z, a) values.

    """
    model_proper0 = np.array(model_proper)
    model_proper0[varIdxs] = model

    ml = mh = 0
    if 0 in varIdxs:
        par = fundam_params[0]
        mh = min(len(par) - 1, np.searchsorted(par, model_proper0[0]))
        ml = mh - 1

    al = ah = 0
    if 1 in varIdxs:
        par = fundam_params[1]
        ah = min(len(par) - 1, np.searchsorted(par, model_proper0[1]))
        al = ah - 1

    return *model_proper0, ml, mh, al, ah

    # DEPRECATED 04/22
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

    # model_proper = [z_model, a_model] + model_proper

    # if (np.array(model_proper)==model_proper1).all() is False:
    #     breakpoint()
    # if ml1 != ml or mh1 != mh or al1 != al or ah1 != ah:
    #     breakpoint()

    return model_proper, ml, mh, al, ah


import copy
import itertools
import numpy as np
from . import obs_clust_prepare
from ..synth_clust import imf, synth_cluster, zaWAverage, move_isochrone,\
    cut_max_mag, mass_interp, completeness_rm
from .bf_common import varPars


def main(pd, clp, td):
    """
    Obtain several parameters needed for the best fit process / synthetic
    clusters generation.
    """

    # Number of stars in the observed cluster
    clp['N_obs_stars'] = len(clp['cl_reg_fit'])

    clp['max_mag_syn'] = np.max(list(zip(*list(zip(
        *clp['cl_reg_fit']))[1:][2]))[0])
    clp['cl_syn_fit'] = copy.deepcopy(clp['cl_reg_fit'])

    # Processed observed cluster.
    clp['obs_clust'] = obs_clust_prepare.main(
        clp['cl_syn_fit'], pd['lkl_method'], pd['lkl_binning'],
        pd['lkl_manual_bins'])

    return clp


def totMassEstim(clp, td, IMF_name, **kwargs):
    """
    Estimate the minimum mass necessary such that the combinations of the
    extreme values for the parameters (z, a, e, d)
    """
    fundam_params = td['fundam_params']
    varIdxs, _, _ = varPars(fundam_params)

    theor_tracks, m_ini_idx, ext_coefs, N_fc =\
    td['theor_tracks'], td['m_ini_idx'], td['ext_coefs'], td['N_fc'],\


    def getNmass(model, mass, ext_unif_rand):
        """
        """
        model_proper, z_model, a_model, ml, mh, al, ah =\
            synth_cluster.properModel(fundam_params, model, varIdxs)
        isochrone = zaWAverage.main(
            theor_tracks, fundam_params, m_ini_idx, z_model, a_model,
            ml, mh, al, ah)
        e, d, M_total, bin_frac, R_V = model_proper
        isoch_moved = move_isochrone.main(
            isochrone, e, d, R_V, ext_coefs, N_fc, ext_unif_rand,
            m_ini_idx)
        isoch_cut = cut_max_mag.main(isoch_moved, clp['max_mag_syn'])

        N_mass = 0
        if isoch_cut.any():
            mass_ini = isoch_cut[m_ini_idx]
            mmin, mmax = mass_ini.min(), mass_ini.max()
            msk_m = (mass >= mmin) & (mass <= mmax)
            N_mass = msk_m.sum()

        return N_mass

    # Combine the extreme values for (z, a, e, d)
    all_models = [
        (fundam_params[0][-1], fundam_params[0][0]),
        (fundam_params[1][-1], fundam_params[1][0]),
        (fundam_params[2][-1], fundam_params[2][0]),
        (fundam_params[3][-1], fundam_params[3][0])]
    models = list(itertools.product(*all_models))

    #
    inv_cdf = imf.invTrnsfSmpl(IMF_name)
    ext_unif_rand = np.random.uniform(0., 1., theor_tracks.shape[-1])

    Max_mass = clp['N_obs_stars']
    while True:
        mass = imf.sampleInv(Max_mass, inv_cdf)

        break_flag = True
        for model in models:
            model = list(model) + [0, 0]

            N_compl_rm_stars = 0
            if clp['completeness'][-1] is True:
                crp = testSynthClust(
                    model, clp['max_mag_syn'], clp['completeness'], varIdxs,
                    mass, ext_unif_rand, True, **td)
                N_compl_rm_stars = int(clp['N_obs_stars'] * crp)
            N_stars = clp['N_obs_stars'] + N_compl_rm_stars

            N_mass = getNmass(model, mass, ext_unif_rand)

            if N_mass < N_stars:
                break_flag = False
                print(model[:-2], round(Max_mass, 0), N_stars, N_mass)
                break
            else:
                print("great success!")
                print(model[:-2], round(Max_mass, 0), N_stars, N_mass)

        if break_flag is False:
            Max_mass = Max_mass * (1 + .1)
        else:
            break

    breakpoint()

    return Max_mass


def completenessPercEstim(pd, clp, td):
    """
    Simulate synthetic cluster to obtain an estimate of the percentage of stars
    that are removed by the 'completeness' function. This value is made to
    depend in the extinction and distance modulus parameter, found to be the
    two that most affect it.
    """
    fundam_params = td['fundam_params']
    varIdxs, _, _ = varPars(fundam_params)

    # Only extinction and distance modulus vary, the remaining parameters are
    # kept fixed.
    z0, a0, m0, b0 = fundam_params[0][0], fundam_params[1][0], 0, 0
    if len(fundam_params[2]) == 1:
        exts = [fundam_params[2][0]]
    else:
        exts = np.linspace(fundam_params[2][0], fundam_params[2][-1], 10)
    if len(fundam_params[3]) == 1:
        dms = [fundam_params[3][0]]
    else:
        dms = np.linspace(fundam_params[3][0], fundam_params[3][-1], 10)
    ed_compl_vals = [exts, dms]

    exts_dms = []
    for ext in exts:
        dms_lst = []
        for dm in dms:
            crp = testSynthClust(
                [z0, a0, ext, dm, m0, b0], clp['max_mag_syn'],
                clp['completeness'], varIdxs, td['st_dist_mass'],
                td['ext_unif_rand'], False, **td)
            dms_lst.append(crp)
        exts_dms.append(dms_lst)
    ed_compl_vals.append(np.array(exts_dms))

    return ed_compl_vals


def testSynthClust(
    model, max_mag_syn, completeness, varIdxs, st_dist_mass, ext_unif_rand,
    mass_estim_flag, fundam_params, theor_tracks,
    m_ini_idx, ext_coefs, N_fc, **kwargs):
    """
    Generate synthetic cluster and return either:
    1. the total number of stars
    generated previous to the removal of stars with the completeness function
     required percentage of
    stars removed by the 'completeness' function. All the stars in the
    'st_dist_mass' array are used.
    """
    model_proper, z_model, a_model, ml, mh, al, ah =\
        synth_cluster.properModel(fundam_params, model, varIdxs)
    isochrone = zaWAverage.main(
        theor_tracks, fundam_params, m_ini_idx, z_model,
        a_model, ml, mh, al, ah)
    e, d, _, _, R_V = model_proper
    isoch_moved = move_isochrone.main(
        isochrone, e, d, R_V, ext_coefs, N_fc,
        ext_unif_rand, m_ini_idx)
    isoch_cut = cut_max_mag.main(isoch_moved, max_mag_syn)

    if not isoch_cut.any():
        return 1

    mass_ini = isoch_cut[m_ini_idx]
    mmin, mmax = mass_ini.min(), mass_ini.max()
    if mass_estim_flag:
        mass = st_dist_mass
    else:
        mass = st_dist_mass[ml][0]
    msk_m = (mass >= mmin) & (mass <= mmax)
    mass_dist = mass[msk_m]

    isoch_mass = mass_interp.main(isoch_cut, mass_ini, mass_dist)
    if not isoch_mass.any():
        return 1
    isoch_compl = completeness_rm.main(isoch_mass, completeness)

    crp = 1 - len(isoch_compl[0]) / len(mass_dist)
    return crp

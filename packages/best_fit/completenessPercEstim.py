
import numpy as np
from ..synth_clust import synth_cluster, zaWAverage, move_isochrone,\
    cut_max_mag, mass_interp, completeness_rm


def main(clp, td):
    """
    Simulate synthetic cluster to obtain an estimate of the percentage of stars
    that are removed by the 'completeness' function. This value is made to
    depend on the age, extinction, and distance modulus parameter, found to be
    the three that most affect it.

    Used to correct the synthetic clusters in 'synth_cluster.py'.
    """
    fundam_params = td['fundam_params']
    # Indexes of "free" parameters (ext, dm)
    varIdxs = [1, 2, 3]

    # Only age, extinction, and distance modulus vary, the remaining parameters
    # are kept fixed.
    if len(fundam_params[1]) == 1:
        ages = [fundam_params[1][0]]
    else:
        ages = np.linspace(fundam_params[1][0], fundam_params[1][-1], 10)
    if len(fundam_params[2]) == 1:
        exts = [fundam_params[2][0]]
    else:
        exts = np.linspace(fundam_params[2][0], fundam_params[2][-1], 10)
    if len(fundam_params[3]) == 1:
        dms = [fundam_params[3][0]]
    else:
        dms = np.linspace(fundam_params[3][0], fundam_params[3][-1], 10)
    ed_compl_vals = [ages, exts, dms]

    age_exts_dms = []
    for age in ages:
        exts_dms = []
        for ext in exts:
            dms_lst = []
            for dm in dms:
                # Only *free* parameters go into the model
                model = [age, ext, dm]
                dms_lst.append(testSynthClust(
                    model, clp['max_mag_syn'], clp['completeness'], varIdxs,
                    **td))
            exts_dms.append(dms_lst)
        age_exts_dms.append(exts_dms)
    ed_compl_vals.append(np.array(age_exts_dms))

    return ed_compl_vals


def testSynthClust(
    model, max_mag_syn, completeness, varIdxs, st_dist_mass, ext_unif_rand,
        fundam_params, theor_tracks, m_ini_idx, ext_coefs, N_fc, **kwargs):
    """
    Generate synthetic cluster and return the percentage of stars removed
    by the 'completeness' function.
    """
    model_proper, z_model, a_model, ml, mh, al, ah =\
        synth_cluster.properModel(fundam_params, model, varIdxs)
    isochrone = zaWAverage.main(
        theor_tracks, fundam_params, m_ini_idx, z_model,
        a_model, ml, mh, al, ah)
    e, d, _, _, R_V = model_proper
    isoch_moved = move_isochrone.main(
        isochrone, e, d, R_V, ext_coefs, N_fc,
        ext_unif_rand[ml], m_ini_idx)
    isoch_cut = cut_max_mag.main(isoch_moved, max_mag_syn)

    if not isoch_cut.any():
        return 1

    mass_ini = isoch_cut[m_ini_idx]
    mmin, mmax = mass_ini.min(), mass_ini.max()
    mass = st_dist_mass[ml][0]
    msk_m = (mass >= mmin) & (mass <= mmax)
    mass_dist = mass[msk_m]

    isoch_mass = mass_interp.main(isoch_cut, mass_ini, mass_dist)
    if not isoch_mass.any():
        return 1
    isoch_compl = completeness_rm.main(isoch_mass, completeness)

    compl_rm_perc = 1 - len(isoch_compl[0]) / msk_m.sum()
    return compl_rm_perc

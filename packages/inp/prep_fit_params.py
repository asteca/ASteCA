
import copy
from . import read_met_files
from ..synth_clust import tracksPrep, extin_coefs
from ..best_fit import dataPrep


def main(pd, clp, clust_name, **kwargs):
    """
    Prepare tracks and other required data
    """
    td = {}
    if pd['best_fit_algor'] == 'n':
        return td

    clPresent = True
    try:
        td['fundam_params'] = copy.deepcopy(
            pd['fundam_params_all'][clust_name])
    except KeyError:
        clPresent = False
        print(("Cluster not found in line 'R3'. Default to "
               "'CLUSTER' values"))
    if clPresent is False:
        td['fundam_params'] = copy.deepcopy(
            pd['fundam_params_all']['CLUSTER'])

    clPresent = True
    try:
        td['priors_mcee'] = copy.deepcopy(
            pd['priors_mcee_all'][clust_name])
    except KeyError:
        clPresent = False
        print(("Cluster not found in line 'B2'. Default to "
               "'CLUSTER' values"))
    if clPresent is False:
        td['priors_mcee'] = copy.deepcopy(pd['priors_mcee_all']['CLUSTER'])

    # Store the number of defined filters and colors.
    td['N_fc'] = [len(pd['filters']), len(pd['colors'])]
    # Index of 'M_ini' (theoretical initial mass), stored in the
    # interpolated isochrones: right after the magnitude and color(s)
    td['m_ini_idx'] = td['N_fc'][0] + td['N_fc'][1]

    # Check and store metallicity files.
    # Keys added: 'isoch_list', 'extra_pars'
    td = read_met_files.check_get(pd, td)

    # Keys added: 'ext_coefs', 'N_fc', 'm_ini_idx',
    # 'st_dist_mass', 'theor_tracks', 'mean_bin_mr', 'err_norm_rand',
    # 'binar_probs', 'ext_unif_rand'
    td = tracksPrep.main(td, **pd)

    # Obtain extinction coefficients.
    td['ext_coefs'] = extin_coefs.main(
        pd['cmd_systs'], pd['filters'], pd['colors'])

    clp = dataPrep.main(pd, clp, td)

    td['Max_mass'] = dataPrep.totMassEstim(clp, td, **pd)

    clp['ed_compl_vals'] = dataPrep.completenessPercEstim(pd, clp, td)

    # Obtain mass distribution using the selected IMF.
    Nmets = len(td['fundam_params'][0])
    td['st_dist_mass'] = imf.main(IMF_name, Nmets, td['Max_mass'])

    td['err_norm_rand'], td['binar_probs'], td['ext_unif_rand'] = randVals(
        bp_vs_mass, td['st_dist_mass'], td['theor_tracks'])

    return td

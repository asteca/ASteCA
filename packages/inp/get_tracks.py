
from . import read_met_files
from ..synth_clust import tracksPrep


def main(pd, clust_name, **kwargs):
    """
    Prepare tracks and other required data
    """
    td = {}
    if pd['best_fit_algor'] != 'n':

        clPresent = True
        try:
            td['fundam_params'] = pd['fundam_params_all'][clust_name]
        except KeyError:
            clPresent = False
            print(("Cluster not found in line 'R1'. Default to "
                   "'CLUSTER' values"))
        if not clPresent:
            td['fundam_params'] = pd['fundam_params_all']['CLUSTER']

        clPresent = True
        try:
            td['priors_mcee'] = pd['priors_mcee_all'][clust_name]
        except KeyError:
            clPresent = False
            print(("Cluster not found in line 'B2'. Default to "
                   "'CLUSTER' values"))
        if not clPresent:
            td['priors_mcee'] = pd['priors_mcee_all']['CLUSTER']

        # Check and store metallicity files.
        # Keys added: 'isoch_list', 'extra_pars'
        td = read_met_files.check_get(pd, td)
        # Keys added: 'ext_coefs', 'binar_flag', 'N_fc', 'm_ini_idx',
        # 'st_dist_mass', 'theor_tracks', 'mean_bin_mr', 'err_norm_rand',
        # 'binar_probs', 'ext_unif_rand'
        td = tracksPrep.main(td, **pd)

    return td

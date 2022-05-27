
import copy
import numpy as np
from packages.inp import readZA
from packages.inp import read_isochs


def main(pd, clust_name, **kwargs):
    """
    Process all the metallicity files and the ages stored in them.
    """
    td = {}
    if pd['best_fit_algor'] == 'n':
        return td

    clPresent = True
    try:
        td['fundam_params'] = [
            pd['fundam_params_all'][clust_name][0],
            pd['fundam_params_all'][clust_name][1]]
        td['fundam_params'] += [
            np.array(_) for _ in pd['fundam_params_all'][clust_name][2:]]
    except KeyError:
        clPresent = False
        print(("Cluster not found in line 'R5'. Default to "
               "'CLUSTER' values"))
    if clPresent is False:
        td['fundam_params'] = [
            pd['fundam_params_all']['CLUSTER'][0],
            pd['fundam_params_all']['CLUSTER'][1]]
        td['fundam_params'] += [
            np.array(_) for _ in pd['fundam_params_all']['CLUSTER'][2:]]

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

    # Print info about tracks.
    nt = '' if len(pd['all_syst_filters']) == 1 else 's'
    print(
        "Processing {} isochrones in the photometric system{}:".format(
            pd['evol_track'], nt))
    for syst in pd['all_syst_filters']:
        print(" * {}".format(syst[0]))

    # Get all metallicity files and their values, and the log(age) values.
    met_files, met_vals_all, age_vals_all, ages_strs = readZA.main(
        td['fundam_params'], pd['iso_paths'], pd['evol_track'])

    # Store the common grid values for the metallicity and age.
    td['fundam_params'][:2] = met_vals_all, age_vals_all

    # Get isochrones and their extra parameters (mass, etc.).
    td['isoch_list'], td['extra_pars'] = read_isochs.main(
        met_files, ages_strs, pd['evol_track'], pd['CMD_extra_pars'],
        pd['all_syst_filters'])

    # Check equality of the initial mass across photometric systems.
    miniCheck(td['extra_pars'], met_vals_all, age_vals_all)

    # print("Grid values")
    # print("z        : {:<5} [{}, {}]".format(
    #     len(met_vals_all), td['fundam_params'][0][0],
    #     td['fundam_params'][0][-1]))
    # print("log(age) : {:<5} [{}, {}]".format(
    #     len(age_vals_all), td['fundam_params'][1][0],
    #     td['fundam_params'][1][-1]))

    return td


def miniCheck(extra_pars, met_vals_all, age_vals_all):
    """
    The extra isochrone parameter 'M_ini' is assumed to be equal across
    photometric systems, for a given metallicity and age. We check here that
    this is the case.

    extra_pars.shape = (#met_vals, #log(ages), #phot_systs)
    """
    extra_pars = np.array(extra_pars, dtype=object)
    # If a single z and log(age) are defined, this array will have a shape
    # (#met_vals, #log(ages), #phot_systs, #stars). Hence the '[:3]'.
    Nz, Na, Ndim = extra_pars.shape[:3]
    if Ndim == 1:
        # Single photometric system defined. Nothing to check.
        return
    else:
        txt = "initial mass values are not equal across the\n" +\
            "photometric systems for the isochrone: z={}, log(age)={}"
        for d in range(1, Ndim):
            for z in range(Nz):
                for a in range(Na):
                    arr0, arrd = extra_pars[z, a, 0], extra_pars[z, a, d]
                    if not np.array_equal(arr0, arrd):
                        # Check across (z, log(age)) values for each
                        # photometric system.
                        raise ValueError(
                            txt.format(met_vals_all[z], age_vals_all[a]))


import sys
from os.path import isdir, join


def check(mypath, pd):
    """
    Check all parameters related to the generation of synthetic clusters.

    Check that all the required photometric systems evolutionary tracks
    are present with checkEvolTracks()
    """

    # If best fit method is set to run.
    pd['fundam_params'] = []
    if pd['bf_flag'] or pd['best_fit_algor'] == 'synth_gen':

        # Check the random seed
        if pd['synth_rand_seed'] in ('n', 'N', 'None', 'none', 'no'):
            pd['synth_rand_seed'] = None
        else:
            try:
                pd['synth_rand_seed'] = int(pd['synth_rand_seed'])
            except ValueError:
                sys.exit("ERROR: random seed is not a valid integer")

        checkParamRanges(pd)

        pd = checkEvolTracks(mypath, pd)

        checkSynthClustParams(pd)

        pd = getParamVals(pd)

    return pd


def checkParamRanges(pd):
    """
    """
    # Check that no parameter range is empty.
    p_names = [
        'metallicity', 'age', 'extinction', 'distance', 'mass', 'binarity']
    for i, p in enumerate(pd['par_ranges']):
        if not p:
            sys.exit("ERROR: Range defined for '{}' parameter is"
                     " empty.".format(p_names[i]))
    m_range = pd['par_ranges'][4]
    if m_range[0] == 0.:
        print("  WARNING: minimum total mass is zero in params_input"
              " file.\n  Changed minimum mass to 10.\n")
        m_range[0] = 10.


def checkEvolTracks(mypath, pd):
    """
    Check photometric systems for all the evolutionary tracks that will be
    used. Store and return their IDs, filter names, and paths.
    """

    # Check selected isochrones set.
    if pd['evol_track'] not in pd['all_evol_tracks'].keys():
        sys.exit("ERROR: the selected isochrones set ('{}') does\n"
                 "not match a valid input.".format(pd['evol_track']))

    all_syst_filters, iso_paths = [], []
    # Remove duplicate filters (if they exist), and combine them into one
    # tuple per photometric system.
    # The resulting list looks like this:
    # [('2', 'T1', 'C'), ('4', 'B', 'V'), ('65', 'J')]
    # where the first element of each tuple points to the photometric
    # system, and the remaining elements are the unique filters in that
    # system.
    all_syst_filters = list(set(pd['filters'] + pd['c_filters']))
    d = {}
    for k, v in all_syst_filters:
        d.setdefault(k, [k]).append(v)
    all_syst_filters = sorted(map(tuple, d.values()))

    # Dictionary of photometric systems defined in the CMD service.
    all_systs = pd['cmd_systs']

    # Fix isochrones location according to the CMD and set selected.
    text1 = pd['all_evol_tracks'][pd['evol_track']][0]
    # Generate correct name for the isochrones path.
    iso_paths = []
    for p_syst in all_syst_filters:
        text2 = all_systs[p_syst[0]][0]
        # Set iso_path according to the above values.
        iso_paths.append(
            join(mypath + 'isochrones/' + text1 + '_' + text2))

    # Check if /isochrones folder exists.
    for iso_path in iso_paths:
        if not isdir(iso_path):
            sys.exit(
                "ERROR: 'Best synthetic cluster fit' function is set to"
                " run but the folder:\n\n {}\n\ndoes not exists."
                .format(iso_path))

    pd['all_syst_filters'], pd['iso_paths'] = all_syst_filters, iso_paths

    return pd


def checkSynthClustParams(pd):
    """
    """
    # Check IMF defined.
    if pd['IMF_name'] not in pd['imf_funcs']:
        sys.exit("ERROR: Name of IMF ({}) is incorrect.".format(
            pd['IMF_name']))

    if not 0. <= pd['bin_mr'] <= 1.:
        sys.exit("ERROR: Binary mass ratio set ('{}') is out of\n"
                 "boundaries. Please select a value in the range [0., 1.]".
                 format(pd['bin_mr']))

    # Check R_V defined.
    if pd['R_V'] <= 0.:
        sys.exit("ERROR: Ratio of total to selective absorption\n"
                 "R_V ({}) must be positive defined.".format(pd['R_V']))

    # Check maximum magnitude limit defined.
    if isinstance(pd['max_mag'], str):
        if pd['max_mag'] != 'max':
            sys.exit("ERROR: Maximum magnitude value selected ({}) is"
                     " not valid.".format(pd['max_mag']))


def getParamVals(pd):
    """
    Properly format parameter ranges to be used by the selected best fit
    method.
    """

    fundam_params = []
    for i, param in enumerate(pd['par_ranges']):
        # If only one value is defined.
        if len(param) == 1:
            fundam_params.append([param[0]])
        # If min == max store single value in array.
        elif param[0] == param[1]:
            fundam_params.append([param[0]])
        else:
            # Store range of values.
            fundam_params.append(param)

    # Check for min/max presence in z,log(age) parameters
    for idx in (0, 1):
        p_rng = fundam_params[idx]
        t = 'RZ' if idx == 0 else 'RA'
        if len(p_rng) == 1:
            try:
                p = float(p_rng[0])
            except ValueError:
                p = p_rng[0]
                if p not in ('min', 'max'):
                    sys.exit("ERROR '{}': unrecognized string '{}'".format(
                        t, p))
            fundam_params[idx] = [p]
        else:
            try:
                pmin = float(p_rng[0])
            except ValueError:
                pmin = p_rng[0]
                if pmin != 'min':
                    sys.exit(
                        ("ERROR '{}': unrecognized string '{}'.\nOnly 'min' " +
                         "string is accepted as the lower range.").format(
                            t, pmin))
            try:
                pmax = float(p_rng[1])
            except ValueError:
                pmax = p_rng[1]
                if pmax != 'max':
                    sys.exit(
                        ("ERROR '{}': unrecognized string '{}'.\nOnly 'max' " +
                         "string is accepted as the upper range.").format(
                            t, pmax))
            fundam_params[idx] = [pmin, pmax]

    pd['fundam_params'] = fundam_params
    return pd

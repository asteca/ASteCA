
import sys
from os.path import isdir, join


def check(mypath, pd):
    # par_ranges, lkl_method, lkl_methods, lkl_binning,
    # evol_track, max_mag, IMF_name, R_V, N_mass, bin_mr, bin_methods,
    # lkl_weight, bin_weights, cmd_evol_tracks, iso_paths, imf_funcs, hmax,
    #  inst_packgs_lst, ntemps,  nwalkers_ptm, nburn_ptm
    """
    Check all parameters related to the search for the best synthetic cluster
    match.
    """
    # If best fit method is set to run.
    if pd['bf_flag']:

        # TODO REMOVE when (if) support for multiple mags/colors is in place.
        if len(pd['filters']) > 1:
            sys.exit("ERROR: more than one filter defined.")
        if len(pd['colors']) > 2:
            sys.exit("ERROR: more than two colors defined.")

        checkParamRanges(pd)

        if pd['best_fit_algor'] == 'boot+GA':
            checkGA(pd)
        elif pd['best_fit_algor'] == 'ptemcee':
            checkPtemcee(pd)

        chechLkl(pd)

        pd = checkEvolTracks(mypath, pd)

        checkSynthClustParams(pd)

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


def checkGA(pd):
    """
    """

    if pd['hperc_btstrp'] < 0. or pd['hperc_btstrp'] > 0.9:
        sys.exit((
            "ERROR: 'pd['hperc_btstrp']' parameter ({}) must be in "
            "the [0, 0.9] range.").format(pd['hperc_btstrp']))

    # First set of params.
    oper_dict0 = {
        'n_pop': pd['N_pop'], 'n_gen': pd['N_gen'], 'n_el': pd['N_el'],
        'n_ei': pd['N_ei'], 'n_es': pd['N_es']}
    for oper in oper_dict0:
        if oper_dict0[oper] < 1:
            sys.exit("ERROR: number must be greater than zero in\n"
                     "'{}' GA parameter; {} is set.".format(
                         oper, oper_dict0[oper]))
    # Second set of params.
    oper_dict1 = {
        'fdif': pd['fit_diff'], 'p_cross': pd['cross_prob'],
        'p_mut': pd['mut_prob']}
    for oper in oper_dict1:
        if oper_dict1[oper] < 0. or oper_dict1[oper] > 1.:
            sys.exit("ERROR: GA '{}' parameter is out of the valid\n"
                     "[0., 1.] range; {} is set.".format(
                         oper, oper_dict1[oper]))
    # Handle separately.
    if pd['cross_sel'] not in ('1P', '2P'):
        sys.exit("ERROR: GA 'cr_sel' operator is not a valid choice;\n"
                 "'{}' is set.".format(pd['cross_sel']))
    # Number of solutions to pass to the nest generation (elitism)
    if pd['N_el'] >= pd['N_pop']:
        sys.exit("ERROR: GA 'n_el' must be smaller than 'n_pop';\n"
                 "'{}' and '{}' are set respectively.".format(
                     pd['N_el'], pd['N_pop']))


def checkPtemcee(pd):
    """
    """
    for pr in pd['priors_ptm']:
        if pr[0] not in pd['bayes_priors']:
            sys.exit("ERROR: one of the selected priors ({}) is not"
                     " allowed.".format(pr))

    if pd['nwalkers_ptm'] % 2 != 0:
        # Number is even
        sys.exit("ERROR: the number of walkers must be even.")
    if pd['nwalkers_ptm'] < 12:
        sys.exit("ERROR: the minimum number of walkers is 12.")
    if pd['nburn_ptm'] < 1:
        sys.exit("ERROR: the minimum number of burn-in samples is 1.")
    if pd['ntemps'] not in ('n', 'none', 'None'):
        if int(float(pd['ntemps'])) < 1:
            sys.exit("ERROR: the minimum number of temperatures is 1.")


def chechLkl(pd):
    """
    """
    # Check likelihood method selected.
    if pd['lkl_method'] not in pd['lkl_methods']:
        sys.exit("ERROR: the selected likelihood method '{}' does not"
                 " match a valid input.".format(pd['lkl_method']))

    # Check binning method selected.
    if pd['lkl_method'] != 'tolstoy' and pd['lkl_binning'] not in\
            pd['bin_methods']:
        sys.exit("ERROR: the selected binning method '{}' for the 'Best"
                 "\nfit' function does not match a valid input."
                 .format(pd['lkl_binning']))

    # Check binning weight method selected.
    if pd['lkl_method'] != 'tolstoy' and pd['lkl_weight'] not in\
            pd['bin_weights']:
        sys.exit("ERROR: the selected weight method '{}' for the 'Best"
                 "\nfit' function does not match a valid input."
                 .format(pd['lkl_weight']))

    # Check mass range selected.
    m_range = pd['par_ranges'][4]
    if pd['lkl_method'] == 'tolstoy':
        if len(m_range) > 1:
            sys.exit("ERROR: 'tolstoy' method was selected but more than"
                     "\none initial mass value is set.")


def checkEvolTracks(mypath, pd):
    """
    """
    # Check selected isochrones set.
    if pd['evol_track'] not in pd['cmd_evol_tracks'].keys():
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
    text1 = pd['cmd_evol_tracks'][pd['evol_track']][0]
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

    m_range = pd['par_ranges'][4]
    # Check mass range.
    try:
        if pd['N_mass'] > 100 and m_range[1] > 1e5:
            print("  WARNING: the number of masses defined is > 100 and\n"
                  "  the max mass is large ({:.0f}). This could cause\n"
                  "  memory issues when sampling the IMF.\n".format(
                      m_range[-1]))
    except IndexError:
        # Single value
        pass

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

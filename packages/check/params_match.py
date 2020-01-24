
import sys


def check(pd):
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

        if pd['best_fit_algor'] == 'boot+GA':
            checkGA(pd)
        if pd['best_fit_algor'] in ('ptemcee', 'emcee'):
            checkMcee(pd)
        if pd['best_fit_algor'] == 'ptemcee':
            checkPtemcee(pd)
        if pd['best_fit_algor'] == 'emcee':
            checkemcee(pd)

        chechLkl(pd)

    return pd


def checkGA(pd):
    """
    """

    if pd['hperc_btstrp'] < 0. or pd['hperc_btstrp'] > 0.9:
        sys.exit((
            "ERROR: 'pd['hperc_btstrp']' parameter ({}) must be in "
            "the [0, 0.9] range.").format(pd['hperc_btstrp']))

    # First set of params.
    oper_dict0 = {
        'n_pop': pd['N_pop'], 'n_el': pd['N_el'],
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


def checkMcee(pd):
    """
    """
    if pd['nwalkers_mcee'] % 2 != 0:
        # Number is even
        sys.exit("ERROR: the number of walkers must be even.")
    if pd['nwalkers_mcee'] < 12:
        sys.exit("ERROR: the minimum number of walkers is 12.")
    if pd['nburn_mcee'] <= 0. or pd['nburn_mcee'] >= 1:
        sys.exit("ERROR: burn-in percentage must be in the range (0., 1.)")

    for pr in pd['priors_mcee']:
        if pr[0] not in pd['bayes_priors']:
            sys.exit("ERROR: one of the selected priors ({}) is not"
                     " allowed.".format(pr))


def checkPtemcee(pd):
    """
    """
    if pd['pt_ntemps'] not in ('n', 'none', 'None'):
        if int(float(pd['pt_ntemps'])) < 1:
            sys.exit("ERROR: the minimum number of temperatures is 1.")

    try:
        float(pd['pt_tmax'])
    except ValueError:
        if pd['pt_tmax'] not in ('n', 'none', 'None', 'inf'):
            sys.exit("ERROR: 'Tmax' parameter is not a valid string.")


def checkemcee(pd):
    """
    """
    if 'emcee' not in pd['inst_packgs_lst']:
        sys.exit("ERROR: 'emcee' method is selected, but the package is\n" +
                 "not installed")


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

    if pd['lkl_binning'] == 'manual':
        for nbin in pd['lkl_manual_bins']:
            try:
                if int(nbin) < 2:
                    sys.exit("ERROR: the number of bins must be >2")
            except ValueError:
                sys.exit("ERROR: bin numbers must be integers")

        N_colors = int(len(pd['id_cols']) / 2.)
        if len(pd['lkl_manual_bins']) != N_colors + 1:
            sys.exit("ERROR: there are {} bin values defined,"
                     " there should be {}.".format(
                         len(pd['lkl_manual_bins']), N_colors + 1))

    # DEPRECATED 10/01/20
    # # Check binning weight method selected.
    # if pd['lkl_method'] != 'tolstoy' and pd['lkl_weight'] not in\
    #         pd['bin_weights']:
    #     sys.exit("ERROR: the selected weight method '{}' for the 'Best"
    #              "\nfit' function does not match a valid input."
    #              .format(pd['lkl_weight']))

    # Check mass range selected.
    m_range = pd['par_ranges'][4]
    if pd['lkl_method'] == 'tolstoy':
        if len(m_range) > 1:
            sys.exit("ERROR: 'tolstoy' method was selected but more than"
                     "\none initial mass value is set.")



def check(pd):
    """
    Check all parameters related to the search for the best synthetic cluster
    match.
    """

    # If best fit method is set to run.
    if pd['best_fit_algor'] != 'n':

        # TODO REMOVE when (if) support for multiple mags/colors is in place.
        if len(pd['filters']) > 1:
            raise ValueError("more than one filter defined.")
        if len(pd['colors']) > 2:
            raise ValueError("more than two colors defined.")

        checkPtemcee(pd)
        chechLkl(pd)

    return pd


def checkPtemcee(pd):
    """
    """
    if pd['nwalkers_mcee'] % 2 != 0:
        # Number is even
        raise ValueError("the number of walkers must be even.")
    if pd['nwalkers_mcee'] < 12:
        raise ValueError("the minimum number of walkers is 12.")
    if pd['nburn_mcee'] <= 0. or pd['nburn_mcee'] >= 1:
        raise ValueError("burn-in percentage must be in the range (0., 1.)")

    for pr in pd['priors_mcee']:
        if pr[0] not in pd['bayes_priors']:
            raise ValueError("one of the selected priors ({}) is not"
                             " allowed.".format(pr))

    if pd['pt_ntemps'] not in ('n', 'none', 'None'):
        if int(float(pd['pt_ntemps'])) < 1:
            raise ValueError("the minimum number of temperatures is 1.")

    try:
        float(pd['pt_tmax'])
    except ValueError:
        if pd['pt_tmax'] not in ('n', 'none', 'None', 'inf'):
            raise ValueError("'Tmax' parameter is not a valid string.")


def chechLkl(pd):
    """
    """
    # Check likelihood method selected.
    if pd['lkl_method'] not in pd['lkl_methods']:
        raise ValueError("the selected likelihood method '{}' does not"
                         " match a valid input.".format(pd['lkl_method']))

    # Check binning method selected.
    if pd['lkl_method'] != 'tolstoy' and pd['lkl_binning'] not in\
            pd['bin_methods']:
        raise ValueError(
            "the selected binning method '{}' for the 'Best\nfit' function"
            " does not match a valid input.".format(pd['lkl_binning']))

    if pd['lkl_binning'] == 'manual':
        for nbin in pd['lkl_manual_bins']:
            try:
                if int(nbin) < 2:
                    raise ValueError("the number of bins must be >2")
            except ValueError:
                raise ValueError("bin numbers must be integers")

        N_colors = int(len(pd['id_cols']))
        if len(pd['lkl_manual_bins']) != N_colors + 1:
            raise ValueError(
                "there are {} bin values defined, there should be {}.".format(
                    len(pd['lkl_manual_bins']), N_colors + 1))

    # DEPRECATED 10/01/20
    # # Check binning weight method selected.
    # if pd['lkl_method'] != 'tolstoy' and pd['lkl_weight'] not in\
    #         pd['bin_weights']:
    #     raise ValueError("the selected weight method '{}' for the 'Best"
    #              "\nfit' function does not match a valid input."
    #              .format(pd['lkl_weight']))

    # Check mass range selected.
    m_range = pd['par_ranges'][4]
    if pd['lkl_method'] == 'tolstoy':
        if len(m_range) > 1:
            raise ValueError("'tolstoy' method was selected but more than"
                             "\none initial mass value is set.")


# DEPRECATED May 2020
# def checkGA(pd):
#     """
#     """

#     if pd['hperc_btstrp'] < 0. or pd['hperc_btstrp'] > 0.9:
#         raise ValueError((
#             "'pd['hperc_btstrp']' parameter ({}) must be in "
#             "the [0, 0.9] range.").format(pd['hperc_btstrp']))

#     # First set of params.
#     oper_dict0 = {
#         'n_pop': pd['N_pop'], 'n_el': pd['N_el'],
#         'n_ei': pd['N_ei'], 'n_es': pd['N_es']}
#     for oper in oper_dict0:
#         if oper_dict0[oper] < 1:
#             raise ValueError(
#                 "number must be greater than zero in\n'{}' GA parameter; {} "
#                 "is set.".format(oper, oper_dict0[oper]))
#     # Second set of params.
#     oper_dict1 = {
#         'fdif': pd['fit_diff'], 'p_cross': pd['cross_prob'],
#         'p_mut': pd['mut_prob']}
#     for oper in oper_dict1:
#         if oper_dict1[oper] < 0. or oper_dict1[oper] > 1.:
#             raise ValueError(
#                 "GA '{}' parameter is out of the valid\n[0., 1.] range; {} is"
#                 " set.".format(oper, oper_dict1[oper]))
#     # Handle separately.
#     if pd['cross_sel'] not in ('1P', '2P'):
#         raise ValueError(
#             "GA 'cr_sel' operator is not a valid choice;\n'{}' is set.".format(
#                 pd['cross_sel']))
#     # Number of solutions to pass to the nest generation (elitism)
#     if pd['N_el'] >= pd['N_pop']:
#         raise ValueError(
#             "GA 'n_el' must be smaller than 'n_pop';\n'{}' and '{}' are set"
#             " respectively.".format(pd['N_el'], pd['N_pop']))

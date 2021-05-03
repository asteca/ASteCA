

def check(pd):
    """
    Check that kinematic data is properly defined.
    """

    if len(pd['id_kinem']) != 8:
        raise ValueError("there should be 8 entries in 'PK' line ({})".format(
            len(pd['id_kinem'])))

    # Read PMs, parallax, and RV data.
    k_cols = ('plx', 'e_plx', 'pmx', 'e_pmx', 'pmy', 'e_pmy', 'rv', 'e_rv')
    for i, ci in enumerate(pd['id_kinem']):
        if ci not in ('n', 'N'):
            try:
                pd[k_cols[i] + '_col'] = ci
            except ValueError:
                raise ValueError(
                    "bad index ('{}') for '{}' column in "
                    "'params_input.dat'".format(ci, k_cols[i]))
        else:
            pd[k_cols[i] + '_col'] = False

    # Check that PMs are either both or none defined.
    if (pd['pmx_col'] is False and pd['pmy_col'] is not False) or\
            (pd['pmy_col'] is False and pd['pmx_col'] is not False):
        raise ValueError("both (or none) PM dimensions must be defined"
                         " in 'params_input.dat'")

    # Check that error columns are present
    for col in ('plx_col', 'pmx_col', 'pmy_col', 'rv_col'):
        if pd[col] is not False and pd['e_' + col] is False:
            raise ValueError("missing error column for '{}' in"
                             "'params_input.dat'".format('e_' + col[:-4]))

    # # DEPRECATED 05/2021
    # if pd['plx_bayes_flag']:
    #     if 'emcee' not in pd['inst_packgs_lst']:
    #         raise ValueError("Plx data is set to be procesed, but 'emcee' is"
    #                          " not installed")

    if pd['plx_chains'] < 4:
        raise ValueError("set a minimum of 4 chains for Plx Bayesian analysis")
    if pd['plx_burn'] <= 0. or pd['plx_burn'] >= 1.:
        raise ValueError("Plx 'nburn' should be in the range (0., 1.)")
    # if pd['PM_KDE_std'] <= 0.:
    #     raise ValueError("'KDE_std' parameter must be greater than 0.")

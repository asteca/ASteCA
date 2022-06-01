

def check(pd):
    """
    Check that kinematic data is properly defined.
    """

    if len(pd['id_kinem']) != 6:
        raise ValueError("there should be 6 entries in 'I4' line ({})".format(
            len(pd['id_kinem'])))

    # Read PMs, parallax, and RV data.
    k_cols = ('plx', 'e_plx', 'pmx', 'e_pmx', 'pmy', 'e_pmy')
    for i, ci in enumerate(pd['id_kinem']):
        if ci not in ('n', 'N'):
            try:
                pd[k_cols[i] + '_col'] = ci
            except ValueError:
                raise ValueError(
                    "bad index ('{}') for '{}' column in "
                    "'asteca.ini'".format(ci, k_cols[i]))
        else:
            pd[k_cols[i] + '_col'] = False

    # Check that PMs are either both or none defined.
    if (pd['pmx_col'] is False and pd['pmy_col'] is not False) or\
            (pd['pmy_col'] is False and pd['pmx_col'] is not False):
        raise ValueError("both (or none) PM dimensions must be defined"
                         " in 'asteca.ini'")

    # Check that error columns are present
    for col in ('plx_col', 'pmx_col', 'pmy_col'):
        if pd[col] is not False and pd['e_' + col] is False:
            raise ValueError("missing error column for '{}' in"
                             "'asteca.ini'".format('e_' + col[:-4]))

    if pd['plx_chains'] < 4:
        raise ValueError("set a minimum of 4 chains for Plx Bayesian analysis")
    if pd['plx_burn'] <= 0. or pd['plx_burn'] >= 1.:
        raise ValueError("Plx 'nburn' should be in the range (0., 1.)")


# from os.path import isfile
# from packages.inp import names_paths


def check(
    cl_files, da_algor, bayesda_runs, bayesda_dflag,
    fld_rem_methods, bin_methods, fld_clean_mode, fld_clean_bin, colors,
        **kwargs):
    """
    Check parameters related to the decontamination algorithm functions.
    """

    if da_algor == 'y':
        # Check Bayesian decontamination algorithm parameters.
        if bayesda_runs < 2:
            raise ValueError("must input 'runs'>=2 for the Bayesian DA.")

        for bw in bayesda_dflag:
            if bw not in ('y', 'n'):
                raise ValueError("'bayes' DA flags must be either 'y' or 'n'")

        # This assumes that there is no maximum number of colors that can be
        # defined
        if 'y' not in bayesda_dflag:
            raise ValueError(
                "at least one 'bayes' DA weight must be set to 'y'.")
        if len(bayesda_dflag) - 4 != len(colors):  # HARDCODED '4'
            raise ValueError(
                "there are {} 'bayes' DA weights defined, there should "
                "be {}.".format(len(bayesda_dflag), 4 + len(colors)))

    # Check 'field stars removal' method selected.
    if fld_clean_mode not in fld_rem_methods:
        raise ValueError(
            "the selected field stars removal method ('{}') does\nnot match "
            "a valid input.".format(fld_clean_mode))
    # Check binning if 'local' method was selected.
    if fld_clean_mode == 'local' and fld_clean_bin not in bin_methods:
        raise ValueError(
            "the selected binning method '{}' for the 'Reduced"
            "\nmembership' function does not match a valid input.".format(
                fld_clean_bin))

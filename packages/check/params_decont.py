
from os.path import isfile
from packages.inp import names_paths


def check(cl_files, da_algor, da_algors_accpt, bayesda_runs, fixedda_port,
          bayesda_dflag, fld_rem_methods, bin_methods, fld_clean_mode,
          fld_clean_bin, colors, **kwargs):
    """
    Check parameters related to the decontamination algorithm functions.
    """

    # Check selected algorithm.
    if da_algor not in da_algors_accpt:
        raise ValueError("the selected decontamination algorithm ({})\n"
                         "is not recognized.".format(da_algor))

    if da_algor == 'bayes':
        # Check Bayesian decontamination algorithm parameters.
        if bayesda_runs < 2:
            raise ValueError("must input 'runs'>=2 for the Bayesian DA.")

        for bw in bayesda_dflag:
            if bw not in ('y', 'n'):
                raise ValueError("'bayes' DA flags must be either 'y' or 'n'")

        # This assumes that there is no maximum number of colors that can be
        # defined
        if len(bayesda_dflag) - 5 != len(colors):
            raise ValueError(
                "there are {} 'bayes' DA weights defined, there should "
                "be {}.".format(len(bayesda_dflag), 5 + len(colors)))

    if da_algor == 'fixed':
        # Check Bayesian decontamination algorithm parameters.
        if fixedda_port > 1 or fixedda_port < 0:
            raise ValueError("'fport' can not be outside the [0, 1] range.")

    # 'Read' mode is set.
    if da_algor == 'read':
        # Check if file exists.
        for cl_file in cl_files:
            # Get file name for membership files.
            memb_file = names_paths.memb_file_name(cl_file)
            if not isfile(memb_file):
                # File does not exist.
                raise ValueError(
                    "'read' mode was set for decontamination algorithm but "
                    "the file:\n\n {}\n\ndoes not exist.".format(memb_file))

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

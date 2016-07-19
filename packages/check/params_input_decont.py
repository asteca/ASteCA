
import sys
from os.path import isfile
import packages.inp.input_params as g
from packages.inp import names_paths


def check(cl_files, bin_methods_dict):
    """
    Check parameters related to the decontamination algorithm functions.
    """

    # Check decontamination algorithm params.
    mode_da = g.da_params[0]
    # Check if 'mode' was correctly set.
    if mode_da not in ['auto', 'manual', 'read', 'skip']:
        sys.exit("ERROR: Wrong name ('{}') for decontamination algorithm "
                 "'mode'.".format(mode_da))

    if mode_da == 'read':
        # Check if file exists.
        for cl_file in cl_files:

            # Get memb file names.
            memb_file = names_paths.main(cl_file)[2]
            if not isfile(memb_file):
                # File does not exist.
                sys.exit("ERROR: 'read' mode was set for decontamination "
                         "algorithm but the file:\n\n {}\n\ndoes not "
                         "exist.".format(memb_file))

    # Check 'Reduced membership' method selected.
    mode_red_memb, local_bin, min_prob = g.rm_params
    if mode_red_memb not in {'local', 'n_memb', 'mp_05', 'top_h', 'man', 'mag',
                             'skip'}:
        sys.exit("ERROR: the selected reduced membership method ('{}')"
                 " does\nnot match a valid input.".format(mode_red_memb))
    # Check binning if 'local' method was selected.
    if mode_red_memb == 'local' and local_bin not in bin_methods_dict:
        sys.exit("ERROR: the selected binning method '{}' for the 'Reduced"
                 "\nmembership' function does not match a valid input."
                 .format(local_bin))

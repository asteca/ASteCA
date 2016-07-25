
import sys
from os.path import join, isfile
import packages.inp.input_params as g


def check(mypath, cl_files):
    """
    Check that the parameters are properly written.
    """

    # Check mode.
    if g.mode not in {'auto', 'semi', 'manual'}:
        sys.exit("ERROR: 'mode' value selected ('{}') is not valid.".format(
            g.mode))

    if g.mode == 'semi':
        # Check if semi_input.dat file exists.
        semi_file = 'semi_input.dat'
        if not isfile(join(mypath, semi_file)):
            # File semi_input.dat does not exist.
            sys.exit("ERROR: 'semi' mode is set but semi_input.dat file does"
                     " not exist.")

    # Check px/deg.
    if g.gd_params[-1] not in {'px', 'deg'}:
        sys.exit("ERROR: the coordinate units given in the input parameters\n"
                 "file ('{}') are incorrect.".format(g.gd_params[-1]))

    # Selected CMD.
    if g.ps_params[1] not in {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13}:
        sys.exit("ERROR: the stored CMD value ({}) does not match a valid"
                 " selection.".format(g.ps_params[1]))

    # Output figure.
    if g.pl_params[0] is True:
        if g.pl_params[1] not in {'png', 'pdf', 'PNG', 'PDF'}:
            sys.exit("ERROR: figure output format selected ('{}') is"
                     " not valid.".format(g.pl_params[1]))

    # 2D positional histogram.
    if g.gh_params[0] not in {'auto', 'manual'}:
        sys.exit("ERROR: mode selected ('{}') for 2D histogram"
                 " is not valid.".format(g.gh_params[0]))

    # Radius finding function.
    if g.cr_params[0] not in {'low', 'mid', 'high'}:
        sys.exit("ERROR: mode selected ('{}') for radius finding"
                 " function is not valid.".format(g.cr_params[0]))

    # Errors function.
    if g.er_params[0] not in {'emax', 'lowexp', 'eyefit'}:
        sys.exit("ERROR: mode selected ('{}') for error rejecting"
                 " function is not valid.".format(g.er_params[0]))
    if g.er_params[0] == 'emax' and len(g.er_params[1:]) < 1:
        sys.exit("ERROR: missing parameters for error rejecting function")
    if g.er_params[0] == 'eyefit' and len(g.er_params[1:]) < 3:
        sys.exit("ERROR: missing parameters for error rejecting function")
    if g.er_params[0] == 'lowexp' and len(g.er_params[1:]) < 4:
        sys.exit("ERROR: missing parameters for error rejecting function")

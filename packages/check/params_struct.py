
import sys
from os.path import join, isfile


def check(mypath, cl_files, mode, coords, pl_params, gh_params,
          cr_params, er_params, **kwargs):
    """
    Check that the parameters are properly written.
    """

    # Check mode.
    if mode not in ('auto', 'semi', 'manual'):
        sys.exit("ERROR: 'mode' value selected ('{}') is not valid.".format(
            mode))

    if mode == 'semi':
        # Check if semi_input.dat file exists.
        semi_file = 'semi_input.dat'
        if not isfile(join(mypath, semi_file)):
            # File semi_input.dat does not exist.
            sys.exit("ERROR: 'semi' mode is set but 'semi_input.dat' file does"
                     " not exist.")

    # Check px/deg.
    if coords not in ('px', 'deg'):
        sys.exit("ERROR: coordinate units '{}' given in the input parameters\n"
                 "file are incorrect.".format(coords))

    # Output figure.
    if pl_params[0]:
        if pl_params[1] not in ('png', 'pdf', 'PNG', 'PDF'):
            sys.exit("ERROR: figure output format selected ('{}') is"
                     " not valid.".format(pl_params[1]))

    # 2D positional histogram.
    if gh_params[0] not in ('auto', 'manual'):
        sys.exit("ERROR: mode selected ('{}') for 2D histogram"
                 " is not valid.".format(gh_params[0]))

    # Radius finding function.
    if cr_params[0] not in ('low', 'mid', 'high'):
        sys.exit("ERROR: mode selected ('{}') for radius finding"
                 " function is not valid.".format(cr_params[0]))

    # Errors function.
    if er_params[0] not in ('emax', 'lowexp', 'eyefit'):
        sys.exit("ERROR: mode selected ('{}') for error rejecting"
                 " function is not valid.".format(er_params[0]))
    if er_params[0] == 'emax' and len(er_params[1:]) < 1:
        sys.exit("ERROR: missing parameters for error rejecting function")
    if er_params[0] == 'eyefit' and len(er_params[1:]) < 3:
        sys.exit("ERROR: missing parameters for error rejecting function")
    if er_params[0] == 'lowexp' and len(er_params[1:]) < 4:
        sys.exit("ERROR: missing parameters for error rejecting function")

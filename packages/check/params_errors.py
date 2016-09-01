
import sys


def check(er_params, **kwargs):
    """
    Check that the error parameters are properly written.
    """
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

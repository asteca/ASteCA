
import sys


def check(center_stddev, radius_method, **kwargs):
    """
    Check that the parameters are properly written.
    """

    if center_stddev <= 0.:
        sys.exit("ERROR: standard deviation value ('{}') must be greater\n"
                 "than zero.".format(center_stddev))

    # Radius finding function.
    if radius_method not in ('low', 'mid', 'high'):
        sys.exit("ERROR: mode selected ('{}') for radius finding"
                 " function is not valid.".format(radius_method))

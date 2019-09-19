
import sys
from collections import Counter


def check(manual_struct, center_bw, radius_method, fdens_method, **kwargs):
    """
    Check that the parameters are properly written.
    """

    dups = [
        _ for _, c in Counter(list(zip(*manual_struct))[0]).items() if c > 1]
    if dups:
        sys.exit((
            "ERROR: duplicated entries found in 'Structure' block:\n" +
            "{}".format(dups)))

    if center_bw < 0.:
        sys.exit("ERROR: KDE bandwidth ('{}') must be greater\n"
                 "than (or equal to) zero.".format(center_bw))

    # Radius finding function.
    if radius_method not in ('low', 'mid', 'high'):
        sys.exit("ERROR: mode selected ('{}') for radius finding"
                 " function is not valid.".format(radius_method))

    # Radius finding function.
    try:
        fd = float(fdens_method)
        if fd < 0.:
            sys.exit("ERROR: field density ('{}') must be"
                     " greater than zero.".format(fd))
    except ValueError:
        if fdens_method != 'auto':
            sys.exit("ERROR: field density mode ('{}') is not"
                     " recognized.".format(fdens_method))


from collections import Counter


def check(manual_struct, center_bw, fdens_method, RDP_rings, **kwargs):
    """
    Check that the parameters are properly written.
    """

    dups = [
        _ for _, c in Counter(list(zip(*manual_struct))[0]).items() if c > 1]
    if dups:
        raise ValueError((
            "duplicated entries found in 'Structure' block:\n" +
            "{}".format(dups)))

    if center_bw < 0.:
        raise ValueError("KDE bandwidth ('{}') must be greater\n"
                         "than (or equal to) zero.".format(center_bw))

    # Radius finding function.
    try:
        fd = float(fdens_method)
        if fd < 0.:
            raise ValueError("field density ('{}') must be"
                             " greater than zero.".format(fd))
    except ValueError:
        if fdens_method[-1] == '%':
            fdens = float(fdens_method[:-1])
            if fdens > 100. or fdens <= 0.:
                raise ValueError(
                    "percentage value in field density must be\n"
                    "in the (0., 100] range")
        elif fdens_method not in ('min', 'last', 'iter'):
            raise ValueError("field density mode ('{}') is not"
                             " recognized.".format(fdens_method))

    if RDP_rings < 5 or RDP_rings > 500:
        raise ValueError("valid range is 5 < RDP_rings < 500")

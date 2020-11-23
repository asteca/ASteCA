
from collections import Counter


def check(
    trim_frame_range, manual_struct, center_bw, fdens_method, kp_ndim,
        kp_nchains, kp_nburn, inst_packgs_lst, **kwargs):
    """
    Check that the parameters are properly written.
    """

    dups = [
        _ for _, c in Counter(list(zip(*trim_frame_range))[0]).items()
        if c > 1]
    if dups:
        raise ValueError((
            "duplicated entries found in 'Trim frame' block:\n" +
            "{}".format(dups)))

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

    if kp_ndim not in (0, 2, 4):
        raise ValueError(
            "Unrecognized value for King profile 'ndim' parameter")
    elif kp_ndim in (2, 4):
        if 'emcee' not in inst_packgs_lst:
            raise ValueError("King profile is selected to run, but 'emcee' is"
                             " not installed")
        if kp_nchains < 10:
            raise ValueError(
                "set a minimum of 10 chains for KP Bayesian analysis")
        if kp_nburn <= 0. or kp_nburn >= 1.:
            raise ValueError("KP 'nburn' should be in the range (0., 1.)")

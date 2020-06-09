

def check(pd):
    """
    Check that the parameters are properly written.
    """

    # Output figure.
    if pd['flag_make_plot']:
        for _ in pd['flag_make_plot']:
            if _ not in [
                    'A1', 'A2', 'A3', 'B1', 'B2', 'B3', 'C1', 'C2', 'C3', 'D0',
                    'D1', 'D2', 'D3', 's']:
                raise ValueError(
                    "unrecognized block ('{}') selected for plotting.".format(
                        _))

    pd['stop_idx'] = 'no'
    if 's' in pd['flag_make_plot']:
        if pd['flag_make_plot'].index('s') > 0:
            stop_idx = pd['flag_make_plot'].index('s')
            pd['stop_idx'] = pd['flag_make_plot'][stop_idx - 1]

    return pd

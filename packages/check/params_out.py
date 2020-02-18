

def check(pd):
    """
    Check that the parameters are properly written.
    """

    # Output figure.
    if pd['flag_make_plot']:
        for _ in pd['flag_make_plot']:
            if _ not in [
                    'A1', 'A2', 'B1', 'B2', 'C1', 'C2', 'C3', 'D0', 'D1', 'D2',
                    'D3', 's']:
                raise ValueError(
                    "unrecognized block ('{}') selected for plotting.".format(
                        _))

    if pd['plot_frmt'] not in ('png', 'pdf', 'PNG', 'PDF'):
        raise ValueError("figure output format selected ('{}') is"
                         " not valid.".format(pd['plot_frmt']))

    pd['stop_idx'] = 'no'
    if 's' in pd['flag_make_plot']:
        if pd['flag_make_plot'].index('s') > 0:
            stop_idx = pd['flag_make_plot'].index('s')
            pd['stop_idx'] = pd['flag_make_plot'][stop_idx - 1]

    return pd

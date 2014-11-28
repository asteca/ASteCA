from math import log10, floor


def round_to_y(x, y=1):
    '''
    Rounds the float x according to the number of decimal places given
    by the integer y (y=1 if not given).

    >>> round_to_y(31.51)
    30.0
    >>> round_to_y(0.0035)
    0.004
    >>> round_sig(31.56, 5)
    31.56
    >>> round_sig(0.0034512, 6)
    0.0034512
    '''
    x_r = 0.0 if x == 0. else round(x, y - int(floor(log10(abs(x)))) - 1)

    return x_r


def round_to_ref(x, y):
    '''
    Takes a float x and rounds it according to the number of decimal
    figures present in the float y.

    >>> round_to_ref(31.56, 0.002)
    31.56
    >>> round_to_ref(31.56, 0.1)
    31.6
    >>> round_to_ref(0.0035, 0.0005)
    0.0035
    >>> round_to_ref(0.0035, 0.05)
    0.0
    '''
    x_r = x if y == 0. else round(x, -int(floor(log10(abs(y)))))

    return x_r


def round_sig_fig(params, errors):
    '''
    Round errors to 1 significant figure and parameter values to the
    corresponding number given by the length of each error.
    '''

    # Round errors to 1 significant figure.
    errors_r = map(round_to_y, errors)

    # Round parameters to the number of significant figures in the erroi.
    params_r = []
    for i, p in enumerate(params):
        if errors_r[i] > 0.:
            params_r.append(round_to_ref(p, errors_r[i]))
        else:
            # In no error was assigned (i.e.: errors_r[i]=-1), don't attempt
            # to call the rounding function. Instead simply write the value
            # using the 'g' string format.
            params_r.append(float('{:g}'.format(p)))

    return params_r, errors_r

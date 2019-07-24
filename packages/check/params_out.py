
import sys


def check(flag_make_plot, plot_frmt, **kwargs):
    """
    Check that the parameters are properly written.
    """

    # Output figure.
    if flag_make_plot:
        for _ in flag_make_plot:
            if _ not in [
                    'A1', 'A2', 'B1', 'B2', 'C1', 'C2', 'C3', 'D1', 'D2',
                    'D3']:
                sys.exit("ERROR: unrecognized block ('{}') selected for"
                         " plotting.".format(_))

    if plot_frmt not in ('png', 'pdf', 'PNG', 'PDF'):
        sys.exit("ERROR: figure output format selected ('{}') is"
                 " not valid.".format(plot_frmt))


import os
import platform
from subprocess import Popen, PIPE


def check(pd):
    """
    Check that the parameters are properly written.
    """

    # Output figure.
    if pd['flag_make_plot']:
        for _ in pd['flag_make_plot']:
            if _ not in pd['plots_names']:
                raise ValueError(
                    "unrecognized block ('{}') selected for plotting.".format(
                        _))

    pd['stop_idx'] = 'no'
    if 's' in pd['flag_make_plot']:
        if pd['flag_make_plot'].index('s') > 0:
            stop_idx = pd['flag_make_plot'].index('s')
            pd['stop_idx'] = pd['flag_make_plot'][stop_idx - 1]

    import matplotlib
    import matplotlib.pyplot as plt

    # Force matplotlib to not use Xwindows backend. This call prevents
    # the code from crashing when used in a cluster. See:
    # http://stackoverflow.com/a/3054314/1391441
    if not X_is_running():
        matplotlib.use('Agg')
        print("(Force matplotlib to not use Xwindows backend)\n")

    # Define the style to use in all the plots.
    if pd['plot_style'] == 'asteca':
        combined_style = [os.getcwd() + '/packages/out/asteca.mplstyle']
    else:
        combined_style = [pd['plot_style']]
    if os.getcwd() in matplotlib.matplotlib_fname():
        # Combine style with 'matplotlibrc' file
        print("Custom 'matplotlibrc' file found\n")
        combined_style += ['matplotlibrc']
    plt.style.use(combined_style)

    return pd


def X_is_running():
    """
    Detect if X11 is available. Source:
    https://stackoverflow.com/a/1027942/1391441
    """
    if platform.system() == 'Linux':
        p = Popen(["xset", "-q"], stdout=PIPE, stderr=PIPE)
        p.communicate()
        return p.returncode == 0
    else:
        # If this is not a Linux system, assume that it is either Mac OS or
        # Windows, and thus assume that a windows system is present.
        return True


from subprocess import check_output
import datetime
import matplotlib.pyplot as plt
from .._version import __version__


def main(x_fix=0.02, y_fix=.996):
    """
    Add version number to top left of plot.
    """

    commit, Nc = '', 0
    if __version__[-3:] == 'dev':
        try:
            logs = str(check_output(['git', 'log', '-1']).decode("UTF-8"))
            Nc = 12
            commit = ' (' + logs[7:Nc] + ')'
        except Exception:
            # No git repository or some error
            pass

    ver_commit = __version__ + commit
    ver = '[ASteCA ' + ver_commit + '] ' + datetime.datetime.now().strftime(
        "%Y-%m-%d %H:%M")
    x_coord = x_fix - (len(ver_commit) - (6 + Nc)) * 0.001
    plt.figtext(x_coord, y_fix, ver, fontsize=9, color='#585858')

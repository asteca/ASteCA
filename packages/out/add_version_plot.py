
import datetime
import matplotlib.pyplot as plt
from .._version import __version__


def main(x_fix=0.02, y_fix=.996):
    """
    Add version number to top left of plot.
    """
    ver = '[ASteCA ' + __version__ + '] ' + datetime.datetime.now().strftime(
        "%Y-%m-%d %H:%M")
    x_coord = x_fix - (len(__version__) - 6) * 0.001
    plt.figtext(x_coord, y_fix, ver, fontsize=9, color='#585858')

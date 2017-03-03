
import matplotlib.pyplot as plt
from .._version import __version__


def main(x_fix=0.957):
    """
    Add version number to top left of plot.
    """
    ver = '[ASteCA ' + __version__ + ']'
    x_coord = x_fix - (len(__version__) - 6) * 0.001
    plt.figtext(x_coord, .993, ver, fontsize=9, color='#585858')

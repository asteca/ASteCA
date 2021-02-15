
import numpy as np


def exp_3p(x, a, b, c):
    """
    Three-parameters exponential function.

    This function is tied to the 'synth_cluster.add_errors' function.
    """
    return a * np.exp(b * x) + c


def exp_2p(x, a, b):
    """
    Two-parameters exponential function.
    """
    return a * np.exp(b * x)

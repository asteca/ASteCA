
import numpy as np


def main(synth_rand_seed):
    """
    Random seed used in the generation of synthetic clusters. Used by:

    * binary
    * imf
    * add_errors
    """
    np.random.seed(synth_rand_seed)

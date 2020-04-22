
import numpy as np


def main(synth_rand_seed):
    """
    Random seed used in the generation of synthetic clusters. Used by:

    * binary
    * imf
    * add_errors
    """
    if synth_rand_seed is None:
        synth_rand_seed = np.random.randint(100000)
        print("Random seed for synthetic clusters (random): {}".format(
            synth_rand_seed))
    else:
        print("Random seed for synthetic clusters (manual): {}".format(
            synth_rand_seed))

    np.random.seed(synth_rand_seed)

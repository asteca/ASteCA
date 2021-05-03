
import numpy as np


def main(synth_rand_seed):
    """
    Random seed used in the generation of synthetic clusters. Used by:

    * imf
    * binary
    * move_isochrone (In place for #174) & add_errors use the random arrays
      generated in tracksPrep
    """
    if synth_rand_seed is None:
        synth_rand_seed = np.random.randint(100000)
        print("\nRandom seed for synthetic clusters (random): {}".format(
            synth_rand_seed))
    else:
        print("\nRandom seed for synthetic clusters (manual): {}".format(
            synth_rand_seed))

    # In place for #505
    # numpy_RNG = np.random.default_rng(seed)
    # return numpy_RNG

    np.random.seed(synth_rand_seed)

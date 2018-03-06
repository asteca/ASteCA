
import numpy as np


def main(isoch_moved, max_mag_syn):
    '''
    Remove stars from isochrone with magnitude values larger that the maximum
    value found in the observation (entire field, not just the cluster
    region).

    This assumes that the entire array is *already* sorted (min, max) according
    to the main magnitude.
    '''

    # Get index of closest mag value to max_mag_syn. The array *must* be
    # ordered for this to work.
    max_indx = np.searchsorted(isoch_moved[0], max_mag_syn)
    # Discard elements beyond max_mag_syn limit.
    isoch_cut = isoch_moved[:, :max_indx]

    return isoch_cut

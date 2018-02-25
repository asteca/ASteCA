
import numpy as np


def main(isoch_moved, max_mag_syn):
    '''
    Remove stars from isochrone with magnitude values larger that the maximum
    value found in the observation (entire field, not just the cluster
    region).

    This assumes that the entire array is *already* sorted (min, max) according
    to the main magnitude.
    '''

    # Get index of closest mag value to max_mag_syn.
    max_indx = np.argmin(np.abs(isoch_moved[0] - max_mag_syn))
    # Discard elements beyond max_mag_syn limit.
    isoch_cut = np.array([d[0:max_indx] for d in isoch_moved])

    return isoch_cut

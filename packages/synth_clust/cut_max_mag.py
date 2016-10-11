
import numpy as np


def main(isoch_moved, max_mag):
    '''
    Remove stars from isochrone with magnitude values larger that the maximum
    value found in the observation (entire field, not just the cluster
    region).
    '''
    # Sort isochrone according to magnitude values (min to max).
    # with timeblock(" cut1"):
    isoch_sort = zip(*sorted(zip(*isoch_moved), key=lambda x: x[0]))
    # Get index of closest mag value to max_mag.
    # with timeblock(" cut2"):
    max_indx = min(range(len(isoch_sort[0])), key=lambda i:
                   abs(isoch_sort[0][i] - max_mag))
    # Discard elements beyond max_mag limit.
    # with timeblock(" cut3"):
    isoch_cut = np.array([d[0:max_indx] for d in isoch_sort])

    return isoch_cut


import numpy as np


def find_closest(key, target):
    '''
    See: http://stackoverflow.com/a/8929827/1391441
    Helper function for locating the mass values in the IMF distribution onto
    the isochrone.
    General: find closest target element for elements in key.

    Returns an array of indexes of the same length as 'target'.
    '''

    # Find indexes where the masses in 'target' should be located in 'key' such
    # that if the masses in 'target' were inserted *before* these indices, the
    # order of 'key' would be preserved. I.e.: pair masses in 'target' with the
    # closest masses in 'key'.
    # key must be sorted in ascending order.
    idx = key.searchsorted(target)

    # Convert indexes in the limits (both left and right) smaller than 1 and
    # larger than 'len(key) - 1' to 1 and 'len(key) - 1', respectively.
    idx = np.clip(idx, 1, len(key) - 1)

    # left < target <= right for each element in 'target'.
    left = key[idx - 1]
    right = key[idx]

    # target - left < right - target is True (or 1) when target is closer to
    # left and False (or 0) when target is closer to right.
    # The indexes stored in 'idx' point, for every element (IMF mass) in
    # 'target' to the closest element (isochrone mass) in 'key'. Thus:
    # target[XX] <closest> key[idx[XX]]
    idx -= target - left < right - target

    return idx


def main(isoch_cut, mass_dist, m_ini):
    '''
    For each mass in the IMF mass distribution, find the star in the isochrone
    with the closest mass value and pass it forward.
    Masses that fall outside of the isochrone's mass range are rejected.
    '''
    # Returns the indices that would sort the isochrone with
    # the minimum mass. This is why we use the 'm_ini' index.
    order = isoch_cut[m_ini, :].argsort()

    # Returns an array with the mass values in the theoretical isochrone
    # ordered from min to max.
    key = isoch_cut[m_ini, order]

    ##########################################################################
    # # Uncomment this block to see how many stars are being discarded
    # # for having mass values outside of the isochrone's range.
    # # Also uncomment the first block in the main function that uses the
    # # entire isochrone (instead of one with a max-mag cut).
    # try:
    #     print 'Min, max mass values in isochrone: {:.3f}, {:.3f}'.format(
    #         min(isoch_cut[m_ini]), max(isoch_cut[m_ini]))
    #     print 'Total mass in mass_dist: {:.2f}'.format(sum(mass_dist))
    #     # Masses out of boundary to the left, ie: smaller masses.
    #     reject_min = mass_dist[(mass_dist < key[0])]
    #     print ("Rejected {} stars with masses in the range [{:.2f}, {:.2f}]"
    #         "\nTotal mass: {:.2f}".format(len(reject_min),
    #             min(reject_min), max(reject_min), sum(reject_min)))
    #     # Masses out of boundary to the right, ie: higher masses.
    #     reject_max = mass_dist[(mass_dist > key[-1])]
    #     print ("Rejected {} stars with masses in the range [{:.2f}, {:.2f}]"
    #         "\nTotal mass: {:.2f}".format(len(reject_max),
    #             min(reject_max), max(reject_max), sum(reject_max)))
    #     print 'Total mass in isochrone: {:0.2f}\n'.format(
    #         sum(mass_dist) - sum(reject_min) - sum(reject_max))
    # except:
    #     pass
    # raw_input()
    ##########################################################################

    # Reject masses in the IMF mass distribution that are located outside of
    # the theoretical isochrone's mass range.
    mass_dist = mass_dist[(mass_dist >= key[0]) & (mass_dist <= key[-1])]

    # Obtain indexes for mass values in the 'mass_dist' array (i.e.: IMF mass
    # distribution) pointing to where these masses should be located within
    # the theoretical isochrone mass distribution.
    closest = find_closest(key, mass_dist)

    # The indexes in 'closest' are used to *replicate* stars in the theoretical
    # isochrone, following the mass distribution given by the IMF.
    # This is: for every mass value in the IMF mass distribution, the star
    # with the closest mass value in the theoretical isochrone is found, and
    # finally these "closest" stars in the isochrone are stored and passed.
    # The "interpolated" isochrone contains as many stars as masses in the
    # IMF distribution were located between the isochrone's mass range.
    isoch_interp = isoch_cut[:, order[closest]]

    return isoch_interp

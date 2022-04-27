
import numpy as np


def main(isoch_cut, mass_ini, mass_dist):
    """
    For each mass in the sampled IMF mass distribution, interpolate its value
    (and those of all the sub-arrays in 'isoch_cut') into the isochrone.

    Masses that fall outside of the isochrone's mass range are rejected.
    """

    # Interpolate sampled masses
    isoch_interp = interp1d(mass_dist, mass_ini, isoch_cut)

    # # This is equivalent to the block above, but slower for more than 5
    # # dimensions, and *very* slightly faster for lower dimensions
    # isoch_interp = np.empty([isoch_cut.shape[0], mass_dist.size])
    # for i, arr in enumerate(isoch_cut):
    #     isoch_interp[i] = np.interp(mass_dist, mass_ini, arr)

    return isoch_interp


def interp1d(x_new, x, y):
    """
    Stripped down version of scipy.interpolate.interp1d. Assumes sorted
    'x' data.

    `x` and `y` are arrays of values used to approximate some function f:
    ``y = f(x)`` using some new values
    y_new = f(x_new)

    x_new.shape --> M
    x.shape     --> N
    y.shape     --> (D, N)
    y_new.T.shape --> (D, M)
    """

    _y = y.T

    # Find where in the original data, the values to interpolate
    # would be inserted.
    # Note: If x_new[n] == x[m], then m is returned by searchsorted.
    x_new_indices = np.searchsorted(x, x_new)

    # Clip x_new_indices so that they are within the range of
    # self.x indices and at least 1. Removes mis-interpolation
    # of x_new[n] = x[0]
    # x_new_indices = x_new_indices.clip(1, len(x) - 1).astype(int)

    # Calculate the slope of regions that each x_new value falls in.
    lo = x_new_indices - 1
    # hi = x_new_indices

    x_lo = x[lo]
    x_hi = x[x_new_indices]
    y_lo = _y[lo]
    y_hi = _y[x_new_indices]

    # Note that the following two expressions rely on the specifics of the
    # broadcasting semantics.
    slope = (y_hi - y_lo) / (x_hi - x_lo)[:, None]

    # Calculate the actual value for each entry in x_new.
    y_new = slope * (x_new - x_lo)[:, None] + y_lo

    return y_new.T


# DEPRECATED FEB 2021

# # Returns the indices that would sort the isochrone with
# # the minimum mass. This is why we use the 'm_ini_idx' index.
# order = isoch_cut[m_ini_idx, :].argsort()

# # Returns an array with the mass values in the theoretical isochrone
# # ordered from min to max.
# key = isoch_cut[m_ini_idx, order]

# ##########################################################################
# # # Uncomment this block to see how many stars are being discarded
# # # for having mass values outside of the isochrone's range.
# # # Also uncomment the first block in the main function that uses the
# # # entire isochrone (instead of one with a max-mag cut).
# # try:
# #     print 'Min, max mass values in isochrone: {:.3f}, {:.3f}'.format(
# #         min(isoch_cut[m_ini_idx]), max(isoch_cut[m_ini_idx]))
# #     print 'Total mass in mass_dist: {:.2f}'.format(sum(mass_dist))
# #     # Masses out of boundary to the left, ie: smaller masses.
# #     reject_min = mass_dist[(mass_dist < key[0])]
# #     print ("Rejected {} stars with masses in the range [{:.2f}, {:.2f}]"
# #         "\nTotal mass: {:.2f}".format(len(reject_min),
# #             min(reject_min), max(reject_min), sum(reject_min)))
# #     # Masses out of boundary to the right, ie: higher masses.
# #     reject_max = mass_dist[(mass_dist > key[-1])]
# #     print ("Rejected {} stars with masses in the range [{:.2f}, {:.2f}]"
# #         "\nTotal mass: {:.2f}".format(len(reject_max),
# #             min(reject_max), max(reject_max), sum(reject_max)))
# #     print 'Total mass in isochrone: {:0.2f}\n'.format(
# #         sum(mass_dist) - sum(reject_min) - sum(reject_max))
# # except:
# #     pass
# # raw_input()
# ##########################################################################

# # Reject masses in the IMF mass distribution that are located outside of
# # the theoretical isochrone's mass range.
# # mass_dist = mass_dist[(mass_dist >= key[0]) & (mass_dist <= key[-1])]
# mass_range_msk = (st_dist_mass[0] >= key[0]) & (st_dist_mass[0] <= key[-1])

# # Obtain indexes for mass values in the 'mass_dist' array (i.e.: IMF mass
# # distribution) pointing to where these masses should be located within
# # the theoretical isochrone mass distribution.
# closest = find_closest(key, st_dist_mass[0])

# # The indexes in 'closest' are used to *replicate* stars in the theoretical
# # isochrone, following the mass distribution given by the IMF.
# # This is: for every mass value in the IMF mass distribution, the star
# # with the closest mass value in the theoretical isochrone is found, and
# # finally these "closest" stars in the isochrone are stored and passed.
# # The "interpolated" isochrone contains as many stars as masses in the
# # IMF distribution were located between the isochrone's mass range.
# isoch_interp = isoch_cut[:, order[closest]]

# def find_closest(key, target):
#     """
#     See: http://stackoverflow.com/a/8929827/1391441
#     Helper function for locating the mass values in the IMF distribution onto
#     the isochrone.
#     General: find closest target element for elements in key.

#     Returns an array of indexes of the same length as 'target'.
#     """

#     # Find indexes where the masses in 'target' should be located in 'key' such
#     # that if the masses in 'target' were inserted *before* these indices, the
#     # order of 'key' would be preserved. I.e.: pair masses in 'target' with the
#     # closest masses in 'key'.
#     # key must be sorted in ascending order.
#     idx = key.searchsorted(target)

#     # Convert indexes in the limits (both left and right) smaller than 1 and
#     # larger than 'len(key) - 1' to 1 and 'len(key) - 1', respectively.
#     idx = np.clip(idx, 1, len(key) - 1)

#     # left < target <= right for each element in 'target'.
#     left = key[idx - 1]
#     right = key[idx]

#     # target - left < right - target is True (or 1) when target is closer to
#     # left and False (or 0) when target is closer to right.
#     # The indexes stored in 'idx' point, for every element (IMF mass) in
#     # 'target' to the closest element (isochrone mass) in 'key'. Thus:
#     # target[XX] <closest> key[idx[XX]]
#     idx -= target - left < right - target

#     return idx

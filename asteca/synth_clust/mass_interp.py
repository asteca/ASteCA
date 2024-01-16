import numpy as np


def main(isoch_cut, m_ini_idx, st_dist_mass):
    """
    For each mass in the sampled IMF mass distribution, interpolate its value
    (and those of all the sub-arrays in 'isoch_cut') into the isochrone.

    Masses that fall outside of the isochrone's mass range have been previously
    rejected.
    """
    # Assumes `mass_ini=isoch_cut[m_ini_idx]` is ordered
    mass_ini = isoch_cut[m_ini_idx]

    # Filter masses in the IMF sampling that are outside of the mass
    # range given by 'isoch_cut' (st_dist_mass[0]: sampled masses from IMF)
    # msk_min = (st_dist_mass[0] >= mass_ini.min())
    # msk_max = (st_dist_mass[0] <= mass_ini.max())
    # (~msk_min).sum(): stars lost below the minimum mass (photometric)
    # (~msk_max).sum(): stars lost above the maximum mass (evolutionary)
    msk_m = (st_dist_mass[0] >= mass_ini.min()) & (st_dist_mass[0] <= mass_ini.max())
    if msk_m.sum() == 0:
        return np.array([])

    mass_dist = st_dist_mass[0][msk_m]
    # # Filter the total mass associated with each sampled mass
    # M_total_arr = st_dist_mass[1][msk_m]

    # Interpolate sampled masses
    isoch_mass = interp1d(mass_ini, mass_dist, isoch_cut)

    # # This is equivalent to the block above, but slower for more than 5
    # # dimensions, and *very* slightly faster for lower dimensions
    # isoch_mass = np.empty([isoch_cut.shape[0], mass_dist.size])
    # for i, arr in enumerate(isoch_cut):
    #     isoch_mass[i] = np.interp(mass_dist, mass_ini, arr)

    return isoch_mass


def interp1d(x, x_new, y):
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

    # _y = y.T

    # # Find where in the original data, the values to interpolate
    # # would be inserted.
    # # Note: If x_new[n] == x[m], then m is returned by searchsorted.
    # x_new_indices = np.searchsorted(x, x_new)

    # # Calculate the slope of regions that each x_new value falls in.
    # lo = x_new_indices - 1
    # # hi = x_new_indices

    # x_lo = x[lo]
    # x_hi = x[x_new_indices]
    # y_lo = _y[lo]
    # y_hi = _y[x_new_indices]

    # # Note that the following two expressions rely on the specifics of the
    # # broadcasting semantics.
    # slope = (y_hi - y_lo) / (x_hi - x_lo)[:, None]

    # # Calculate the actual value for each entry in x_new.
    # y_new = slope * (x_new - x_lo)[:, None] + y_lo

    # return y_new.T

    # This is a slightly optimized version of the function above
    x_new_indices = np.searchsorted(x, x_new)
    lo = x_new_indices - 1

    x_lo = x[lo]
    x_hi = x[x_new_indices]

    y_lo = y[:, lo]
    y_hi = y[:, x_new_indices]

    slope = (y_hi - y_lo) / (x_hi - x_lo)

    x_diff = x_new - x_lo
    y_diff = slope * x_diff

    y_new = y_diff + y_lo

    return y_new

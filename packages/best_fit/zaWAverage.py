
import numpy as np


def main(theor_tracks, fundam_params, varIdxs, model):
    """
    Average a new (weighted) isochrone from the four closest points in the
    (z, a) grid.

    theor_tracks = [m1, m2, .., mN]
    mX = [age1, age2, ..., ageM]
    ageX = [f1,.., c1, c2,.., f1b,.., c1b, c2b,.., bp, mb, m_ini,.., m_bol]
    where:
    fX:  individual filters (mags)
    cX:  colors
    fXb: filters with binary data
    cXb: colors with the binary data
    bp:  binary probabilities
    mb:  binary masses
    m_ini,..., m_bol: six extra parameters.

    WARNING: currently only reading the 'm_ini' extra parameter in
    'isochs_format()'

    """

    model_proper, z_model, a_model, ml, mh, al, ah = properModel(
        fundam_params, model, varIdxs)

    # If (z, a) are both fixed, just return the single processed isochrone
    if ml == al == mh == ah == 0:
        return theor_tracks[ml][al], model_proper

    # The four points in the (z, age) grid that define the box that contains
    # the model value (z_model, a_model)
    z1, z2 = fundam_params[0][ml], fundam_params[0][mh]
    a1, a2 = fundam_params[1][al], fundam_params[1][ah]
    pts = np.array([(z1, a1), (z1, a2), (z2, a1), (z2, a2)])

    # Order: (z1, a1), (z1, a2), (z2, a1), (z2, a2)
    isochs = np.array([
        theor_tracks[ml][al], theor_tracks[ml][ah], theor_tracks[mh][al],
        theor_tracks[mh][ah]])

    # Distances between the (z, a) points in the 'model', and the four
    # points in the (z, a) grid that contain the model point.
    # Fast euclidean distance: https://stackoverflow.com/a/47775357/1391441
    a_min_b = np.array([(z_model, a_model)]) - pts
    dist = np.sqrt(np.einsum('ij,ij->i', a_min_b, a_min_b))

    # If the model has a 0. distance in (z, a) to the closest isochrone,
    # then just return that isochrone.
    try:
        idx = np.where(dist == 0.)[0][0]
        return isochs[idx], model_proper
    except IndexError:
        pass

    # Inverse of the distance.
    inv_d = 1. / dist

    # Weighted average. This way is faster than using 'np.average()'.
    # Scale weights so they add up to 1.
    weights = inv_d / np.sum(inv_d)
    isochrone = isochs[0] * weights[0] + isochs[1] * weights[1] +\
        isochs[2] * weights[2] + isochs[3] * weights[3]

    # The above averaging assumes that the four isochrones are mass-aligned.
    # This is strictly *not true*, but since the four isochrones have
    # relatively similar values, it should be close enough to true.
    # As more points are interpolated into the theoretical isochrones, this
    # will work better.
    #
    # In my tests, attempting an alignment worsens the final averaged
    # isochrone.

    # nn = np.random.randint(0, 100)
    # if nn == 50:
    #     print(model)
    #     print(model_proper)
    #     import matplotlib.pyplot as plt
    #     plt.subplot(131)
    #     plt.scatter(*pts.T, c='r')
    #     # plt.scatter(*pts[isoch_idx][1:].T, c='g')
    #     plt.scatter(model[0], model[1], marker='x', c='g')
    #     # First color
    #     plt.subplot(132)
    #     plt.plot(theor_tracks[ml][al][1], theor_tracks[ml][al][0], c='b')
    #     plt.plot(theor_tracks[ml][ah][1], theor_tracks[ml][ah][0], c='r')
    #     plt.plot(theor_tracks[mh][al][1], theor_tracks[mh][al][0], c='cyan')
    #     plt.plot(theor_tracks[mh][ah][1], theor_tracks[mh][ah][0], c='orange')
    #     plt.plot(isochrone[1], isochrone[0], c='g', ls='--')
    #     plt.gca().invert_yaxis()
    #     # Second color
    #     plt.subplot(133)
    #     plt.plot(theor_tracks[ml][al][-1], theor_tracks[ml][al][0], c='b')
    #     plt.plot(theor_tracks[ml][ah][-1], theor_tracks[ml][ah][0], c='r')
    #     plt.plot(theor_tracks[mh][al][-1], theor_tracks[mh][al][0], c='cyan')
    #     plt.plot(theor_tracks[mh][ah][-1], theor_tracks[mh][ah][0], c='orange')
    #     plt.plot(isochrone[-1], isochrone[0], c='g', ls='--')
    #     plt.gca().invert_yaxis()
    #     plt.show()

    return isochrone, model_proper


def properModel(fundam_params, model, varIdxs):
    """
    Define the 'proper' model with values for (z, a) taken from its grid,
    and filled values for those parameters that are fixed.

    Parameters
    ----------
    model : array
      Array of *free* fundamental parameters only (ie: not in varIdxs).

    Returns
    -------
    model_proper : list
      Stores (E_BV, dm, Mass, b_fr) including the fixed parameters that are
      missing from 'model'.
    z_model, a_model : floats
      The (z, a) values for this model's isochrone.
    ml, mh, al, ah : ints
      Indexes of the (z, a) values in the grid that define the box that enclose
      the (z_model, a_model) values.

    """

    model_proper, j = [], 0
    for i, par in enumerate(fundam_params):
        # Check if this parameter is one of the 'free' parameters.
        if i in varIdxs:
            # If it is the parameter metallicity.
            if i == 0:
                # Select the closest value in the array of allowed values.
                mh = min(len(par) - 1, np.searchsorted(par, model[i - j]))
                ml = mh - 1
                # Define the model's z value
                z_model = model[i - j]
            # If it is the parameter log(age).
            elif i == 1:
                # Select the closest value in the array of allowed values.
                ah = min(len(par) - 1, np.searchsorted(par, model[i - j]))
                al = ah - 1
                a_model = model[i - j]
            else:
                model_proper.append(model[i - j])
        else:
            if i == 0:
                ml = mh = 0
                z_model = fundam_params[0][0]
            elif i == 1:
                al = ah = 0
                a_model = fundam_params[1][0]
            else:
                model_proper.append(par[0])
            j += 1

    return model_proper, z_model, a_model, ml, mh, al, ah

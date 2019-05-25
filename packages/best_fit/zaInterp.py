
import numpy as np
import warnings


def main(theor_tracks, fundam_params, varIdxs, model):
    """
    Interpolate a new isochrone from the four closest points in the (z, a)
    grid.

    The mass value is not interpolated.

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
    This means that '-6' is the index of the initial masses.

    """

    # Select the proper values for the model (ie: that exist in the grid), and
    # choose the (z, a) grid indexes to interpolate the isochrone.
    model_proper, j = [], 0
    for i, par in enumerate(fundam_params):
        # If this parameter is one of the 'free' parameters.
        if i in varIdxs:
            # If it is the parameter metallicity.
            if i == 0:
                # Select the closest value in the array of allowed values.
                mh = min(len(par) - 1, np.searchsorted(par, model[i - j]))
                ml = mh - 1
                model_proper.append(par[mh])
            elif i == 1:
                # Select the closest value in the array of allowed values.
                ah = min(len(par) - 1, np.searchsorted(par, model[i - j]))
                al = ah - 1
                model_proper.append(par[ah])
            elif i == 4:
                # Select the closest value in the array of allowed values for
                # the masses.
                model_proper.append(min(
                    par, key=lambda x: abs(x - model[i - j])))
            else:
                model_proper.append(model[i - j])
        else:
            if i == 0:
                ml = mh = 0
            elif i == 1:
                al = ah = 0
            model_proper.append(par[0])
            j += 1

    # # Minimum and maximum initial mass for each of the four isochrones.
    # mmin = np.min(isochs[:, -6, :], axis=1)
    # mmax = np.max(isochs[:, -6, :], axis=1)

    # Values of the four points in the (z, age) grid that contain the model
    # value (model[0], model[1])
    z1, z2 = fundam_params[0][ml], fundam_params[0][mh]
    a1, a2 = fundam_params[1][al], fundam_params[1][ah]
    pts = np.array([(z1, a1), (z1, a2), (z2, a1), (z2, a2)])

    # Define weights for the average, based on the inverse of the distance
    # to the grid (z, age) points.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # Fast euclidean distance: https://stackoverflow.com/a/47775357/1391441
        a_min_b = np.array([(model[0], model[1])]) - pts
        inv_d = 1. / np.sqrt(np.einsum('ij,ij->i', a_min_b, a_min_b))

        # Order: (z1, a1), (z1, a2), (z2, a1), (z2, a2)
        isochs = np.array([
            theor_tracks[ml][al], theor_tracks[ml][ah], theor_tracks[mh][al],
            theor_tracks[mh][ah]])

        # Sort according to smaller distance to the (model[0], model[1]) point.
        isoch_idx = np.argsort(inv_d)[::-1]

        # Masked weighted average. Source:
        # https://stackoverflow.com/a/35758345/1391441
        x1, x2, x3, x4 = isochs[isoch_idx]
        for x in (x2, x3, x4):
            # Maximum mass difference allowed
            msk = abs(x1[-6] - x[-6]) > .01
            # If the distance in this array is larger than the maximum allowed,
            # mask with the values taken from 'x1'.
            # x[:, msk] = x1[:, msk]
            np.copyto(x, x1, where=msk)

        # # Weighted average with mass "alignment".
        # isochrone = np.average(
        #     np.array([x1, x2, x3, x4]), weights=inv_d[isoch_idx], axis=0)
        # Scale weights so they add up to 1, then add based on them
        weights = inv_d[isoch_idx] / np.sum(inv_d[isoch_idx])
        # If the model has a 0.distance in (z,a) to the closest isochrone, then
        # 'inv_d' will be 'inf', and 'weights' will be 'nan' (all other values
        # will be zero) Replace that weight with 1.
        weights[np.isnan(weights)] = 1.
        isochrone = x1 * weights[0] + x2 * weights[1] + x3 * weights[2] +\
            x4 * weights[3]

    # This way is *marginally* faster
    # wgts = D * inv_d[isoch_idx]
    # x1, x2, x3, x4 = isochs[isoch_idx]
    # for x in (x2, x3, x4):
    #     # Maximum mass difference allowed
    #     msk = abs(x1[-6] - x[-6]) > .01
    #     x[:, msk] = x1[:, msk]
    # isochrone = np.sum(np.array([
    #     x1 * wgts[0], x2 * wgts[1], x3 * wgts[2], x4 * wgts[3]]), 0)

    # Weighted version, no mass alignment.
    # isochrone = np.average(np.array([
    #     theor_tracks[ml][al], theor_tracks[ml][ah], theor_tracks[mh][al],
    #     theor_tracks[mh][ah]]), weights=D * inv_d, axis=0)

    # Simpler mean version, no weights.
    # isochrone = np.mean([
    #     theor_tracks[ml][al], theor_tracks[ml][ah], theor_tracks[mh][al],
    #     theor_tracks[mh][ah]], axis=0)

    # nn = np.random.randint(0, 100)
    # if nn == 50:
    #     print(model)
    #     print(model_proper)
    #     import matplotlib.pyplot as plt
    #     plt.subplot(131)
    #     plt.scatter(*pts[isoch_idx][0], c='r')
    #     plt.scatter(*pts[isoch_idx][1:].T, c='g')
    #     plt.scatter(model[0], model[1], marker='x')
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
    #     plt.plot(theor_tracks[ml][al][2], theor_tracks[ml][al][0], c='b')
    #     plt.plot(theor_tracks[ml][ah][2], theor_tracks[ml][ah][0], c='r')
    #     plt.plot(theor_tracks[mh][al][2], theor_tracks[mh][al][0], c='cyan')
    #     plt.plot(theor_tracks[mh][ah][2], theor_tracks[mh][ah][0], c='orange')
    #     plt.plot(isochrone[2], isochrone[0], c='g', ls='--')
    #     plt.gca().invert_yaxis()
    #     plt.show()

    return isochrone, model_proper

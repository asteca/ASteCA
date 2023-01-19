
import numpy as np


def main(
    theor_tracks, fundam_params, m_ini_idx, z_model, a_model,
        ml, mh, al, ah):
    """
    Generate a new "weighted" isochrone from the four closest points in the
    (z, a) grid.

    theor_tracks = [m1, m2, .., mN]
    mX = [age1, age2, ..., ageM]
    ageX = [f1,.., c1, c2,.., M_ini, f1b,.., c1b, c2b,.., M_b]
    where:
    fX:  individual filters (mags)
    cX:  colors
    M_ini: initial mass
    fXb: filters with binary data
    cXb: colors with the binary data
    M_b:  binary masses

    It is important that the returned isochrone is a *copy* of the
    theor_tracks[mx][ax] array if no averaging is done. Otherwise the
    'theor_tracks' are modified.
    """

    # If (z, a) are both fixed, just return the single processed isochrone
    if ml == al == mh == ah == 0:
        # The np.array() is important to avoid overwriting 'theor_tracks'
        return np.array(theor_tracks[ml][al])

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
    # Don't take the square root, it's not necessary
    dist = np.einsum('ij,ij->i', a_min_b, a_min_b)

    # If the model has a 0. distance in (z, a) to the closest isochrone,
    # then just return that isochrone.
    try:
        idx = np.where(dist == 0.)[0][0]
        return isochs[idx]
    except IndexError:
        pass

    # Weighted average by the (inverse) distance to the four (z, a) grid
    # points. This way is faster than using 'np.average()'.
    # Inverse of the distance.
    inv_d = 1. / dist
    weights = inv_d / sum(inv_d)
    isochrone = isochs[0] * weights[0] + isochs[1] * weights[1] +\
        isochs[2] * weights[2] + isochs[3] * weights[3]

    # DO NOT average the masses or their distribution will be lost. We use the
    # closest isochrone values.
    idx = np.argmin(dist)
    isochrone[m_ini_idx] = isochs[idx][m_ini_idx]
    # Now for the secondary masses
    isochrone[-1] = isochs[idx][-1]

    # isochrone = theor_tracks[ml][al] * weights[0] +\
    #     theor_tracks[ml][ah] * weights[1] +\
    #     theor_tracks[mh][al] * weights[2] +\
    #     theor_tracks[mh][ah] * weights[3]

    # #
    # isochrone = np.average(isochs, weights=weights, axis=0)

    # idxs = (weights * mass_idx).astype(int)
    # # Order: (z1, a1), (z1, a2), (z2, a1), (z2, a2)
    # isochrone = np.concatenate([
    #     theor_tracks[ml][al][:, :idxs[0]],
    #     theor_tracks[ml][ah][:, idxs[0]:idxs[0] + idxs[1]],
    #     theor_tracks[mh][al][:, idxs[0] + idxs[1]:idxs[0] + idxs[1] + idxs[2]],
    #     theor_tracks[mh][ah][:, idxs[0] + idxs[1] + idxs[2]:mass_idx]],
    #     axis=1)

    # tt = np.sum([t1, t2, t3, t4, t5, t6, t7])
    # print(t1/tt, t2/tt, t3/tt, t4/tt, t5/tt, t6/tt, t7/tt)

    # # if np.random.randint(0, 100) == 10:
    # import matplotlib.pyplot as plt
    # plt.subplot(121)
    # plt.scatter(*pts.T, c='r')
    # plt.scatter(z_model, a_model, marker='x', c='g')
    # # First color
    # plt.subplot(122)
    # plt.scatter(
    #     theor_tracks[ml][al][1], theor_tracks[ml][al][0], c='b', alpha=.5)
    # plt.scatter(
    #     theor_tracks[ml][ah][1], theor_tracks[ml][ah][0], c='r', alpha=.5)
    # plt.scatter(
    #     theor_tracks[mh][al][1], theor_tracks[mh][al][0], c='cyan', alpha=.5)
    # plt.scatter(
    #     theor_tracks[mh][ah][1], theor_tracks[mh][ah][0], c='orange', alpha=.5)
    # plt.scatter(isochrone[1], isochrone[0], c='g')
    # plt.gca().invert_yaxis()
    # # Second color
    # # plt.subplot(133)
    # # plt.scatter(theor_tracks[ml][al][2], theor_tracks[ml][al][0], c='b')
    # # plt.scatter(theor_tracks[ml][ah][2], theor_tracks[ml][ah][0], c='r')
    # # plt.scatter(theor_tracks[mh][al][2], theor_tracks[mh][al][0], c='cyan')
    # # plt.scatter(theor_tracks[mh][ah][2], theor_tracks[mh][ah][0], c='orange')
    # # plt.scatter(isochrone[2], isochrone[0], c='g', ls='--')
    # # plt.gca().invert_yaxis()
    # plt.show()
    # print("zaWavrg")
    # breakpoint()

    return isochrone

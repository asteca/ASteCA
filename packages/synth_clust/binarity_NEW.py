
import numpy as np
# from .mass_interp import find_closest
from .. import update_progress


def main(isoch_mass, bin_frac, m_ini_idx, N_fc):
    """
    Update the randomly selected fraction of binary stars.
    """

    # If the binary fraction is zero, skip the whole process.
    if bin_frac > 0.:

        # Select a fraction of stars to be binaries, according to the random
        # probabilities assigned before.
        # bin_indxs = isoch_mass[m_ini_idx - 2] <= bin_frac
        bin_indxs = randomPairs(isoch_mass.shape[1])
        import pdb; pdb.set_trace()  # breakpoint e7bdec85 //


        # Index of the first binary magnitude, stored in the theoretical
        # isochrones list.
        mag_ini = N_fc[0] + N_fc[1]

        # Update array with new values of magnitudes, colors, and masses.
        # Binary magnitudes.
        for i in range(N_fc[0]):
            isoch_mass[i][bin_indxs] = mag_combine(
                isoch_mass[i], isoch_mass[i][bin_indxs])
            # isoch_mass[mag_ini + i][bin_indxs]
        # Binary colors.
        for i in range(N_fc[1]):
            isoch_mass[N_fc[0] + i][bin_indxs] =\
                isoch_mass[mag_ini + N_fc[0] + i][bin_indxs]
        # Binary masses.
        isoch_mass[m_ini_idx][bin_indxs] = isoch_mass[m_ini_idx - 1][bin_indxs]
        # IDs of binary systems start with a '2.'.
        isoch_mass[m_ini_idx - 2][bin_indxs] = 2.

    return isoch_mass


def binarGen(mass_interp_isochs, N_fc):
    """
    Called by isoch_params().

    0. Assign N unique indexes by dividing the range of stars by their total
       number.

    For each theoretical isochrone defined.

        1. Draw random secondary masses for *all* stars.
        2. Find stars in the isochrone with the closest mass.
        3. Calculate the magnitudes for the binary system.
        4. Calculate the colors for the binary system.
        5. Calculate the masses for the binary system.
        7. Assign unique probability values for each interpolated star.

    """
    # NOT USED ANYMORE
    # bin_mass_ratio

    # met_probs, age_probs = metageProbs()
    # unq_probs = randVals(N_mass_interp)

    fc_0 = N_fc[0] + N_fc[1]
    fc_1 = int(fc_0 + 2 * N_fc[1])

    mags_theor, extra_pars = mass_interp_isochs[:, :, 0:N_fc[0], :],\
        mass_interp_isochs[:, :, -1:, :]
    mags_cols_theor = mass_interp_isochs[:, :, fc_0:fc_1, :]

    mags_binar, cols_binar, mass_binar = [], [], []
    # For each metallicity defined.
    N_mags_theor = len(mags_theor)
    for mx, _ in enumerate(mags_theor):
        # Find closest metallicity
        # iz = np.searchsorted(met_probs, all_met_vals[mx])

        mag_bin, col_bin, mass_bin = [], [], []
        # For each age defined.
        for ax, isoch in enumerate(_):
            # Find closest age
            # ia = np.searchsorted(age_probs, all_age_vals[ax])

            # # This ensures that the same (z, a) pair points to the same
            # # 'unq_probs' values (for the same random seed), no matter
            # # the ranges used for these parameters.
            # idx = (iz + ia) % len(unq_probs)

            # Extract initial masses for this isochrone.
            mass_ini = extra_pars[mx][ax][0]

            # # Calculate random secondary masses of these binary stars
            # # between bin_mass_ratio*m1 and m1, where m1 is the primary
            # # mass.
            # # m2 = np.random.uniform(bin_mass_ratio * mass_ini, mass_ini)
            # m2 = fracs * mass_ini

            # # If any secondary mass falls outside of the lower isochrone's
            # # mass range, change its value to the min value.
            # # This is faster than np.clip()
            # m2 = np.maximum(np.min(mass_ini), m2)

            # # Obtain indexes for mass values in the 'm2' array pointing to
            # # the closest mass in the theoretical isochrone.
            # bin_m_close = find_closest(mass_ini, m2)

            bin_idxs = randomPairs(len(mass_ini))

            # Calculate unresolved binary magnitude for each
            # filter/magnitude defined.
            temp = []
            for mag in isoch:
                temp.append(mag_combine(mag, mag[bin_idxs]))
            mag_bin.append(temp)

            # Calculate unresolved color for each color defined.
            filt_colors = mags_cols_theor[mx][ax]

            # Filters composing the colors, i.e.: C = (f1 - f2).
            f1, f2 = [], []
            # The [::2] slice indicates even positions, starting from 0.
            for m in filt_colors[::2]:
                f1.append(mag_combine(m, m[bin_idxs]))
            # The [1::2] slice indicates odd positions.
            for m in filt_colors[1::2]:
                f2.append(mag_combine(m, m[bin_idxs]))

            # Create the colors affected by binarity.
            temp = []
            for f_1, f_2 in zip(*[f1, f2]):
                temp.append(f_1 - f_2)
            col_bin.append(temp)

            # Add masses to obtain the binary system's mass.
            mass_bin.append([mass_ini + mass_ini[bin_idxs]])

            # # Binary system probability
            # probs = np.zeros(len(mass_ini))
            # prob_bin.append([probs])

            # import matplotlib.pyplot as plt
            # print(mx, ax)
            # print(min(isoch[0]), max(isoch[0]))
            # plt.scatter(filt_colors[0] - filt_colors[1], isoch[0], c='g')
            # plt.scatter(col_bin[0], mag_bin[0], c='r')
            # plt.gca().invert_yaxis()
            # plt.show()

        # Store for each metallicity value.
        mags_binar.append(mag_bin)
        cols_binar.append(col_bin)
        # probs_binar.append(prob_bin)
        mass_binar.append(mass_bin)

        update_progress.updt(N_mags_theor, mx + 1)

    return np.array(mags_binar), np.array(cols_binar), np.array(mass_binar)


# def metageProbs(
#     zmin=0., zmax=0.06, amin=6., amax=10.5, N_mets=50000,
#         N_ages=50000):
#     """
#     Process the required random values making sure that they are reproducible
#     to the maximum possible extent. In the N_mets-->inf, N_ages-->inf limit
#     all (z, a) pairs have a unique index assigned and hence always point to
#     the exact same position in 'unq_probs'.

#     HARDCODED
#     zmin, zmax, amin, amax : full range for each parameter
#     N_mets, N_ages : number of elements in each array
#     """
#     met_probs = np.linspace(zmin, zmax, N_mets)
#     age_probs = np.linspace(amin, amax, N_ages)

#     return met_probs, age_probs


# def randVals(N_mass_interp, N_unq_probs=10000):
#     """
#     The 'N_unq_probs' value is not that important. In the limit
#     N_unq_probs-->inf every (z, a) pair has a unique array of probabilities
#     assigned.

#     N_unq_probs : number of elements in array
#     """
#     Nh = N_mass_interp / 2.

#     # All theoretical isochrones are interpolated with the same length,
#     # assign unique binarity probabilities to each star randomly.
#     b_probs = np.arange(Nh) / float(Nh)

#     # As long as the length of 'b_probs' stays the same (i.e.:
#     # as long as N_mass_interp stays the same), this will produce
#     # the same results for the same random seed.
#     unq_probs = []
#     for _ in range(N_unq_probs):
#         # Shuffle binarity probabilities.
#         np.random.shuffle(b_probs)
#         unq_probs.append(list(b_probs))

#     return unq_probs


def randomPairs(N):
    """
    https://stackoverflow.com/q/65361557/1391441
    """

    # # Faster?
    # result = [None] * N
    # A  = np.random.choice(range(N), N, replace=False)
    # for a, b in zip(A[::2], A[1::2]):
    #     result[a] = b
    #     result[b] = a
    # return np.array(result)

    result = [None] * N
    if N & 1 == 0:
        # A = random.sample(range(N), N // 2)
        A = np.random.choice(range(N), N // 2, replace=False)
        # B = random.sample(list(set(range(N)) - set(A)), N // 2)
        B = np.random.choice(
            list(set(range(N)) - set(A)), N // 2, replace=False)
        for a, b in zip(A, B):
            result[a] = b
            result[b] = a
    return np.array(result)


def mag_combine(m1, m2):
    """
    Combine two magnitudes. This is a faster re-ordering of the standard
    formula:

    -2.5 * np.log10(10 ** (-0.4 * m1) + 10 ** (-0.4 * m2))

    """

    c = 10 ** -.4
    mbin = -2.5 * (-.4 * m1 + np.log10(1. + c ** (m2 - m1)))

    return mbin


import numpy as np
from .. import update_progress


def main(isoch_mass, bin_frac, m_ini_idx, N_fc):
    """
    Update the randomly selected fraction of binary stars.

    The returned array is of shape:

    isoch_mass.shape = [mag, col1, (col2), m_ini, binar_idxs]

    or simply

    isoch_mass.shape = [mag, col1, (col2), m_ini]

    if bin_frac = 0.
    """

    # If the binary fraction is zero, skip the whole process.
    if bin_frac > 0.:

        # Select a fraction of stars to be binaries, according to the random
        # probabilities assigned before.
        bin_indxs = isoch_mass[-1] <= bin_frac

        # Index of the first binary magnitude.
        # mag_binar = m_ini_idx + 1
        # Update array with new values of magnitudes and colors.
        for i in range(N_fc[0] + N_fc[1]):
            isoch_mass[i][bin_indxs] = isoch_mass[m_ini_idx + 1 + i][bin_indxs]

        # Update masses of single systems with binary masses (index '-2').
        isoch_mass[m_ini_idx][bin_indxs] = isoch_mass[-2][bin_indxs]

        # Update the (now useless) binary magnitude array so that the
        # binary systems are identified with a '-99.' value.
        isoch_mass[m_ini_idx + 1][bin_indxs] = -99.

        # Discard data not used anymore. This changes the shape of the returned
        # array
        isoch_mass = isoch_mass[:m_ini_idx + 2, :]

    return isoch_mass


def binarGen(
    data_tracks, m_ini_idx, interp_tracks, mags_cols_intp, all_met_vals,
    all_age_vals, bin_mass_ratio, N_fc, zmin=0., zmax=0.1, amin=6.,
        amax=13, N_met_age=100000, N_bprobs=500000):
    """
    0. Assign N unique indexes by dividing the range of stars by their total
       number.

    For each theoretical isochrone defined.

        1. Draw random secondary masses for *all* stars.
        2. Find stars in the isochrone with the closest mass.
        3. Calculate the magnitudes for the binary system.
        4. Calculate the colors for the binary system.
        5. Calculate the masses for the binary system.
        7. Assign unique probability values for each interpolated star.

    HARDCODED
    zmin, zmax, amin, amax : full range for each parameter
    N_met_age : controls the number of indexes generated for the (z, a) values.
      It should be large enough so that all possible (z, a) values from the
      grid (metallicity files) are assigned a unique index.
    N_bprobs  : number of random elements [0., 1.] to generate in 'b_probs'

    """
    from .set_rand_seed import np

    # If binary_fraction=0 don't bother obtaining the binary magnitudes,
    # colors, etc.
    print("Generating binary data (b_mr=[{:.2f}, 1.])".format(
        bin_mass_ratio))

    # Non mass-sampled interpolated tracks, i.e tracks as read from the
    # input files
    # Shape: (Nd, Nz, Na, Np)
    mags_track, _, extra_pars_track, mags_cols_track = data_tracks

    # interp_tracks.shape = (Nz, Na, Nd, Np)
    tot_tracks = interp_tracks.shape[0]
    N_mass_interp = interp_tracks.shape[-1]
    extra_pars = interp_tracks[:, :, m_ini_idx:m_ini_idx + 1, :]

    # Fractions for second mass
    fracs = np.random.uniform(bin_mass_ratio, 1., N_mass_interp)

    # unq_probs = randVals(N_mass_interp)

    # In the N_zap-->inf limit all (z, a) pairs have a unique index
    # assigned and hence always point to the exact same position in
    # 'b_probs'.
    b_probs = np.random.uniform(0., 1., 500000)
    met_probs = np.linspace(zmin, zmax, N_met_age)
    age_probs = np.linspace(amin, amax, N_met_age)

    # For each metallicity defined.
    for mx, _ in enumerate(interp_tracks):
        # For each age defined.
        for ax, isoch in enumerate(_):

            # Extract initial masses for this isochrone. Assumes that
            # the initial masses are in the '0' index.
            mass_ini = extra_pars[mx][ax][0]
            mass_ini_tracks = extra_pars_track[0][mx][ax]

            # Calculate random secondary masses of these binary stars
            # between bin_mass_ratio*m1 and m1, where m1 is the primary
            # mass.
            m2 = fracs * mass_ini

            # Calculate unresolved binary magnitude for each
            # filter/magnitude defined.
            for i, mag in enumerate(isoch[:N_fc[0]]):
                mag_m2 = np.interp(
                    m2, mass_ini_tracks, mags_track[i][mx][ax])
                interp_tracks[mx][ax][m_ini_idx + 1 + i] = mag_combine(
                    mag, mag_m2)

            # Calculate unresolved color for each color defined.
            # The [::2] slice indicates even positions, starting from 0.
            even_f_cols = mags_cols_intp[mx][ax][::2]
            even_f_cols_tracks = mags_cols_track[::2]
            # The [1::2] slice indicates odd positions.
            odd_f_cols = mags_cols_intp[mx][ax][1::2]
            odd_f_cols_tracks = mags_cols_track[1::2]

            # Filters composing the colors, i.e.: C = (f1 - f2).
            f1, f2 = [], []
            for i, m in enumerate(even_f_cols):
                f_m2 = np.interp(
                    m2, mass_ini_tracks, even_f_cols_tracks[i][mx][ax])
                f1.append(mag_combine(m, f_m2))
            for i, m in enumerate(odd_f_cols):
                f2_m2 = np.interp(
                    m2, mass_ini_tracks, odd_f_cols_tracks[i][mx][ax])
                f2.append(mag_combine(m, f2_m2))

            # Create the colors affected by binarity.
            for i, (f_1, f_2) in enumerate(zip(*[f1, f2])):
                interp_tracks[mx][ax][m_ini_idx + 1 + N_fc[0] + i] =\
                    f_1 - f_2

            # import matplotlib.pyplot as plt
            # print(mx, ax)
            # print(min(isoch[0]), max(isoch[0]))
            # # plt.subplot(121)
            # plt.title("Gaia eDR3")
            # plt.scatter(mags_cols_intp[mx][ax][0] - mags_cols_intp[mx][ax][1], isoch[0], c='g')
            # # mag_bin0 = interp_tracks[mx][ax][m_ini_idx + 1]
            # # col_bin0 = interp_tracks[mx][ax][m_ini_idx + 1 + N_fc[0] + 0]
            # # plt.scatter(col_bin0, mag_bin0, c='r', alpha=.5)
            # plt.gca().invert_yaxis()
            # # # Second color
            # # plt.subplot(122)
            # # col_bin1 = interp_tracks[mx][ax][m_ini_idx + 1 + N_fc[0] + 1]
            # # plt.scatter(mags_cols_intp[mx][ax][2] - mags_cols_intp[mx][ax][3], isoch[0], c='g')
            # # plt.scatter(col_bin1, mag_bin0, c='r', alpha=.5)
            # # plt.gca().invert_yaxis()
            # plt.show()

            # Add masses to obtain the binary system's mass.
            interp_tracks[mx][ax][m_ini_idx + 1 + N_fc[0] + N_fc[1]] =\
                mass_ini + m2

            # Find closest met & age values
            iz = np.searchsorted(met_probs, all_met_vals[mx])
            ia = np.searchsorted(age_probs, all_age_vals[ax])
            # This ensures that the same (z, a) pair points to the same
            # 'unq_probs' values (for the same random seed), no matter
            # the ranges used for these parameters. Will only work if
            # 'N_unq_probs' is large enough.
            # idx = (iz + ia) % len(unq_probs)
            # interp_tracks[mx][ax][m_ini_idx + 1 + N_fc[0] + N_fc[1] + 1] =\
            #     unq_probs[idx]

            # idx0 = iz + ia
            # idx1 = idx0 + N_mass_interp
            # if idx1 >= N_bprobs:

            # random.Random(iz).shuffle(b_probs)
            # random.Random(ia).shuffle(b_probs)

            # print(all_met_vals[mx], all_age_vals[ax], iz, ia, b_probs[iz + ia:iz + ia + N_mass_interp][:3])

            interp_tracks[mx][ax][m_ini_idx + 1 + N_fc[0] + N_fc[1] + 1] =\
                b_probs[iz + ia:iz + ia + N_mass_interp]
            # interp_tracks[mx][ax][m_ini_idx + 1 + N_fc[0] + N_fc[1] + 1] =\
            #     np.random.uniform(0., 1., N_mass_interp)

        update_progress.updt(tot_tracks, mx + 1)

    return interp_tracks


def randVals(N_mass_interp, N_unq_probs=1000):
    """
    Process the required random values making sure that they are reproducible
    to the maximum possible extent.

    In the limit N_unq_probs-->inf every (z, a) pair has a unique array of
    probabilities assigned.

    HARDCODED
    N_unq_probs : number of XXXX
    """

    # All theoretical isochrones are interpolated with the same length,
    # assign unique binarity probabilities to each star randomly.
    b_probs = np.arange(N_mass_interp) / float(N_mass_interp)

    # As long as the length of 'b_probs' stays the same (i.e.:
    # as long as N_mass_interp stays the same), this will produce
    # the same results for the same random seed.

    unq_probs = []
    for _ in range(N_unq_probs):
        # Shuffle binarity probabilities.
        np.random.shuffle(b_probs)
        unq_probs.append(list(b_probs))

    return unq_probs


def mag_combine(m1, m2):
    """
    Combine two magnitudes. This is a faster re-ordering of the standard
    formula:

    -2.5 * np.log10(10 ** (-0.4 * m1) + 10 ** (-0.4 * m2))

    """

    c = 10 ** -.4
    mbin = -2.5 * (-.4 * m1 + np.log10(1. + c ** (m2 - m1)))

    return mbin

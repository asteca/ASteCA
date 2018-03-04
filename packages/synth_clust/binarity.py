
import numpy as np
from mass_interp import find_closest


def mag_combine(m1, m2):
    """
    Combine two magnitudes. This is a faster re-ordering of the standard
    formula:

    -2.5 * np.log10(10 ** (-0.4 * m1) + 10 ** (-0.4 * m2))

    """
    c = 10 ** -.4
    # This catches an overflow warning issued because some Marigo isochrones
    # contain huge values in the U filter. Not sure if this happens with
    # other systems/filters. See issue #375
    np.warnings.filterwarnings('ignore')
    mbin = -2.5 * (-.4 * m1 + np.log10(1. + c ** (m2 - m1)))

    return mbin


def binarGen(
        binar_fracs, N_interp, mags_theor, mags_cols_theor, extra_pars,
        bin_mass_ratio):
    '''

    0. Assign N unique indexes by dividing range of stars by their total
       number.

    For each theoretical isochrone defined.

        1. Draw random secondary masses for *all* stars.
        2. Find stars in the isochrone with the closest mass.
        3. Calculate the magnitudes for the binary system.
        4. Calculate the colors for the binary system.
        5. Calculate the masses for the binary system.
        7. Assign unique probability values for each interpolated star.

    '''

    # If binary_fraction = 0. don't bother obtaining the binary magnitudes,
    # colors, etc.
    if binar_fracs.any():

        print("Generating binary data (b_mr={:.2f})\n".format(bin_mass_ratio))
        # All theoretical isochrones are interpolated with the same length,
        # assign unique binarity probabilities to each star randomly.
        unq_b_probs = np.arange(N_interp) / float(N_interp)

        mags_binar, cols_binar, probs_binar, mass_binar = [], [], [], []
        # For each metallicity defined.
        for mx, _ in enumerate(mags_theor):

            mag_bin, col_bin, prob_bin, mass_bin = [], [], [], []
            # For each age defined.
            for ax, isoch in enumerate(_):

                # Extract initial masses for this isochrone.
                mass_ini = extra_pars[mx][ax][0]

                # Calculate random secondary masses of these binary stars
                # between bin_mass_ratio*m1 and m1, where m1 is the primary
                # mass.
                m2 = np.random.uniform(bin_mass_ratio * mass_ini, mass_ini)
                # If any secondary mass falls outside of the lower isochrone's
                # mass range, change its value to the min value.
                m2 = np.maximum(np.min(mass_ini), m2)

                # Obtain indexes for mass values in the 'm2' array pointing to
                # the closest mass in the theoretical isochrone.
                bin_m_close = find_closest(mass_ini, m2)

                # Calculate unresolved binary magnitude for each
                # filter/magnitude defined.
                temp = []
                for mag in isoch:
                    temp.append(mag_combine(mag, mag[bin_m_close]))
                mag_bin.append(temp)

                # Calculate unresolved color for each color defined.
                filt_colors = mags_cols_theor[mx][ax]

                # Filters composing the colors, i.e.: C = (f1 - f2).
                f1, f2 = [], []
                # The [::2] slice indicates even positions, starting from 0.
                for m in filt_colors[::2]:
                    f1.append(mag_combine(m, m[bin_m_close]))
                # The [1::2] slice indicates odd positions.
                for m in filt_colors[1::2]:
                    f2.append(mag_combine(m, m[bin_m_close]))

                # Create the colors affected by binarity.
                temp = []
                for f_1, f_2 in zip(*[f1, f2]):
                    temp.append(f_1 - f_2)
                col_bin.append(temp)

                # import matplotlib.pyplot as plt
                # print(mx, ax)
                # print(min(isoch[0]), max(isoch[0]))
                # plt.scatter(filt_colors[0] - filt_colors[1], isoch[0], c='g')
                # plt.scatter(col_bin[0], mag_bin[0], c='r')
                # plt.gca().invert_yaxis()
                # plt.show()

                # Add masses to obtain the binary system's mass.
                mass_bin.append([mass_ini + mass_ini[bin_m_close]])

                # Shuffle binarity probabilities.
                np.random.shuffle(unq_b_probs)
                prob_bin.append([unq_b_probs])

            # Store for each metallicity value.
            mags_binar.append(mag_bin)
            cols_binar.append(col_bin)
            probs_binar.append(prob_bin)
            mass_binar.append(mass_bin)
    else:
        mags_binar, cols_binar, mass_binar =\
            np.zeros(np.shape(mags_theor)),\
            np.zeros(np.shape(mags_cols_theor)),\
            np.array(extra_pars)[:, :, :1, :]
        probs_binar = np.zeros(np.shape(mass_binar))

    return mags_binar, cols_binar, probs_binar, mass_binar


def main(isoch_mass, bin_frac, m_ini, N_fc):
    '''
    Update the randomly selected fraction of binary stars.
    '''

    # If the binary fraction is zero, skip the whole process.
    if bin_frac > 0.:

        # Select a fraction of stars to be binaries, according to the random
        # probabilities assigned before.
        bin_indxs = isoch_mass[m_ini - 2] <= bin_frac

        # Index of the first binary magnitude, stored in the theoretical
        # isochrones list. The 5 is for the remaining 5 extra parameters.
        mag_ini = N_fc[0] + N_fc[1]

        # Update array with new values of magnitudes, colors, and masses.
        # New magnitudes.
        for i in range(N_fc[0]):
            isoch_mass[i][bin_indxs] = isoch_mass[mag_ini + i][bin_indxs]
        # New colors.
        for i in range(N_fc[1]):
            isoch_mass[N_fc[0] + i][bin_indxs] =\
                isoch_mass[mag_ini + N_fc[0] + i][bin_indxs]

        # Update masses.
        isoch_mass[m_ini][bin_indxs] = isoch_mass[m_ini - 1][bin_indxs]

        # Update binary systems to a '2.'.
        isoch_mass[m_ini - 2][bin_indxs] = 2.

    return isoch_mass


import numpy as np
import random
import mass_interp


def mag_combine(m1, m2):
    """
    Combine two magnitudes.
    """
    return -2.5 * np.log10(10 ** (-0.4 * m1) + 10 ** (-0.4 * m2))


def main(isoch_mass, isoch_cut, bin_frac, bin_mass_ratio, N_fc):
    '''
    Randomly select a fraction of stars to be binaries.
    '''
    # If the binary fraction is zero, skip the whole process.
    if bin_frac > 0.:
        # Indexes of the randomly selected stars in isoch_mass.
        bin_indxs = random.sample(range(len(isoch_mass[0])),
                                  int(bin_frac * len(isoch_mass[0])))

        # Calculate the secondary masses of these binary stars between
        # bin_mass_ratio*m1 and m1, where m1 is the primary mass.
        # Index of m_ini, stored in the theoretical isochrones.
        m_ini = N_fc[0] + N_fc[1] + 2 * N_fc[1]
        # Primary masses.
        m1 = np.asarray(isoch_mass[m_ini][bin_indxs])
        # Secondary masses.
        mass_bin0 = np.random.uniform(bin_mass_ratio * m1, m1)
        # If any secondary mass falls outside of the lower isochrone's mass
        # range, change its value to the min value.
        mass_bin = np.maximum(min(isoch_mass[m_ini]), mass_bin0)

        # Find color and magnitude values for each secondary star. This will
        # slightly change the values of the masses, since they will be
        # assigned to the closest value found in the interpolated isochrone.
        bin_isoch = mass_interp.main(isoch_cut, mass_bin, N_fc)

        # Calculate unresolved binary magnitude for each filter/magnitude
        # defined.
        mag_bin = []
        for i, m in enumerate(isoch_mass[:N_fc[0]]):
            mag_bin.append(mag_combine(m[bin_indxs], bin_isoch[i]))

        # Calculate unresolved color for each color defined.
        # Lower/upper index where filters that make up colors start/finish.
        fcl = N_fc[0] + N_fc[1]
        fcu = fcl + 2 * N_fc[1]
        # Extract only the portion that contains the color filters.
        filt_colors, bin_f_cols = isoch_mass[fcl:fcu], bin_isoch[fcl:fcu]

        # Filters composing the colors, i.e.: C = (f1 - f2).
        f1, f2 = [], []
        # The [::2] slice indicates even positions, starting from 0.
        for m, bin_m in zip(*[filt_colors[::2], bin_f_cols[::2]]):
            f1.append(mag_combine(m[bin_indxs], bin_m))
        # The [1::2] slice indicates odd positions.
        for m, bin_m in zip(*[filt_colors[1::2], bin_f_cols[1::2]]):
            f2.append(mag_combine(m[bin_indxs], bin_m))

        # Create the colors affected by binarity.
        col_bin = []
        for f_1, f_2 in zip(*[f1, f2]):
            col_bin.append(f_1 - f_2)

        # Add masses to obtain the binary system's mass.
        mass_bin = isoch_mass[m_ini][bin_indxs] + bin_isoch[m_ini]

        # Update array with new values of magnitudes, colors, and masses.
        for i in range(N_fc[0]):
            for indx, j in enumerate(bin_indxs):
                isoch_mass[i][j] = mag_bin[i][indx]

        for i in range(N_fc[1]):
            for indx, j in enumerate(bin_indxs):
                isoch_mass[N_fc[0] + i][j] = col_bin[i][indx]

        for indx, j in enumerate(bin_indxs):
            isoch_mass[m_ini][j] = mass_bin[indx]

    return isoch_mass

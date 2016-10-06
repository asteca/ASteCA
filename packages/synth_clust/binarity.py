
import numpy as np
import random


def main(isoch_mass, isoch_cut, bin_frac, bin_mass_ratio):
    '''
    Randomly select a fraction of stars to be binaries.
    '''
    # Indexes of the randomly selected stars in isoch_mass.
    bin_indxs = random.sample(range(len(isoch_mass[0])),
                              int(bin_frac * len(isoch_mass[0])))

    # Calculate the secondary masses of these binary stars between
    # bin_mass_ratio*m1 and m1, where m1 is the primary mass.
    # Primary masses.
    m1 = np.asarray(isoch_mass[2][bin_indxs])
    # Secondary masses.
    mass_bin0 = np.random.uniform(bin_mass_ratio * m1, m1)
    # This prevents a rare error where apparently mass_bin0 is a float.
    if type(mass_bin0) is not float:

        # If any secondary mass falls outside of the lower isochrone's mass
        # range, change its value to the min value.
        mass_bin = np.maximum(min(isoch_mass[2]), mass_bin0)

        # Find color and magnitude values for each secondary star. This will
        # slightly change the values of the masses since they will be
        # assigned to the closest value found in the interpolated isochrone.
        bin_isoch = mass_interp(isoch_cut, mass_bin)

        # Obtain color, magnitude and masses for each binary system.
        # Transform color to the second magnitude before obtaining
        # the new binary magnitude.
        if cmd_sel in {2, 5, 9}:
            # E.g.: V vs (V-I)
            mag2_iso = isoch_mass[1][bin_indxs] - isoch_mass[0][bin_indxs]
            mag2_bin = bin_isoch[1] - bin_isoch[0]
        else:
            # E.g.: V vs (B-V)
            mag2_iso = isoch_mass[0][bin_indxs] + isoch_mass[1][bin_indxs]
            mag2_bin = bin_isoch[0] + bin_isoch[1]
        col_mag_bin = -2.5 * np.log10(10 ** (-0.4 * mag2_iso) +
                                      10 ** (-0.4 * mag2_bin))
        # Magnitude in y axis.
        mag_bin = -2.5 * np.log10(10 ** (-0.4 * isoch_mass[1][bin_indxs]) +
                                  10 ** (-0.4 * bin_isoch[1]))
        # Transform back first filter's magnitude into color.
        if cmd_sel in {2, 5, 9}:
            col_bin = mag_bin - col_mag_bin
        else:
            col_bin = col_mag_bin - mag_bin

        # Add masses to obtain the system's mass.
        mass_bin = isoch_mass[2][bin_indxs] + bin_isoch[2]

        # Update array with new values of color, magnitude and masses.
        for indx, i in enumerate(bin_indxs):
            isoch_mass[0][i] = col_bin[indx]
            isoch_mass[1][i] = mag_bin[indx]
            isoch_mass[2][i] = mass_bin[indx]

    return isoch_mass

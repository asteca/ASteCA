
import numpy as np
import mass_interp


def mag_combine(m1, m2):
    """
    Combine two magnitudes. This is a faster re-ordering of the standard
    formula:

    -2.5 * np.log10(10 ** (-0.4 * m1) + 10 ** (-0.4 * m2))

    """
    c = 10 ** -.4
    return -2.5 * (-.4 * m1 + np.log10(1. + c ** (m2 - m1)))


def main(isoch_mass, isoch_cut, bin_frac, bin_mass_ratio, N_fc):
    '''
    Randomly select a fraction of stars to be binaries.
    '''
    bin_indxs = []
    # If the binary fraction is zero, skip the whole process.
    if bin_frac > 0.:

        # import time
        # t1, t2, t3, t4, t5, t6, t7, t8 = 0., 0., 0., 0., 0., 0., 0., 0.

        # Indexes of the randomly selected stars (without replacement) in
        # 'isoch_mass' to be converted to binary systems.
        # s = time.clock()
        bin_indxs = np.random.choice(
            len(isoch_mass[0]), int(bin_frac * len(isoch_mass[0])),
            replace=False)
        # t1 += time.clock() - s

        # s = time.clock()
        # Calculate the secondary masses of these binary stars between
        # bin_mass_ratio*m1 and m1, where m1 is the primary mass.
        # Index of m_ini (theoretical initial mass),
        # stored in the theoretical isochrones.
        m_ini = N_fc[0] + N_fc[1] + 2 * N_fc[1]
        # Draw the random 'bin_indxs' primary masses.
        m1 = np.asarray(isoch_mass[m_ini][bin_indxs])
        # Generate the random 'bin_indxs' secondary masses.
        m2 = np.random.uniform(bin_mass_ratio * m1, m1)
        # If any secondary mass falls outside of the lower isochrone's mass
        # range, change its value to the min value.
        m2 = np.maximum(np.min(isoch_mass[m_ini]), m2)
        # t2 += time.clock() - s

        # s = time.clock()
        # Find color and magnitude values for each secondary star. This will
        # slightly change the values of the masses, since they will be
        # assigned to the closest value found in the interpolated isochrone.
        bin_isoch = mass_interp.main(isoch_cut, m2, N_fc)
        # t3 += time.clock() - s

        # s = time.clock()
        # Calculate unresolved binary magnitude for each filter/magnitude
        # defined.
        mag_bin = []
        for i, m in enumerate(isoch_mass[:N_fc[0]]):
            mag_bin.append(mag_combine(m[bin_indxs], bin_isoch[i]))
        # t4 += time.clock() - s

        # s = time.clock()
        # Calculate unresolved color for each color defined.
        # Lower/upper index where filters that make up colors start/finish.
        fcl = N_fc[0] + N_fc[1]
        fcu = fcl + 2 * N_fc[1]
        # Extract only the portion that contains the color filters.
        filt_colors, bin_f_cols = isoch_mass[fcl:fcu], bin_isoch[fcl:fcu]
        # t5 += time.clock() - s

        # s = time.clock()
        # Filters composing the colors, i.e.: C = (f1 - f2).
        f1, f2 = [], []
        # The [::2] slice indicates even positions, starting from 0.
        for m, bin_m in zip(*[filt_colors[::2], bin_f_cols[::2]]):
            f1.append(mag_combine(m[bin_indxs], bin_m))
        # The [1::2] slice indicates odd positions.
        for m, bin_m in zip(*[filt_colors[1::2], bin_f_cols[1::2]]):
            f2.append(mag_combine(m[bin_indxs], bin_m))
        # t6 += time.clock() - s

        # s = time.clock()
        # Create the colors affected by binarity.
        col_bin = []
        for f_1, f_2 in zip(*[f1, f2]):
            col_bin.append(f_1 - f_2)
        # t7 += time.clock() - s

        # s = time.clock()
        # Update array with new values of magnitudes, colors, and masses.
        # New magnitudes.
        for i in range(N_fc[0]):
            isoch_mass[i][bin_indxs] = mag_bin[i]
        # New colors.
        for i in range(N_fc[1]):
            isoch_mass[N_fc[0] + i][bin_indxs] = col_bin[i]

        # Add masses to obtain the binary system's mass.
        isoch_mass[m_ini][bin_indxs] += bin_isoch[m_ini]
        # t8 += time.clock() - s

    return isoch_mass, bin_indxs
    # return np.array([t1, t2, t3, t4, t5, t6, t7, t8])


if __name__ == '__main__':

    N = 10000
    times = np.array([0., 0., 0., 0., 0., 0., 0., 0.])
    for _ in range(N):

        times += main()

    print("{:7.2f}:  {}".format(times.sum(), "  ".join(map(str, times))))

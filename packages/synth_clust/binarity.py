
import numpy as np
from .mass_interp import find_closest
from .. import update_progress


def main(isoch_mass, bin_frac, m_ini_idx, N_fc):
    """
    Update the randomly selected fraction of binary stars.
    """

    # If the binary fraction is zero, skip the whole process.
    if bin_frac > 0.:

        # Select a fraction of stars to be binaries, according to the random
        # probabilities assigned before.
        bin_indxs = isoch_mass[m_ini_idx - 2] <= bin_frac

        # Index of the first binary magnitude, stored in the theoretical
        # isochrones list.
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
        isoch_mass[m_ini_idx][bin_indxs] = isoch_mass[m_ini_idx - 1][bin_indxs]

        # Update binary systems to a '2.'.
        isoch_mass[m_ini_idx - 2][bin_indxs] = 2.

    return isoch_mass


def binarGen(
    binar_fracs, N_mass_interp, mags_theor, cols_theor, mags_cols_theor,
    extra_pars, all_met_vals, all_age_vals, bin_mass_ratio,
        synth_rand_seed):
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
    # First set the random seed in numpy
    from . import set_rand_seed
    set_rand_seed.main(synth_rand_seed)
    # Now import it
    from .set_rand_seed import np

    # If binary_fraction=0 don't bother obtaining the binary magnitudes,
    # colors, etc.
    if len(binar_fracs) > 1 or binar_fracs[0] > 0.:

        print("Generating binary data (b_mr={:.2f})".format(bin_mass_ratio))

        met_probs, age_probs, unq_probs, fracs = randVals(
            N_mass_interp, bin_mass_ratio)

        mags_binar, cols_binar, probs_binar, mass_binar = [], [], [], []
        # For each metallicity defined.
        N_mags_theor = len(mags_theor)
        for mx, _ in enumerate(mags_theor):

            mag_bin, col_bin, prob_bin, mass_bin = [], [], [], []
            # For each age defined.
            for ax, isoch in enumerate(_):

                # Extract initial masses for this isochrone. Assumes that
                # the initial masses are in the '0' index.
                mass_ini = extra_pars[mx][ax][0]

                # Calculate random secondary masses of these binary stars
                # between bin_mass_ratio*m1 and m1, where m1 is the primary
                # mass.
                # m2 = np.random.uniform(bin_mass_ratio * mass_ini, mass_ini)
                m2 = fracs * mass_ini

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

                # Find closest met & age values
                iz = np.searchsorted(met_probs, all_met_vals[mx])
                ia = np.searchsorted(age_probs, all_age_vals[ax])
                # This ensures that the same (z, a) pair points to the same
                # 'unq_probs' values (for the same random seed), no matter
                # the ranges used for these parameters.
                idx = (iz + ia) % len(unq_probs)
                probs = unq_probs[idx]
                prob_bin.append([probs])

            # Store for each metallicity value.
            mags_binar.append(mag_bin)
            cols_binar.append(col_bin)
            probs_binar.append(prob_bin)
            mass_binar.append(mass_bin)

            update_progress.updt(N_mags_theor, mx + 1)
    else:
        return None

    return mags_binar, cols_binar, probs_binar, mass_binar


def randVals(
    N_mass_interp, bin_mass_ratio, zmin=0., zmax=0.06, amin=6., amax=10.5,
        N_mets=50000, N_ages=50000, N_unq_probs=10000):
    """
    Process the required random values making sure that they are reproducible
    to the maximum possible extent. In the N_mets-->inf, N_ages-->inf limit
    all (z, a) pairs have a unique index assigned and hence always point to
    the exact same position in 'unq_probs'.

    The 'N_unq_probs' value is not that important. In the limit
    N_unq_probs-->inf every (z, a) pair has a unique array of probabilities
    assigned.

    HARDCODED
    zmin, zmax, amin, amax : full range for each parameter
    N_mets, N_ages, N_unq_probs : number of elements in each array
    """

    # All theoretical isochrones are interpolated with the same length,
    # assign unique binarity probabilities to each star randomly.
    b_probs = np.arange(N_mass_interp) / float(N_mass_interp)

    # As long as the length of 'b_probs' stays the same (i.e.:
    # as long as N_mass_interp stays the same), this will produce
    # the same results for the same random seed.

    met_probs = np.linspace(zmin, zmax, N_mets)
    age_probs = np.linspace(amin, amax, N_ages)
    unq_probs = []
    for _ in range(N_unq_probs):
        # Shuffle binarity probabilities.
        np.random.shuffle(b_probs)
        unq_probs.append(list(b_probs))

    # Fractions for second mass
    fracs = np.random.uniform(bin_mass_ratio, 1., N_mass_interp)

    return met_probs, age_probs, unq_probs, fracs


def mag_combine(m1, m2):
    """
    Combine two magnitudes. This is a faster re-ordering of the standard
    formula:

    -2.5 * np.log10(10 ** (-0.4 * m1) + 10 ** (-0.4 * m2))

    """

    c = 10 ** -.4
    mbin = -2.5 * (-.4 * m1 + np.log10(1. + c ** (m2 - m1)))

    return mbin

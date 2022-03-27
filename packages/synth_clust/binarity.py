
import numpy as np
from scipy import stats
from .. import update_progress


def main(isoch_mass, bin_frac, m_ini_idx, N_fc, binar_probs):
    """
    Update the randomly selected fraction of binary stars.
    """

    # If the binary fraction is zero, skip the whole process.
    if bin_frac > 0.:

        #
        Nm = isoch_mass[m_ini_idx].size
        b_p = binar_mass(isoch_mass[m_ini_idx])
        bin_indxs = np.random.choice(range(Nm), int(Nm * bin_frac), False, b_p)

        # Select a fraction of stars to be binaries, according to the random
        # probabilities assigned before.
        # bin_indxs = binar_probs[:isoch_mass.shape[-1]] <= bin_frac

        # Index of the first binary magnitude.
        # mag_binar = m_ini_idx + 1
        # Update array with new values of magnitudes and colors.
        for i in range(N_fc[0] + N_fc[1]):
            isoch_mass[i][bin_indxs] = isoch_mass[m_ini_idx + 1 + i][bin_indxs]

        # Update the binary systems' masses so that the secondary masses for
        # SINGLE systems are identified with a '-99.' value.
        # isoch_mass[-1][~bin_indxs] = -99.

        sing_indxs = list(set(range(Nm)) - set(bin_indxs))
        isoch_mass[-1][sing_indxs] = -99.

    return isoch_mass


def binar_mass(x):
    """
    https://stackoverflow.com/a/29359275/1391441
    """
    from scipy.interpolate import interp1d

    # xx = (0.08, .29, .96, 2.4, 7.65, 28.5, 70)
    # yy = (.22, .26, .46, .52, .64, .82, 1)
    # logx = np.log10(xx)
    # logy = np.log10(yy)
    logx = (-1.097, -0.538, -0.018, 0.380, 0.884, 1.455, 1.845)
    logy = (-0.658, -0.585, -0.337, -0.284, -0.194, -0.086, 0.)
    lin_interp = interp1d(logx, logy)
    b_p = np.power(10.0, lin_interp(np.log10(x)))
    return b_p / b_p.sum() # <-- IS THIS NORMALIZATION REASONABLE


def binarGen(
    min_bmass_ratio, m_ini_idx, N_fc, interp_tracks, mags_cols_intp,
        all_met_vals, all_age_vals):
    """
    For each theoretical isochrone defined.
        1. Draw random secondary masses for *all* stars.
        2. Interpolate magnitude values for these masses
        3. Calculate the magnitudes for the binary system.
        4. Calculate the colors for the binary system.
    """
    from .set_rand_seed import np

    print("Generating binary data (b_mr=[{:.2f}, 1.])".format(
        min_bmass_ratio))

    # interp_tracks.shape = (Nz, Na, Nd, Np)
    N_mets, Na, Nd, N_mass = interp_tracks.shape

    # Extend to accommodate binary data
    interp_tracks = np.concatenate((
        interp_tracks, np.zeros([N_mets, Na, Nd, N_mass])), axis=2)

    # Fractions for second mass, one per metallicity
    mass_ratios, mean_bin_mr = randBinarFracs(min_bmass_ratio, N_mass, N_mets)

    # For each metallicity defined.
    for mx, _ in enumerate(interp_tracks):
        # For each age defined.
        for ax, isoch in enumerate(_):

            # Extract initial masses for this isochrone.
            mass_ini = isoch[m_ini_idx]

            # Calculate random secondary masses of these binary stars
            # between bin_mass_ratio*m1 and m1, where m1 is the primary
            # mass.
            m2 = mass_ratios[mx] * mass_ini

            # Calculate unresolved binary magnitude for each
            # filter/magnitude defined.
            for i, mag in enumerate(isoch[:N_fc[0]]):
                mag_m2 = np.interp(m2, mass_ini, mag)
                interp_tracks[mx][ax][m_ini_idx + 1 + i] = mag_combine(
                    mag, mag_m2)

            # Calculate unresolved color for each color defined.
            # The [::2] slice indicates even positions, starting from 0.
            even_f_cols = mags_cols_intp[mx][ax][::2]
            # The [1::2] slice indicates odd positions.
            odd_f_cols = mags_cols_intp[mx][ax][1::2]

            # Filters composing the colors, i.e.: C = (f1 - f2).
            f1, f2 = [], []
            for i, m in enumerate(even_f_cols):
                f_m2 = np.interp(m2, mass_ini, m)
                f1.append(mag_combine(m, f_m2))
            for i, m in enumerate(odd_f_cols):
                f2_m2 = np.interp(m2, mass_ini, m)
                f2.append(mag_combine(m, f2_m2))

            # Create the colors affected by binarity.
            for i, (f_1, f_2) in enumerate(zip(*[f1, f2])):
                interp_tracks[mx][ax][m_ini_idx + 1 + N_fc[0] + i] =\
                    f_1 - f_2

            # Secondary masses
            interp_tracks[mx][ax][-1] = m2

            # import matplotlib.pyplot as plt
            # # plt.subplot(121)
            # plt.title("Gaia eDR3")
            # plt.scatter(mags_cols_intp[mx][ax][0] - mags_cols_intp[mx][ax][1],
            #             isoch[0], c='g')
            # mag_bin0 = interp_tracks[mx][ax][m_ini_idx + 1]
            # col_bin0 = interp_tracks[mx][ax][m_ini_idx + 1 + N_fc[0] + 0]
            # plt.scatter(col_bin0, mag_bin0, c='r', alpha=.5)
            # plt.gca().invert_yaxis()
            # # # Second color
            # # plt.subplot(122)
            # # col_bin1 = interp_tracks[mx][ax][m_ini_idx + 1 + N_fc[0] + 1]
            # # plt.scatter(mags_cols_intp[mx][ax][2] - mags_cols_intp[mx][ax][3],
            # #             isoch[0], c='g')
            # # plt.scatter(col_bin1, mag_bin0, c='r', alpha=.5)
            # # plt.gca().invert_yaxis()
            # plt.show()

        update_progress.updt(N_mets, mx + 1)

    return interp_tracks, mean_bin_mr


def randBinarFracs(mbr, N_mass, N_mets):
    """
    IN PLACE FOR #496

    Define the mass ratio for the secondary masses
    """

    mbr = 'raghavan'
    # mbr = 'fisher'

    def fQ(xk, pk):
        """
        Discrete function
        """
        pk /= pk.sum()
        fq = stats.rv_discrete(a=0., b=1., name='custm', values=(xk, pk))
        # Average minimum mass fraction for binary systems
        mean_bin_mr = fq.mean()

        return fq, mean_bin_mr

    # Fisher's distribution
    if mbr == 'fisher':
        # Fisher, Schröder & Smith (2005), 10.1111/j.1365-2966.2005.09193.x;
        # Fig 6 (solid line)
        xk = np.arange(0.05, 1.01, 0.1)
        pk = np.array([
            0.07320166320166316, 0.0811413721413721, 0.0924109494109494,
            0.08995495495495492, 0.07918849618849616, 0.07049757449757446,
            0.07344906444906443, 0.08970755370755368, 0.11137075537075536,
            0.24156202356202355])
        fq, mean_bin_mr = fQ(xk, pk)
    elif mbr == 'raghavan':
        # Raghavan et al. (2010), 10.1088/0067-0049/190/1/1; Fig 16 (left)
        xk = np.arange(0.025, 1.01, 0.05)
        pk = np.array([
            0.5263157894736867, 2.6105263157894765, 0.5263157894736867,
            4.673684210526318, 7.810526315789476, 3.6421052631578963,
            9.894736842105264, 5.7052631578947395, 4.694736842105266,
            5.726315789473686, 4.673684210526318, 6.757894736842108,
            5.747368421052634, 5.726315789473686, 2.6105263157894765,
            5.7052631578947395, 4.715789473684213, 5.7052631578947395,
            3.6421052631578963, 12.989473684210529])
        fq, mean_bin_mr = fQ(xk, pk)
    else:
        mean_bin_mr = (mbr + 1.) / 2.

    mass_ratios = []
    for _ in range(N_mets):
        if mbr in ('fisher', 'raghavan'):
            # 'ppf' is the inverse CDF
            dist = fq.ppf(np.random.uniform(0., 1., N_mass))
        else:
            dist = np.random.uniform(mbr, 1., N_mass)
        mass_ratios.append(dist)

    # import matplotlib.pyplot as plt
    # plt.hist(dist, 20);plt.show()

    return mass_ratios, mean_bin_mr


def mag_combine(m1, m2):
    """
    Combine two magnitudes. This is a faster re-ordering of the standard
    formula:

    -2.5 * np.log10(10 ** (-0.4 * m1) + 10 ** (-0.4 * m2))

    """

    c = 10 ** -.4
    mbin = -2.5 * (-.4 * m1 + np.log10(1. + c ** (m2 - m1)))

    return mbin


# DEPRECATED FEB 2021
# def randVals(N_mass_interp, N_unq_probs=1000):
#     """
#     Process the required random values making sure that they are reproducible
#     to the maximum possible extent.

#     In the limit N_unq_probs-->inf every (z, a) pair has a unique array of
#     probabilities assigned.

#     HARDCODED
#     N_unq_probs : number of XXXX
#     """

#     # All theoretical isochrones are interpolated with the same length,
#     # assign unique binarity probabilities to each star randomly.
#     b_probs = np.arange(N_mass_interp) / float(N_mass_interp)

#     # As long as the length of 'b_probs' stays the same (i.e.:
#     # as long as N_mass_interp stays the same), this will produce
#     # the same results for the same random seed.

#     unq_probs = []
#     for _ in range(N_unq_probs):
#         # Shuffle binarity probabilities.
#         np.random.shuffle(b_probs)
#         unq_probs.append(list(b_probs))

#     return unq_probs

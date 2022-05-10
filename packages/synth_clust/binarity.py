
import numpy as np


def main(isoch_mass, alpha, beta, m_ini_idx, N_fc, rand_unif_vals):
    """
    Select a fraction of stars to be binaries, given a chosen method.
    """

    Ns = isoch_mass.shape[-1]
    mass = isoch_mass[m_ini_idx]

    b_p = binarProbsF(mass, alpha, beta)
    # Stars (masses) with the largest binary probabilities are selected
    # proportional to their probability
    bin_indxs = b_p > rand_unif_vals[:Ns]

    # Index of the first binary magnitude.: mag_binar = m_ini_idx + 1
    # Update array with new values of magnitudes and colors for the binary
    # systems.
    if bin_indxs.any():
        for i in range(N_fc[0] + N_fc[1]):
            isoch_mass[i][bin_indxs] = isoch_mass[m_ini_idx + 1 + i][bin_indxs]

    # Update the binary systems' masses so that the secondary masses for
    # SINGLE systems are identified with a '0.' value.
    isoch_mass[-1][~bin_indxs] = 0.

    return isoch_mass


def binarGen(
    gamma, m_ini_idx, N_fc, interp_tracks, mags_cols_intp,
        all_met_vals, all_age_vals):
    """
    For each theoretical isochrone defined.
        1. Draw random secondary masses for *all* stars.
        2. Interpolate magnitude values for these masses
        3. Calculate the magnitudes for the binary system.
        4. Calculate the colors for the binary system.
    """
    from .set_rand_seed import np

    # interp_tracks.shape = (Nz, Na, Nd, Np)
    N_mets, Na, Nd, N_mass = interp_tracks.shape

    # Extend to accommodate binary data
    interp_tracks = np.concatenate((
        interp_tracks, np.zeros([N_mets, Na, Nd, N_mass])), axis=2)

    # For each metallicity defined.
    for mx, _ in enumerate(interp_tracks):
        # For each age defined.
        for ax, isoch in enumerate(_):

            # Extract initial masses for this isochrone.
            mass_ini = isoch[m_ini_idx]

            # Mass-ratio distribution
            mass_ratios = qDistribution(mass_ini, gamma)

            # Calculate random secondary masses of these binary stars
            # between bin_mass_ratio*m1 and m1, where m1 is the primary
            # mass.
            m2 = mass_ratios * mass_ini

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

    return interp_tracks


def qDistribution(M1, gamma):
    """
    D&K      : Distribution of mass-ratios versus primary masses
    (Duchene & Kraus 2013). Mass dependent.
    powerlaw : Power-law distribution with shape parameter 'gamma'. Not mass
    dependent.


    Use 'gamma + 1' in the power-law distribution below because in D&K this
    distribution is defined as f(q)~q^gamma, while numpy's distribution is
    defined as ~a*x^(a-1).
    """
    try:
        gamma = float(gamma)
        # mass_ratios = powerlaw.rvs(gamma, size=M1.size)
        mass_ratios = np.random.power(gamma + 1, M1.size)
    except ValueError:
        msk1, gamma1 = M1 <= 0.1, 4.2
        msk2, gamma2 = (M1 > 0.1) & (M1 <= 0.6), 0.4
        msk3, gamma3 = (M1 > 0.6) & (M1 <= 1.4), 0.3
        msk4, gamma4 = (M1 > 1.4) & (M1 <= 6.5), -0.5
        msk5, gamma5 = (M1 > 6.5) & (M1 <= 16), 0.0  # <- Not sure. Use uniform
        msk6, gamma6 = M1 > 16, 0.0  # <- Not sure. Use uniform

        mass_ratios = np.zeros(M1.size)
        for msk, gamma in (
                (msk1, gamma1), (msk2, gamma2), (msk3, gamma3), (msk4, gamma4),
                (msk5, gamma5), (msk6, gamma6)):
            q = np.random.power(gamma + 1, msk.sum())
            mass_ratios[msk] = q

    mass_ratios = np.clip(mass_ratios, a_min=0, a_max=1)

    return mass_ratios


def binarProbsF(x, alpha=None, beta=None):
    """
    Distribution of probability of binarity (multiplicity fraction) versus
    primary masses.

    D&K: from Duchene and Kraus 2013.
    Source: https://stackoverflow.com/a/29359275/1391441
    """
    # if bp_vs_mass == "D&K":
    #     # xx = (0.08, .29, .96, 2.4, 7.65, 28.5, 151)
    #     # yy = (.22, .26, .46, .52, .64, .82, 1)
    #     # logx = np.log10(xx)
    #     # logy = np.log10(yy)
    #     logx = (-1.097, -0.538, -0.018, 0.380, 0.884, 1.455, 2.18)
    #     logy = (-0.658, -0.585, -0.337, -0.284, -0.194, -0.086, 0.)
    #     lin_interp = interp1d(logx, logy)
    #     b_p = np.power(10.0, lin_interp(np.log10(x)))

    # elif bp_vs_mass == "logfit":
    b_p = alpha + beta * np.log(x + 1)
    b_p = np.clip(b_p, a_min=0, a_max=1)

    return b_p


def mag_combine(m1, m2):
    """
    Combine two magnitudes. This is a faster re-ordering of the standard
    formula:

    -2.5 * np.log10(10 ** (-0.4 * m1) + 10 ** (-0.4 * m2))

    """

    c = 10 ** -.4
    mbin = -2.5 * (-.4 * m1 + np.log10(1. + c ** (m2 - m1)))

    return mbin


# def randBinarFracs(mbr, N_mass, N_mets):
#     """
#     IN PLACE FOR #496

#     Define the mass ratio for the secondary masses
#     """

#     mbr = 'raghavan'
#     # mbr = 'fisher'

#     def fQ(xk, pk):
#         """
#         Discrete function
#         """
#         pk /= pk.sum()
#         fq = stats.rv_discrete(a=0., b=1., name='custm', values=(xk, pk))
#         # Average minimum mass fraction for binary systems
#         mean_bin_mr = fq.mean()

#         return fq, mean_bin_mr

#     # Fisher's distribution
#     if mbr == 'fisher':
#         # Fisher, Schröder & Smith (2005), 10.1111/j.1365-2966.2005.09193.x;
#         # Fig 6 (solid line)
#         xk = np.arange(0.05, 1.01, 0.1)
#         pk = np.array([
#             0.07320166320166316, 0.0811413721413721, 0.0924109494109494,
#             0.08995495495495492, 0.07918849618849616, 0.07049757449757446,
#             0.07344906444906443, 0.08970755370755368, 0.11137075537075536,
#             0.24156202356202355])
#         fq, mean_bin_mr = fQ(xk, pk)
#     elif mbr == 'raghavan':
#         # Raghavan et al. (2010), 10.1088/0067-0049/190/1/1; Fig 16 (left)
#         xk = np.arange(0.025, 1.01, 0.05)
#         pk = np.array([
#             0.5263157894736867, 2.6105263157894765, 0.5263157894736867,
#             4.673684210526318, 7.810526315789476, 3.6421052631578963,
#             9.894736842105264, 5.7052631578947395, 4.694736842105266,
#             5.726315789473686, 4.673684210526318, 6.757894736842108,
#             5.747368421052634, 5.726315789473686, 2.6105263157894765,
#             5.7052631578947395, 4.715789473684213, 5.7052631578947395,
#             3.6421052631578963, 12.989473684210529])
#         fq, mean_bin_mr = fQ(xk, pk)
#     else:
#         mean_bin_mr = (mbr + 1.) / 2.

#     mass_ratios = []
#     for _ in range(N_mets):
#         if mbr in ('fisher', 'raghavan'):
#             # 'ppf' is the inverse CDF
#             dist = fq.ppf(np.random.uniform(0., 1., N_mass))
#         else:
#             dist = np.random.uniform(mbr, 1., N_mass)
#         mass_ratios.append(dist)

#     # import matplotlib.pyplot as plt
#     # plt.hist(dist, 20);plt.show()

#     return mass_ratios, mean_bin_mr


import numpy as np


def ccm_model(mw):
    '''
    Cardelli, Clayton, and Mathis (1989 ApJ. 345, 245) model for extinction
    coefficients with updated coefficients for near-UV from O'Donnell (1994).

    ccm_coef = a + b / Rv

    Implementation taken from:

    http://idlastro.gsfc.nasa.gov/ftp/pro/astro/ccm_unred.pro

    There appears to be an error in the Far-UV range in the original IDL
    routine where the maximum inverse wavelength is 11 and it should be 10
    according to Cardelli et al. 1989 (pag 251, Eq (5,a,b)).
    '''

    if 0.3 <= mw < 1.1:
        # Infrared.
        a, b = 0.574 * (mw ** 1.61), -0.527 * (mw ** 1.61)

    elif 1.1 <= mw < 3.3:
        # Optical/NIR.
        # Original coefficients from CCM89
        # c1 = [1., 0.17699, -0.50447, -0.02427, 0.72085, 0.01979, -0.77530,
        #       0.32999]
        # c2 = [0., 1.41338, 2.28305, 1.07233, -5.38434, -0.62251, 5.30260,
        #       -2.09002]
        # New coefficients from O'Donnell (1994)
        c1 = [1., 0.104, -0.609, 0.701, 1.137, -1.718, -0.827, 1.647, -0.505]
        c2 = [0., 1.952, 2.908, -3.989, -7.985, 11.102, 5.491, -10.805, 3.347]
        y = mw - 1.82
        # Reverse because polyval starts from the highest degree.
        c1.reverse(), c2.reverse()
        a, b = np.polyval(c1, y), np.polyval(c2, y)

    elif 3.3 <= mw < 8.:
        # Mid-UV
        F_a, F_b = 0., 0.
        if mw >= 5.9:
            y = mw - 5.9
            F_a = -0.04473 * y ** 2 - 0.009779 * y ** 3
            F_b = 0.2130 * y ** 2 + 0.1207 * y ** 3
        a = 1.752 - 0.316 * mw - (0.104 / ((mw - 4.67) ** 2 + 0.341)) + F_a
        b = -3.090 + 1.825 * mw + (1.206 / ((mw - 4.62) ** 2 + 0.263)) + F_b

    elif 8. <= mw <= 10.:
        # Far-UV
        c1 = [-1.073, -0.628, 0.137, -0.070]
        c2 = [13.670, 4.257, -0.420, 0.374]
        y = mw - 8.
        c1.reverse(), c2.reverse()
        a, b = np.polyval(c1, y), np.polyval(c2, y)
    else:
        raise ValueError(
            "The effective wavelength is {} [1/micron], beyond "
            "the CCM model limit (10 [1/micron]).".format(mw))

    return a, b


def main(cmd_systs, filters, colors, ext_shape):
    """
    Obtain extinction coefficients for all the observed filters and colors,
    in the order in which they are stored in theor_tracks.
    """
    a_ext_coefs, b_ext_coefs = [], []

    # For filters.
    for f in filters:
        # f[0]: photometric system; f[1]: name of filter
        # Index of filter
        fi = cmd_systs[f[0]][1].index(f[1])
        # Effective wavelength in Armstrong.
        eff_wave = cmd_systs[f[0]][2][fi]
        # CCM coefficient for this filter. Use the effective wavelength in
        # inverse microns.
        a, b = ccm_model(10000. / eff_wave)
        a_ext_coefs.append(a)
        b_ext_coefs.append(b)

    # For colors.
    for c in colors:
        # c[0]: photometric system; c[1]: color index
        # Index of filters.
        ci1 = cmd_systs[c[0]][1].index(c[1].split(',')[0])
        ci2 = cmd_systs[c[0]][1].index(c[1].split(',')[1])
        # Effective wavelength in Armstrong.
        eff_wave1 = cmd_systs[c[0]][2][ci1]
        eff_wave2 = cmd_systs[c[0]][2][ci2]
        # CCM coefficients for this color.
        a1, b1 = ccm_model(10000. / eff_wave1)
        a2, b2 = ccm_model(10000. / eff_wave2)
        a_ext_coefs.append(a1 - a2)
        b_ext_coefs.append(b1 - b2)

    # Double for binary mags and colors.
    a_ext_coefs = a_ext_coefs + a_ext_coefs
    b_ext_coefs = b_ext_coefs + b_ext_coefs
    # Add zeros to account for the shape of the isochrones and the extra
    # parameters present.
    a_ext_coefs = a_ext_coefs + [0.] * (ext_shape - len(a_ext_coefs))
    b_ext_coefs = b_ext_coefs + [0.] * (ext_shape - len(b_ext_coefs))

    # Arrays of zeros except for the positions of the magnitudes, to use
    # for multiplying with the distance modulus and move the magnitudes.
    dm_arr = np.zeros(len(a_ext_coefs))
    dm_arr[0], dm_arr[len(filters) + len(colors)] = 1., 1.

    # Store in the proper format for move_isochrone().
    ext_coefs = [np.vstack(np.array(_)) for _ in [
                 a_ext_coefs, b_ext_coefs, dm_arr]]

    return ext_coefs

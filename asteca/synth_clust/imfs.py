import numpy as np


def get_imf(IMF_name, m_star):
    """
    Define any number of IMFs.

    The package https://github.com/keflavich/imf has some more (I think,
    24-09-2019).
    """
    if IMF_name == "salpeter_1955":
        # Salpeter (1955)  IMF.
        # https://ui.adsabs.harvard.edu/abs/1955ApJ...121..161S/
        imf_val = m_star**-2.35

    elif IMF_name == "kroupa_1993":
        # Kroupa, Tout & Gilmore. (1993) piecewise IMF.
        # http://adsabs.harvard.edu/abs/1993MNRAS.262..545K
        # Eq. (13), p. 572 (28)
        alpha = [-1.3, -2.2, -2.7]
        m0, m1, m2 = [0.08, 0.5, 1.0]
        factor = [0.035, 0.019, 0.019]
        if m0 < m_star <= m1:
            i = 0
        elif m1 < m_star <= m2:
            i = 1
        elif m2 < m_star:
            i = 2
        imf_val = factor[i] * (m_star ** alpha[i])

    elif IMF_name == "kroupa_2001":
        # Kroupa (2001), Mon. Not. R. Astron. Soc. 322, 231-246 (2001); Eq. 2
        alpha = [-0.3, -1.3, -2.3]
        m0, m1, m2 = [0.01, 0.08, 0.5]
        if m0 < m_star <= m1:
            i = 0
        elif m1 < m_star <= m2:
            i = 1
        elif m2 < m_star:
            i = 2
        imf_val = m_star ** alpha[i]

    elif IMF_name == "chabrier_2001_log":
        # Chabrier (2001) lognormal form of the IMF.
        # http://adsabs.harvard.edu/abs/2001ApJ...554.1274C
        # Eq (7)
        imf_val = (
            (1.0 / (np.log(10) * m_star))
            * 0.141
            * np.exp(-((np.log10(m_star) - np.log10(0.1)) ** 2) / (2 * 0.627**2))
        )

    elif IMF_name == "chabrier_2001_exp":
        # Chabrier (2001) exponential form of the IMF.
        # http://adsabs.harvard.edu/abs/2001ApJ...554.1274C
        # Eq (8)
        imf_val = 3.0 * m_star ** (-3.3) * np.exp(-((716.4 / m_star) ** 0.25))

    elif IMF_name == "popescu_2009":
        # Kroupa (2002) Salpeter (1995) piecewise IMF taken from Popescu & Hanson
        # (2009; MASSCLEAN), Eq. (2) & (3), p. 1725
        alpha = [-0.3, -1.3, -2.3]
        m0, m1, m2 = [0.01, 0.08, 0.5]
        factor = [
            (1.0 / m1) ** alpha[0],
            (1.0 / m1) ** alpha[1],
            ((m2 / m1) ** alpha[1]) * ((1.0 / m2) ** alpha[2]),
        ]
        if m0 < m_star <= m1:
            i = 0
        elif m1 < m_star <= m2:
            i = 1
        elif m2 < m_star:
            i = 2
        imf_val = factor[i] * (m_star ** alpha[i])

    return imf_val

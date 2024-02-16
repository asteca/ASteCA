import numpy as np
from scipy.interpolate import make_interp_spline


def invTrnsfSmpl(IMF_name, m_low=0.08, m_high=150, mass_step=0.05):
    """
    IMF inverse transform sampling.

    Asked here: https://stackoverflow.com/q/21100716/1391441
    """
    # IMF mass interpolation step and grid values.
    mass_values = np.arange(m_low, m_high, mass_step)

    def IMF_func(m_star, IMF_name):
        return get_imf(IMF_name, m_star)

    # CDF_samples = []
    # for m in mass_values:
    #     x = quad(IMF_func, m_low, m, args=(IMF_name))[0]
    #     CDF_samples.append(x)
    # # Normalize values
    # CDF_samples = np.array(CDF_samples) / max(CDF_samples)
    # # These are (0, 1)
    # # CDF_min, CDF_max = CDF_samples.min(), CDF_samples.max()
    # # Inverse CDF
    # from scipy.interpolate import interp1d
    # inv_cdf = interp1d(CDF_samples, mass_values)  # Almost deprecated

    # x_old, CDF_samples = -np.inf, []
    # for m in mass_values:
    #     x = quad(IMF_func, m_low, m, args=(IMF_name))[0]
    #     # This ensures that CDF_samples is monotonically increasing, bypassing
    #     # rounding errors
    #     if x <= x_old:
    #         x = x_old + 1e-15
    #     CDF_samples.append(x)
    #     x_old = x

    # # Normalize values
    # CDF_samples = np.array(CDF_samples) / max(CDF_samples)
    # # These are (0, 1)
    # # CDF_min, CDF_max = CDF_samples.min(), CDF_samples.max()

    # # k=1 is important otherwise the interpolator can become unstable for values
    # # close to 1 and return huge masses and even negative ones
    # inv_cdf = make_interp_spline(CDF_samples, mass_values, k=1)

    # The lower mass region needs to be sampled more accurately
    mass_values = list(np.linspace(m_low, .5, 50))
    mass_values += list(np.linspace(.501, 10, 250))
    mass_values += list(np.linspace(10.01, m_high, 250))

    CDF_samples = []
    IMF_old, m_old, area_CDF = IMF_func(m_low, IMF_name), m_low, 0.
    for m in mass_values[1:]:
        # Approximate integral with rectangular area, and add to previous total area
        IMF_new = IMF_func(m, IMF_name)
        area_CDF += .5*(IMF_new + IMF_old) * (m - m_old)
        CDF_samples.append(area_CDF)
        IMF_old, m_old = IMF_new, m
    CDF_samples = np.array(CDF_samples) / max(CDF_samples)

    # k=1 is important otherwise the interpolator can become unstable for values
    # close to 1 and return huge masses and even negative ones
    inv_cdf = make_interp_spline(CDF_samples, mass_values[1:], k=1)

    return inv_cdf


def sampleInv(Max_mass, inv_cdf):
    """
    Sample the inverse CDF up to `Max_mass`
    """
    def sampled_inv_cdf(N):
        mr = np.random.rand(N)
        return inv_cdf(mr)

    # Sample in chunks until the maximum defined mass is reached.
    N_chunk = max(100, int(2.5 * Max_mass / 100.0))
    Mass_tot, mass_samples = 0, []
    while Mass_tot < Max_mass:
        masses = sampled_inv_cdf(N_chunk).tolist()
        mass_samples += masses
        Mass_tot += sum(masses)
    sampled_IMF = np.array(mass_samples)

    return sampled_IMF


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

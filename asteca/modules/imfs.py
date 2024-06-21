import numpy as np
from scipy.interpolate import make_interp_spline


def invTrnsfSmpl(IMF_name, m_low=0.08, m_high=100):
    """
    IMF inverse transform sampling.
    """

    # The lower mass region needs to be sampled more accurately
    mass_values = list(np.linspace(m_low, 0.5, 100))
    mass_values += list(np.linspace(0.501, 2, 75))
    mass_values += list(np.linspace(2.01, 10, 50))
    mass_values += list(np.linspace(10.01, m_high, 25))

    CDF_samples = []
    IMF_old, m_old, area_CDF = get_imf(IMF_name, m_low), m_low, 0.0
    for m in mass_values[1:]:
        # Approximate integral with rectangular area, and add to previous total area
        IMF_new = get_imf(IMF_name, m)
        area_CDF += 0.5 * (IMF_new + IMF_old) * (m - m_old)
        CDF_samples.append(area_CDF)
        IMF_old, m_old = IMF_new, m
    CDF_samples = np.array(CDF_samples) / max(CDF_samples)

    # k=1 is important otherwise the interpolator can become unstable for values
    # close to 1 and return huge masses and even negative ones
    inv_cdf = make_interp_spline(CDF_samples, mass_values[1:], k=1)

    return inv_cdf


def sampleInv(ijseed, Max_mass, inv_cdf):
    """
    Sample the inverse CDF up to `Max_mass`
    """

    def sampled_inv_cdf(ijkseed, N):
        # mr = np.random.rand(N)
        # The internal seed can not be 'self.seed' because otherwise all the
        # floats returned are equal
        mr = np.random.default_rng(ijkseed).uniform(0.0, 1.0, N)
        return inv_cdf(mr)

    # Sample in chunks until the maximum defined mass is reached.
    N_chunk = max(100, int(2.5 * Max_mass / 100.0))
    k, Mass_tot, mass_samples = 0, 0, []
    int_seed = np.random.default_rng(ijseed).integers(1)
    while Mass_tot < Max_mass:
        masses = sampled_inv_cdf(k + int_seed + ijseed, N_chunk).tolist()
        mass_samples += masses
        Mass_tot += sum(masses)
        k += 1
    sampled_IMF = np.array(mass_samples)

    return sampled_IMF


def get_imf(IMF_name, m_star):
    """
    Define any number of IMFs.

    The package https://github.com/keflavich/imf has some more (I think).
    """
    if IMF_name == "salpeter_1955":
        # Salpeter (1955)  IMF.
        # https://ui.adsabs.harvard.edu/abs/1955ApJ...121..161S/
        imf_val = m_star**-2.35

    elif IMF_name == "kroupa_2001":
        # Kroupa (2001), Mon. Not. R. Astron. Soc. 322, 231-246 (2001); Eq. 2
        # https://ui.adsabs.harvard.edu/abs/2001MNRAS.322..231K/
        alpha = [-0.3, -1.3, -2.3]
        m0, m1, m2 = [0.01, 0.08, 0.5]
        # Continuity factor taken from Popescu & Hanson (2009; MASSCLEAN),
        # Eq. (2) & (3), p. 1725
        factor = [
            (1.0 / m1) ** alpha[0],
            (1.0 / m1) ** alpha[1],
            ((m2 / m1) ** alpha[1]) * ((1.0 / m2) ** alpha[2]),
        ]
        if m0 <= m_star <= m1:
            i = 0
        elif m1 < m_star <= m2:
            i = 1
        elif m2 < m_star:
            i = 2
        imf_val = factor[i] * (m_star ** alpha[i])

    elif IMF_name == "chabrier_2014":
        # Chabrier et al. (2014)
        # https://ui.adsabs.harvard.edu/abs/2014ApJ...796...75C/ ; Eq (34)
        nc, mc = 11, 0.18
        m0 = nc * mc
        Ah, x = 0.649, 1.35
        Al = Ah * nc ** (x / 2)
        c = 0.434294 / m_star  # np.log10(e)/m ; This is the transformation from
        # dN/log(m) --> dN/dm
        sigma_2 = np.log10(nc) / (x * np.log(10))
        if m_star <= m0:
            imf_val = (
                c
                * Al
                * m0 ** (-x)
                * np.exp(-((np.log10(m_star) - np.log10(mc)) ** 2) / (2 * sigma_2))
            )
        elif m_star > m0:
            imf_val = c * Ah * m_star ** (-x)

    return imf_val

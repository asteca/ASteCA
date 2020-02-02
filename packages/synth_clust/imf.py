
import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d


def main(IMF_name, masses, m_high=150.):
    """
    Returns the number of stars per interval of mass for the selected IMF.

    Parameters
    ----------
    IMF_name : str
      Name of the IMF to be used.
    masses: array
      Array of floats containing the range of masses defined.
    m_high : float
      Maximum mass value to be sampled.

    Returns
    -------
    st_dist_mass : dict
      Dictionary that contains the N defined masses, one array per mass value.
      Each array contains a given number of stars sampled from the selected
      IMF, that (approximately) sum to the associated total mass.

    """
    from .set_rand_seed import np

    print("Sampling selected IMF ({})".format(IMF_name))

    # Low mass limits for each IMF. Defined slightly larger to avoid sampling
    # issues.
    imfs_dict = {
        'chabrier_2001_exp': 0.011, 'chabrier_2001_log': 0.011,
        'kroupa_1993': 0.081, 'kroupa_2002': 0.011, 'salpeter_1955': 0.31}
    # IMF low mass limit.
    m_low = imfs_dict[IMF_name]

    sampled_IMF = invTrnsfSmpl(m_low)

    st_dist_mass = (sampled_IMF, np.cumsum(sampled_IMF))

    return st_dist_mass


def IMF_func(m_star, IMF_name):
    return imfs(IMF_name, m_star)


def imfs(IMF_name, m_star):
    """
    Define any number of IMFs.

    The package https://github.com/keflavich/imf has some more (I think,
    24-09-2019).
    """

    if IMF_name == 'kroupa_1993':
        # Kroupa, Tout & Gilmore. (1993) piecewise IMF.
        # http://adsabs.harvard.edu/abs/1993MNRAS.262..545K
        # Eq. (13), p. 572 (28)
        alpha = [-1.3, -2.2, -2.7]
        m0, m1, m2 = [0.08, 0.5, 1.]
        factor = [0.035, 0.019, 0.019]
        if m0 < m_star <= m1:
            i = 0
        elif m1 < m_star <= m2:
            i = 1
        elif m2 < m_star:
            i = 2
        imf_val = factor[i] * (m_star ** alpha[i])

    elif IMF_name == 'kroupa_2002':
        # Kroupa (2002) Salpeter (1995) piecewise IMF taken from MASSCLEAN
        # article, Eq. (2) & (3), p. 1725
        alpha = [-0.3, -1.3, -2.3]
        m0, m1, m2 = [0.01, 0.08, 0.5]
        factor = [(1. / m1) ** alpha[0], (1. / m1) ** alpha[1],
                  ((m2 / m1) ** alpha[1]) * ((1. / m2) ** alpha[2])]
        if m0 < m_star <= m1:
            i = 0
        elif m1 < m_star <= m2:
            i = 1
        elif m2 < m_star:
            i = 2
        imf_val = factor[i] * (m_star ** alpha[i])

    elif IMF_name == 'chabrier_2001_log':
        # Chabrier (2001) lognormal form of the IMF.
        # http://adsabs.harvard.edu/abs/2001ApJ...554.1274C
        # Eq (7)
        imf_val = (1. / (np.log(10) * m_star)) * 0.141 * \
            np.exp(-((np.log10(m_star) - np.log10(0.1)) ** 2) /
                   (2 * 0.627 ** 2))

    elif IMF_name == 'chabrier_2001_exp':
        # Chabrier (2001) exponential form of the IMF.
        # http://adsabs.harvard.edu/abs/2001ApJ...554.1274C
        # Eq (8)
        imf_val = 3. * m_star ** (-3.3) * np.exp(-(716.4 / m_star) ** 0.25)

    elif IMF_name == 'salpeter_1955':
        # Salpeter (1955)  IMF.
        # https://ui.adsabs.harvard.edu/abs/1955ApJ...121..161S/
        imf_val = m_star ** -2.35

    return imf_val


def invTrnsfSmpl(masses, IMF_name, m_low, m_high):
    """
    IMF inverse transform sampling.

    Asked here: https://stackoverflow.com/q/21100716/1391441
    """

    # Obtain normalization constant (k = \int_{m_low}^{m_up} \xi(m) dm). This
    # makes the IMF behave like a PDF.
    norm_const = quad(IMF_func, m_low, m_high, args=(IMF_name))[0]

    # IMF mass interpolation step and grid values.
    mass_step = 0.05
    mass_values = np.arange(m_low, m_high, mass_step)

    # The CDF is defined as: $F(m)= \int_{m_low}^{m} PDF(m) dm$
    # Sample the CDF
    CDF_samples = []
    for m in mass_values:
        CDF_samples.append(quad(IMF_func, m_low, m, args=(IMF_name))[0])

    # Normalize values
    CDF_samples = np.array(CDF_samples) / norm_const
    CDF_min, CDF_max = CDF_samples.min(), CDF_samples.max()

    # Inverse CDF
    inv_cdf = interp1d(CDF_samples, mass_values)

    def sampled_inv_cdf(N):
        mr = np.random.rand(N)
        mr = mr[(mr >= CDF_min) & (mr <= CDF_max)]
        return inv_cdf(mr)

    # Sample in chunks of 100 stars until the maximum defined mass is reached.
    mass_samples = []
    while np.sum(mass_samples) < masses[-1]:
        mass_samples += sampled_inv_cdf(100).tolist()
    sampled_IMF = np.array(mass_samples)

    return sampled_IMF

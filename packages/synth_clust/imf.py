
from ..core_imp import np
from scipy.integrate import quad
from scipy.interpolate import interp1d


def main(IMF_name, m_high, masses):
    """
    Returns the number of stars per interval of mass for the selected IMF.

    Parameters
    ----------
    IMF_name : str
      Name of the IMF to be used.
    m_high : float
      Maximum mass value to be sampled.
    masses: array
      Array of floats containing the masses defined within the limiting
      values and with the given mass step.

    Returns
    -------
    st_dist_mass : dict
      Dictionary that contains the N defined masses, one array per mass value.
      Each array contains a given number of stars sampled from the selected
      IMF, that (approximately) sum to the associated total mass.

    """

    print("Sampling selected IMF ({})".format(IMF_name))

    # Low mass limits are defined for each IMF to avoid numerical
    # issues when integrating.
    imfs_dict = {
        'chabrier_2001_exp': 0.01, 'chabrier_2001_log': 0.01,
        'kroupa_1993': 0.081, 'kroupa_2002': 0.011, 'salpeter_1955': 0.3}
    # IMF low mass limit.
    m_low = imfs_dict[IMF_name]

    # IMF mass interpolation step.
    mass_step = 0.05

    ##################
    # Obtain normalization constant. This is equivalent to 'k/M_total' in
    # Eq. (7) of Popescu & Hanson 2009 (138:1724-1740; PH09).
    # For m_high > 100 Mo the differences in the resulting normalization
    # constant are negligible. This is because th IMF drops very rapidly for
    # high masses.
    norm_const = 1. / quad(integral_IMF_M, m_low, m_high, args=(IMF_name))[0]

    # Obtain number of stars in each mass interval. Equivalent to the upper
    # fraction of Eq. (8) in PH09, without the total mass.
    mass_up, N_dist = [], []
    m_upper = m_low
    while m_upper < m_high:
        m_lower = m_upper
        m_upper += mass_step
        # Number of stars in the (m_lower, m_upper) interval.
        N_stars = quad(integral_IMF_N, m_lower, m_upper, args=(IMF_name))[0]
        mass_up.append(m_upper)
        N_dist.append(N_stars)

    # Normalize number of stars by constant.
    N_dist = np.asarray(N_dist) * norm_const

    # st_dist_mass = massDist(m_low, mass_up, N_dist, masses)
    ##################

    ##################
    # Obtain normalization constant (k = \int_{m_low}^{m_up} \xi(m) dm). This
    # makes the IMF behave like a PDF.
    norm_const = quad(integral_IMF_N, m_low, m_high, args=(IMF_name))[0]

    # The CDF is defined as: $F(m)= \int_{m_low}^{m} PDF(m) dm$
    # Sample the CDF
    mass_values = np.arange(m_low, m_high, mass_step)
    CDF_samples = []
    for m in mass_values:
        CDF_samples.append(quad(integral_IMF_N, m_low, m, args=(IMF_name))[0])
    CDF_samples = np.array(CDF_samples) / norm_const

    inv_cdf = interp1d(CDF_samples, mass_values)
    CDF_min, CDF_max = CDF_samples.min(), CDF_samples.max()

    def sampled_inv_cdf(N):
        mr = np.random.rand(N)
        mr = mr[(mr >= CDF_min) & (mr <= CDF_max)]
        return inv_cdf(mr)

    # sampled_IMF = sampled_inv_cdf(10000)

    # st_dist_mass2 = (sampled_IMF, np.cumsum(sampled_IMF))
    ##################

    def return_intersection(hist_1, hist_2):
        minima = np.minimum(hist_1, hist_2)
        intersection = np.true_divide(np.sum(minima), np.sum(hist_2))
        return intersection

    import matplotlib.pyplot as plt
    dd = []
    for mass in np.random.uniform(100., 1000., 1000):
        data_1 = massDist(m_low, mass_up, N_dist, [mass])
        data_1 = data_1[mass]
        sampled_IMF = sampled_inv_cdf(10000)
        data_2 = sampled_IMF[:np.searchsorted(np.cumsum(sampled_IMF), mass)]
        hist_1, _ = np.histogram(data_1, bins=100, range=[0., 150])
        hist_2, _ = np.histogram(data_2, bins=100, range=[0., 150])
        rr = return_intersection(hist_1, hist_2)
        dd.append([mass, rr])
        if rr < .7:
            print(mass, rr, len(data_1), len(data_2))
            plt.hist(data_1, 100, alpha=.5, density=True, label='old')
            plt.hist(data_2, 100, alpha=.5, density=True, label='new')
            plt.xscale('log');plt.yscale('log')
            plt.legend();plt.show()

    dd = np.array(dd).T
    plt.scatter(*dd)
    plt.show()

    return st_dist_mass


def integral_IMF_M(m_star, IMF_name):
    '''
    Returns mass values: $$
    '''
    return m_star * imfs(IMF_name, m_star)


def integral_IMF_N(m_star, IMF_name):
    '''
    Returns number of stars: $$
    '''
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


def massDist(m_low, mass_up, N_dist, masses):
    """
    Given the mass distribution from the sampled IMF, obtain for each mass
    defined: the total number of stars, and the scale and base factors
    that will distribute them properly.
    """
    st_dist_mass = {}

    for M_total in masses:
        # Normalize number of stars per interval of mass according to total
        # mass. This is equivalent to Eq. (8) in PH09.
        N_stars = N_dist * M_total
        # Distribute the N_stars for this total mass.
        m_low_i, m_up_i, N_stars_i = starsInterv(m_low, mass_up, N_stars)

        # Store mass distribution parameters.
        N_stars_total = np.sum(N_stars_i)
        base = np.repeat(m_low_i, N_stars_i)
        scale = np.repeat(m_up_i, N_stars_i) - base

        # DEPRECATED May 2019, all models need to be identical for the same
        # input parameters. This mode could be useful for the 'synth_mode'
        # mode to generate synthetic clusters (#239).
        #
        # # Either pass the values so that each synthetic cluster can generate
        # # its own mass distribution, or generate it once here.
        # if m_sample_flag is True:
        #     # Add flag to the list to avoid having to pass 'm_sample_flag'
        #     # across the code.
        #     st_dist_mass[M_total] = [base, scale, N_stars_total, True]
        # else:
        st_dist_mass[M_total] =\
            np.random.random(N_stars_total) * scale + base

    return st_dist_mass


def starsInterv(m_low, mass_up, N_stars):
    """
    Distribute the N_stars in each mass interval, so that each interval
    contains at least one star.
    """
    m_low_i, m_up_i, N_stars_i = [], [], []
    N_st_add = 0.
    for m_up, N_st in zip(*[mass_up, N_stars]):
        N_st += N_st_add
        # If the number of stars in the interval is less than 1, combine as
        # many adjacent intervals as necessary until it reaches at least one
        # star and then generate that star(s) with a random mass in the
        # m_low, m_up interval.
        if N_st < 1.:
            # Store this fraction to be added with the next one.
            N_st_add = N_st
        else:
            # The N_st stars will be randomly distributed between
            # m_low and m_up limits.
            m_low_i.append(m_low)
            m_up_i.append(m_up)
            N_stars_i.append(int(round(N_st)))
            # Reset parameters and move on to the next interval.
            N_st_add, m_low = 0., m_up

    return m_low_i, m_up_i, N_stars_i

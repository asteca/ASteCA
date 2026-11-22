import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd

# from scipy.integrate import quad
from scipy.interpolate import make_interp_spline


def main():
    """ """
    qmass_plot()
    qunif_plot()


def qmass_plot(N_m=50000):
    # fig = plt.figure(figsize=(8, 5))
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(8, 6))

    mbr = "Fisher (2005; stepped)"
    q3 = fisher_raghavan(N_m, mbr)
    y, x = np.histogram(q3, 10)
    # Correct bin edge placement
    x = [(a + x[i + 1]) / 2.0 for i, a in enumerate(x[0:-1])]
    hist = pd.Series(y, x)
    # Plot previously histogrammed data
    w = abs(hist.index[1]) - abs(hist.index[0])
    ax1.bar(
        hist.index,
        hist.values / N_m + 1,
        bottom=-1,
        width=w,
        alpha=0.5,
        align="center",
        fill=False,
        ls=":",
        lw=2,
        edgecolor="red",
        label=mbr,
    )
    # plt.hist(q3, density=False, histtype="step", ls=':', lw=2, label=mbr)

    mbr = "Fisher (2005; peaked)"
    q3 = fisher_raghavan(N_m, mbr)
    y, x = np.histogram(q3, 10)
    # Correct bin edge placement
    x = [(a + x[i + 1]) / 2.0 for i, a in enumerate(x[0:-1])]
    hist = pd.Series(y, x)
    # Plot previously histogrammed data
    w = abs(hist.index[1]) - abs(hist.index[0])
    ax1.bar(
        hist.index,
        hist.values / N_m + 1,
        bottom=-1,
        width=w,
        alpha=0.5,
        align="center",
        fill=False,
        ls="-.",
        lw=2,
        edgecolor="blue",
        label=mbr,
    )
    # plt.hist(q3, density=False, histtype="step", ls=':', lw=2, label=mbr)

    ax1.set_ylim(0, 0.26)
    ax1.set_ylabel(r"$Frequency;\;f(q)$")
    ax1.legend()

    mbr = "Raghavan (2010)"
    q3 = fisher_raghavan(N_m, mbr)
    y, x = np.histogram(q3, 20)
    # Correct bin edge placement
    x = [(a + x[i + 1]) / 2.0 for i, a in enumerate(x[0:-1])]
    hist = pd.Series(y, x)
    # Plot previously histogrammed data
    w = abs(hist.index[1]) - abs(hist.index[0])
    ax2.bar(
        hist.index,
        hist.values / N_m + 1,
        bottom=-1,
        width=w,
        alpha=0.5,
        align="center",
        fill=False,
        ls="--",
        lw=2,
        edgecolor="green",
        label=mbr,
    )
    # plt.hist(q3, 20, density=False, histtype="step", ls='--', lw=2, label=mbr)

    mbr = "Duchene & Kraus (2013)"
    # Sample IMF
    inv_cdf = invTrnsfSmpl()
    M1 = inv_cdf(np.random.rand(N_m))
    q3 = DK_gamma(M1)
    y, x = np.histogram(q3, 10)
    # Correct bin edge placement
    x = [(a + x[i + 1]) / 2.0 for i, a in enumerate(x[0:-1])]
    hist = pd.Series(y, x)
    # Plot previously histogrammed data
    w = abs(hist.index[1]) - abs(hist.index[0])
    ax2.bar(
        hist.index,
        hist.values / N_m + 1,
        bottom=-1,
        width=w,
        alpha=0.5,
        align="center",
        fill=False,
        ls="--",
        lw=2,
        edgecolor="k",
        label=mbr,
    )
    # plt.hist(
    #     q2, 50, density=False, histtype="step", color='k', lw=2, label="Duchene & Kraus (2013)")

    plt.ylim(0, 0.26)
    plt.xlabel(r"$q=m_2/m_1$")
    plt.ylabel(r"$Frequency;\;f(q)$")
    plt.legend()
    plt.savefig("qdist_mass.webp", dpi=300, bbox_inches="tight")
    plt.close()


def qunif_plot(N_m=50000):
    # q = np.linspace(0.01, 1, 10)
    # gamma = (-.5, 0, 1, 4)
    # plt.title(r"$f(q)=q^{\gamma}$")
    # for g in gamma:
    #     y = f_gamma(q, g)
    #     plt.plot(q, y, label=r"$\gamma=${:.2f}".format(g))
    plt.figure(figsize=(8, 3))

    gamma = -0.25
    print(gamma)
    q3 = f_gamma(N_m, gamma)
    y, x = np.histogram(q3, 10)
    # Correct bin edge placement
    x = [(a + x[i + 1]) / 2.0 for i, a in enumerate(x[0:-1])]
    hist = pd.Series(y, x)
    w = abs(hist.index[1]) - abs(hist.index[0])
    plt.bar(
        hist.index,
        hist.values / N_m + 1,
        bottom=-1,
        width=w,
        alpha=0.5,
        align="center",
        fill=False,
        ls="--",
        lw=2,
        edgecolor="red",
        label=r"$\gamma=${:.2f}".format(gamma),
    )

    gamma = 0.0
    print(gamma)
    q3 = f_gamma(N_m, gamma)
    y, x = np.histogram(q3, 10)
    # Correct bin edge placement
    x = [(a + x[i + 1]) / 2.0 for i, a in enumerate(x[0:-1])]
    hist = pd.Series(y, x)
    # Plot previously histogrammed data
    w = abs(hist.index[1]) - abs(hist.index[0])
    plt.bar(
        hist.index,
        hist.values / N_m + 1,
        bottom=-1,
        width=w,
        alpha=0.5,
        align="center",
        fill=False,
        ls="-.",
        lw=2,
        edgecolor="blue",
        label=r"$\gamma=${:.2f}".format(gamma),
    )

    gamma = 2
    print(gamma)
    q3 = f_gamma(N_m, gamma)
    y, x = np.histogram(q3, 10)
    # Correct bin edge placement
    x = [(a + x[i + 1]) / 2.0 for i, a in enumerate(x[0:-1])]
    hist = pd.Series(y, x)
    # Plot previously histogrammed data
    w = abs(hist.index[1]) - abs(hist.index[0])
    plt.bar(
        hist.index,
        hist.values / N_m + 1,
        bottom=-1,
        width=w,
        alpha=0.5,
        align="center",
        fill=False,
        ls="--",
        lw=2,
        edgecolor="green",
        label=r"$\gamma=${:.2f}".format(gamma),
    )

    plt.ylim(0, 0.29)
    plt.xlabel(r"$q=m_2/m_1$")
    plt.ylabel(r"$Frequency;\;f(q)$")
    plt.legend()
    plt.savefig("qdist_unif.webp", dpi=300, bbox_inches="tight")
    plt.close()

    # plt.subplot(222)
    # plt.hist(M1, 100, density=True, histtype="step", color='k', lw=2)
    # plt.xscale("log")
    # plt.xlabel("m1")
    # plt.ylabel("N")

    # plt.subplot(223)
    # for g in gamma:
    #     q1 = np.random.power(g + 1, 1000)
    #     plt.hist(q1, density=True, histtype="step", lw=2, label="g={:.2f}".format(g))
    # plt.xlabel("q=m2/m1")
    # plt.ylabel("f(q)")
    # plt.legend()


def f_gamma(N, gamma):
    def fQ(xk, pk):
        """
        Discrete function
        """
        pk /= pk.sum()
        fq = stats.rv_discrete(a=0.0, b=1.0, values=(xk, pk))
        return fq

    xk = np.linspace(0.01, 1.0, 10)
    pk = xk**gamma
    # f_q = q**(gamma)
    # f_q /= max(f_q)
    fq = fQ(xk, pk)
    # 'ppf' is the inverse CDF
    mass_ratios = fq.ppf(np.random.uniform(0.0, 1.0, N))

    return mass_ratios


def DK_gamma(M1):
    """
    Use 'gamma + 1' in the power-law distribution below because in D&K this
    distribution is defined as f(q)~q^gamma, while numpy's distribution is
    defined as a*x^(a-1).
    """
    msk1, gamma1 = M1 <= 0.1, 4.2
    msk2, gamma2 = (M1 > 0.1) & (M1 <= 0.6), 0.4
    msk3, gamma3 = (M1 > 0.6) & (M1 <= 1.4), 0.3
    msk4, gamma4 = (M1 > 1.4) & (M1 <= 6.5), -0.5
    msk5, gamma5 = (M1 > 6.5) & (M1 <= 16), 0.0  # <- Not sure. Use uniform
    msk6, gamma6 = M1 > 16, 0.0  # <- Not sure. Use uniform

    mass_ratios = np.zeros(M1.size)
    for msk, gamma in (
        (msk1, gamma1),
        (msk2, gamma2),
        (msk3, gamma3),
        (msk4, gamma4),
        (msk5, gamma5),
        (msk6, gamma6),
    ):
        q = np.random.power(gamma + 1, msk.sum())
        mass_ratios[msk] = q

    return mass_ratios


def invTrnsfSmpl(m_low=0.08, m_high=150):
    """
    IMF inverse transform sampling.
    """
    # The lower mass region needs to be sampled more accurately
    mass_values = list(np.linspace(m_low, 0.5, 100))
    mass_values += list(np.linspace(0.501, 2, 75))
    mass_values += list(np.linspace(2.01, 10, 50))
    mass_values += list(np.linspace(10.01, m_high, 25))

    def IMF_func(m_star):
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

    CDF_samples = []
    IMF_old, m_old, area_CDF = IMF_func(m_low), m_low, 0.0
    for m in mass_values[1:]:
        # Approximate integral with rectangular area, and add to previous total area
        IMF_new = IMF_func(m)
        area_CDF += 0.5 * (IMF_new + IMF_old) * (m - m_old)
        CDF_samples.append(area_CDF)
        IMF_old, m_old = IMF_new, m
    CDF_samples = np.array(CDF_samples) / max(CDF_samples)

    # k=1 is important otherwise the interpolator can become unstable for values
    # close to 1 and return huge masses and even negative ones
    inv_cdf = make_interp_spline(CDF_samples, mass_values[1:], k=1)

    return inv_cdf


def fisher_raghavan(N, mbr):
    def fQ(xk, pk):
        """
        Discrete function
        """
        pk /= pk.sum()
        fq = stats.rv_discrete(a=0.0, b=1.0, values=(xk, pk))
        return fq

    # Fisher's distribution
    if mbr == "Fisher (2005; stepped)":
        # Fisher, Schröder & Smith (2005), 10.1111/j.1365-2966.2005.09193.x;
        # Table 3, stepped
        xk = np.linspace(0.0, 1.0, 10)
        pk = np.array([29.0, 29.0, 30.0, 32.0, 31.0, 32.0, 36.0, 45.0, 27.0, 76.0])

    elif mbr == "Fisher (2005; peaked)":
        # Fisher, Schröder & Smith (2005), 10.1111/j.1365-2966.2005.09193.x;
        # Table 3, peaked
        xk = np.linspace(0.0, 1.0, 10)
        pk = np.array([27.0, 30.0, 34.0, 33.0, 29.0, 26.0, 27.0, 33.0, 41.0, 89.0])

    elif mbr == "Raghavan (2010)":
        # Raghavan et al. (2010), 10.1088/0067-0049/190/1/1; Fig 16 (left)
        xk = np.linspace(0.0, 1.0, 20)
        pk = np.array(
            [
                0.53,
                2.61,
                0.53,
                4.67,
                7.81,
                3.64,
                9.89,
                5.71,
                4.69,
                5.73,
                4.67,
                6.76,
                5.75,
                5.73,
                2.61,
                5.71,
                4.72,
                5.71,
                3.64,
                12.99,
            ]
        )

    fq = fQ(xk, pk)
    # 'ppf' is the inverse CDF
    mass_ratios = fq.ppf(np.random.uniform(0.0, 1.0, N))

    return mass_ratios


if __name__ == "__main__":
    main()

import matplotlib.pyplot as plt
import numpy as np


def main():
    fig = plt.figure(figsize=(8, 4))

    x_vals = np.array(
        list(np.arange(0.02, 1, 0.01)) + [1] + [1, 2, 3, 5, 10, 15, 20, 50, 100]
    )
    for imf_f in ("Salpeter (1955)", "Kroupa (2001)", "Chabrier et al. (2014)"):
        print(imf_f)
        y = []
        for x in x_vals:
            y.append(get_imf(imf_f, x))
        # breakpoint()
        idx = np.argmin(abs(x_vals - 1))
        y = np.array(y) / y[idx]
        # plt.plot(np.log10(x_vals), y, label=imf_f)
        plt.plot(x_vals, y, label=imf_f)
        plt.legend()

    plt.xlim(0.018, 20)
    plt.ylim(0.0009, 700)
    plt.xscale("log")
    plt.xlabel(r"$m\,[M_{\odot}]$")
    plt.ylabel(r"$\xi\,(m)$")
    plt.yscale("log")
    plt.savefig("IMFs.png", dpi=300, bbox_inches="tight")


def get_imf(IMF_name, m_star):
    """
    Define any number of IMFs.

    The package https://github.com/keflavich/imf has some more (I think,
    24-09-2019).
    """
    if IMF_name == "Salpeter (1955)":
        # Salpeter (1955)  IMF.
        # https://ui.adsabs.harvard.edu/abs/1955ApJ...121..161S/
        imf_val = m_star**-2.35

    elif IMF_name == "Kroupa (2001)":
        # Kroupa (2001), Mon. Not. R. Astron. Soc. 322, 231-246 (2001); Eq. 2
        # https://ui.adsabs.harvard.edu/abs/2001MNRAS.322..231K/abstract
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

    elif IMF_name == "Chabrier et al. (2014)":
        # Chabrier et al. (2014)
        # https://ui.adsabs.harvard.edu/abs/2014ApJ...796...75C/abstract ; Eq (34)
        nc, mc = 11, 0.18
        m0 = nc * mc
        Ah, x = 0.649, 1.35
        Al = Ah * nc ** (x / 2)
        c = 0.434294 / m_star  # np.log10(e)/m
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


if __name__ == "__main__":
    main()

import matplotlib.pyplot as plt
import numpy as np


def main():
    plt.figure(figsize=(8, 4))

    x_vals = np.array(
        list(np.arange(0.02, 1, 0.01)) + [1] + [1, 2, 3, 5, 10, 15, 20, 50, 100]
    )
    for imf_f in ("Salpeter (1955)", "Kroupa (2001)", "Chabrier et al. (2014)"):
        print(imf_f)
        y = get_imf(imf_f, x_vals)
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
    # plt.show()
    plt.savefig("IMFs.webp", dpi=300, bbox_inches="tight")


def get_imf(IMF_name: str, m_star_array: np.ndarray) -> np.ndarray:
    """
    Define any number of IMFs.

    The package https://github.com/keflavich/imf has some more (I think).
    """
    IMF_name = (
        IMF_name.lower()
        .replace(" et al.", "")
        .replace(" ", "_")
        .replace("(", "")
        .replace(")", "")
    )
    print(IMF_name)

    if IMF_name == "salpeter_1955":
        # Salpeter (1955)  IMF.
        # https://ui.adsabs.harvard.edu/abs/1955ApJ...121..161S/
        imf_vals = m_star_array**-2.35

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

        # Conditions for m_star_array <= m0 and m_star_array > m0
        lower_msk = (m0 <= m_star_array) & (m_star_array <= m1)
        middle_msk = (m1 < m_star_array) & (m_star_array <= m2)
        upper_msk = m_star_array > m2

        # Initialize output array
        imf_vals = np.zeros_like(m_star_array)
        imf_vals[lower_msk] = factor[0] * (m_star_array[lower_msk] ** alpha[0])
        imf_vals[middle_msk] = factor[1] * (m_star_array[middle_msk] ** alpha[1])
        imf_vals[upper_msk] = factor[2] * (m_star_array[upper_msk] ** alpha[2])

    elif IMF_name == "chabrier_2014":
        # Chabrier et al. (2014)
        # https://ui.adsabs.harvard.edu/abs/2014ApJ...796...75C/ ; Eq (34)
        nc, mc = 11, 0.18
        m0 = nc * mc
        Ah, x = 0.649, 1.35
        Al = Ah * nc ** (x / 2)
        sigma_2 = np.log10(nc) / (x * np.log(10))
        c_array = 0.434294 / m_star_array  # np.log10(e)/m ; This is the transformation
        # from dN/log(m) --> dN/dm

        # Conditions for m_star_array <= m0 and m_star_array > m0
        lower_msk = m_star_array <= m0
        upper_msk = m_star_array > m0

        # Initialize output array
        imf_vals = np.zeros_like(m_star_array)
        # Compute imf_val for m_star <= m0
        imf_vals[lower_msk] = (
            c_array[lower_msk]
            * Al
            * m0 ** (-x)
            * np.exp(
                -((np.log10(m_star_array[lower_msk]) - np.log10(mc)) ** 2)
                / (2 * sigma_2)
            )
        )
        # Compute imf_val for m_star > m0
        imf_vals[upper_msk] = c_array[upper_msk] * Ah * m_star_array[upper_msk] ** (-x)

    return imf_vals


if __name__ == "__main__":
    main()

import matplotlib.pyplot as plt
import numpy as np


def main():
    """

    The coefficients should not be extrapolated outside the extinction and
    temperature ranges used, the latter corresponding to
    -0.06< (GBP − GRP )0 < 2.5 mag
    """
    eff_wave = np.linspace(1000, 10000, 1000)  # Angstrom
    eff_micron = eff_wave / 10000

    fig, axes = plt.subplots(1, 2, figsize=(10, 5))

    axes[0].set_title("CCMO")

    for Rv in (2.5, 3.1, 5):
        y = []
        for ew in eff_wave:
            a, b = ccmo_model(10000.0 / ew)
            y.append(a + b / Rv)
        axes[0].plot(1 / eff_micron, y, label=f"Rv={Rv:.1f}")

    axes[0].set_ylabel(r"$A(\lambda)/A_V$")
    axes[0].set_xlabel(r"$\lambda^{-1} \; [\mu m^{-1}$]")
    axes[0].legend()

    # If this model is used the first color is always expected to be BP-RP
    BP_RP = np.linspace(-0.05, 2.5, 100)
    for Av in (0.5, 1, 2):
        ec_mag, ec_col1 = gaiadr3_ext_law(BP_RP, Av)
        # Ex1 = ec_col1 * Av
        # BP_RP_Av = BP_RP + Ex1
        plt.plot(BP_RP, ec_mag, label=f"Av={Av:.1f}")

    axes[1].set_title("Gaia (E)DR3")
    axes[1].set_xlabel("BP-RP")
    axes[1].set_ylabel(r"$A_G/A_V$")
    axes[1].legend()

    fig.tight_layout()
    plt.savefig("../ext_law.webp", dpi=300, bbox_inches="tight")


def ccmo_model(mw: float) -> tuple[float, float]:
    """Cardelli, Clayton, and Mathis (1989 ApJ. 345, 245) model for extinction
    coefficients with updated coefficients for near-UV from O'Donnell
    (1994, ApJ, 422, 158).

    ccm_coef = a + b / Rv

    Implementation taken from:

    http://idlastro.gsfc.nasa.gov/ftp/pro/astro/ccm_unred.pro

    There appears to be an error in the Far-UV range in the original IDL
    routine where the maximum inverse wavelength is 11 and it should be 10
    according to Cardelli et al. 1989 (pag 251, Eq (5,a,b)).

    :param mw: Wavelength in inverse microns.
    :type mw: float

    :raises ValueError: If the effective wavelength is beyond the CCM model limit

    :returns: Extinction coefficients a and b.
    :rtype: tuple[float, float]
    """

    if 0.3 <= mw < 1.1:
        # Infrared.
        a, b = 0.574 * (mw**1.61), -0.527 * (mw**1.61)

    elif 1.1 <= mw < 3.3:
        # Optical/NIR.
        # Original coefficients from CCM89
        # c1 = [1., 0.17699, -0.50447, -0.02427, 0.72085, 0.01979, -0.77530,
        #       0.32999]
        # c2 = [0., 1.41338, 2.28305, 1.07233, -5.38434, -0.62251, 5.30260,
        #       -2.09002]
        # New coefficients from O'Donnell (1994)
        c1 = [1.0, 0.104, -0.609, 0.701, 1.137, -1.718, -0.827, 1.647, -0.505]
        c2 = [0.0, 1.952, 2.908, -3.989, -7.985, 11.102, 5.491, -10.805, 3.347]
        y = mw - 1.82
        # Reverse because polyval starts from the highest degree.
        c1.reverse()
        c2.reverse()
        a, b = np.polyval(c1, y), np.polyval(c2, y)

    elif 3.3 <= mw < 8.0:
        # Mid-UV
        F_a, F_b = 0.0, 0.0
        if mw >= 5.9:
            y = mw - 5.9
            F_a = -0.04473 * y**2 - 0.009779 * y**3
            F_b = 0.2130 * y**2 + 0.1207 * y**3
        a = 1.752 - 0.316 * mw - (0.104 / ((mw - 4.67) ** 2 + 0.341)) + F_a
        b = -3.090 + 1.825 * mw + (1.206 / ((mw - 4.62) ** 2 + 0.263)) + F_b

    elif 8.0 <= mw <= 10.0:
        # Far-UV
        c1 = [-1.073, -0.628, 0.137, -0.070]
        c2 = [13.670, 4.257, -0.420, 0.374]
        y = mw - 8.0
        c1.reverse()
        c2.reverse()
        a, b = np.polyval(c1, y), np.polyval(c2, y)
    else:
        raise ValueError(
            "The effective wavelength is {} [1/micron], beyond "
            "the CCM model limit (10 [1/micron]).".format(mw)
        )

    return float(a), float(b)


def gaiadr3_ext_law(
    X_: np.ndarray, Av_dr: float | np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """
    The 'coeffs' values are the main sequence values taken from:
    https://www.cosmos.esa.int/web/gaia/edr3-extinction-law

    The order of the coefficients is:
    Intercept   X   X2  X3  A   A2  A3  XA  AX2 XA2

    :param X_: Array of BP-RP colors.
    :type X_: np.ndarray
    :param Av_dr: Total absorption (eventually containing differential reddening).
    :type Av_dr: float | np.ndarray

    :returns: Extinction coefficients for G and BP-RP.
    :rtype: tuple[np.ndarray, np.ndarray]
    """
    coeffs = {
        "G": (
            0.995969721536602,
            -0.159726460302015,
            0.0122380738156057,
            0.00090726555099859,
            -0.0377160263914123,
            0.00151347495244888,
            -2.52364537395142e-05,
            0.0114522658102451,
            -0.000936914989014318,
            -0.000260296774134201,
        ),
        "BP": (
            1.15363197483424,
            -0.0814012991657388,
            -0.036013023976704,
            0.0192143585568966,
            -0.022397548243016,
            0.000840562680547171,
            -1.31018008013549e-05,
            0.00660124080271006,
            -0.000882247501989453,
            -0.000111215755291684,
        ),
        "RP": (
            0.66320787941067,
            -0.0179847164933981,
            0.000493769449961458,
            -0.00267994405695751,
            -0.00651422146709376,
            3.30179903473159e-05,
            1.57894227641527e-06,
            -7.9800898337247e-05,
            0.000255679812110045,
            1.10476584967393e-05,
        ),
    }

    X_2 = X_ * X_
    X_3 = X_2 * X_
    Av_2 = Av_dr * Av_dr
    Av_3 = Av_2 * Av_dr

    def ext_coeff(k):
        """
        https://www.cosmos.esa.int/web/gaia/edr3-extinction-law
        """
        c = coeffs[k]
        ay = (
            c[0]
            + c[1] * X_
            + c[2] * X_2
            + c[3] * X_3
            + c[4] * Av_dr
            + c[5] * Av_2
            + c[6] * Av_3
            + c[7] * X_ * Av_dr
            + c[8] * X_2 * Av_dr
            + c[9] * X_ * Av_2  # This index not a mistake
        )
        return ay

    ec_G = ext_coeff("G")
    ec_BPRP = ext_coeff("BP") - ext_coeff("RP")

    return ec_G, ec_BPRP


if __name__ == "__main__":
    main()

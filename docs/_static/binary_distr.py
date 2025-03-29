import matplotlib.pyplot as plt
import numpy as np


def main():
    """
    Find the best function to fit Duchene & Kraus (2013) Multiplicity
    Frequency distribution. This is used by ASteCA to assign mass-dependent
    binary probabilities.
    """

    # # Multiplicity frequency; Duchene & Kraus (2013) Fig 1
    # mass_DK = np.array([0.08, .29, .96, 2.4, 7.65, 28.5, 151])
    # multfreq_DK = (.22, .26, .46, .52, .64, .82, 1)
    # plt.scatter(mass_DK, multfreq_DK, marker='s', c='orange', label="D&K (2013)")

    # # Multiplicity frequency; Arenou (2011)
    # mass_AR = np.array([0.0297, 0.0691, 0.198, 0.392, 0.994, 2.979, 19.841])
    # multfreq_AR = np.array([0.088, 0.126, 0.149, 0.279, 0.558, 0.801, 0.851])
    # plt.scatter(mass_AR, multfreq_AR, marker='^', c='b', label="Arenou (2011)")

    # Multiplicity frequency; Offner et al. (2022) Fig 1 (left)
    mass_OF = np.array(
        [
            0.033,
            0.062,
            0.104,
            0.210,
            0.424,
            0.861,
            1.040,
            1.133,
            1.935,
            3.803,
            6.289,
            11.519,
            29.118,
        ]
    )
    multfreq_OF = np.array(
        [
            0.077,
            0.149,
            0.187,
            0.229,
            0.299,
            0.419,
            0.469,
            0.499,
            0.679,
            0.806,
            0.888,
            0.928,
            0.956,
        ]
    )
    x = np.linspace(0.029, 32, 100)
    y = 0.09 + 0.94 * (x / (1.4 + x))

    plt.figure(figsize=(8, 4))
    plt.scatter(
        mass_OF,
        multfreq_OF,
        marker="+",
        s=150,
        lw=3,
        c="b",
        label="Offner et al. (2023)",
        zorder=4,
    )
    plt.plot(x, y, "--k", label=r"$\alpha=0.09;\;\beta=0.94$")

    plt.ylim(0, 1)
    plt.xscale("log")
    plt.xlabel(r"$m_1\,[M_{\odot}]$")
    plt.ylabel(r"$P_b\,(m_1)$")
    plt.legend()
    plt.savefig("binar_distr.webp", dpi=300, bbox_inches="tight")
    breakpoint()

    # x = np.log(mass_OF)

    # y = 0.112 + 0.563*np.arctan(0.765*multfreq_OF)
    # plt.scatter(x, y, marker='v', c='k')
    # slope, intercept, r_value, _, _ = scipy.stats.linregress(x, y)
    # poly1d_fn = np.poly1d([slope, intercept])
    # y_fit = poly1d_fn(x)
    # MSE = np.square(y_fit - y).sum()
    # plt.plot(x, y_fit, '--k', label="R^2={:.3f}, MSE={:.3f}".format(r_value**2, MSE))

    # y = 0.091 + 0.939*(multfreq_OF/(1.4+multfreq_OF))
    # plt.scatter(x, y, marker='v', c='g')
    # slope, intercept, r_value, _, _ = scipy.stats.linregress(x, y)
    # poly1d_fn = np.poly1d([slope, intercept])
    # y_fit = poly1d_fn(x)
    # MSE = np.square(y_fit - y).sum()
    # plt.plot(x, y_fit, '--g', label="R^2={:.3f}, MSE={:.3f}".format(r_value**2, MSE))

    # y = 0.673 + .153*np.log(.367*multfreq_OF)
    # plt.scatter(x, y, marker='v', c='b')
    # slope, intercept, r_value, _, _ = scipy.stats.linregress(x, y)
    # poly1d_fn = np.poly1d([slope, intercept])
    # y_fit = poly1d_fn(x)
    # MSE = np.square(y_fit - y).sum()
    # plt.plot(x, y_fit, '--b', label="R^2={:.3f}, MSE={:.3f}".format(r_value**2, MSE))

    # plt.legend()
    # plt.show()
    # breakpoint()

    # Find the optimal (alpha, beta) fit
    mass_fit, multfreq_fit = mass_OF, multfreq_OF

    beta_opt, alpha_opt, c_opt, delta_old = np.nan, np.nan, np.nan, np.inf
    N = 50
    c = 1.4
    for alpha in np.linspace(0.0, 0.25, N):
        for beta in np.linspace(0.75, 1.25, N):
            # for c in np.linspace(1.3, 1.5, N):
            y = betaFunc(mass_fit, alpha, beta, c)
            delta = np.square(multfreq_fit - y).sum()
            if delta < delta_old:
                beta_opt, alpha_opt, c_opt = beta, alpha, c
                delta_old = delta

    plt.title("MSS: {:.4f}".format(delta_old))
    x = np.linspace(0.001, 40, 100)
    y = betaFunc(x, alpha_opt, beta_opt, c_opt)
    plt.plot(
        x,
        y,
        ls=":",
        label="alpha={:.3f}, beta={:.3f}, c={:.3f}".format(alpha_opt, beta_opt, c_opt),
        c="magenta",
    )

    plt.xlim(0.01, 40)
    plt.ylim(0, 1)
    # plt.xscale('log')
    plt.legend()
    plt.show()


def betaFunc(x, alpha, beta, c, mcut=0.5):
    # return alpha + beta*np.tanh(x)        # 0.141
    # return alpha - np.exp(-beta*x)        # 0.071
    # return alpha*np.tanh(x*beta)          # 0.068
    # return alpha + 1/(1+beta*np.exp(-x))  # 0.063
    # return alpha + beta*np.log(x)         # 0.056
    # return 1/(alpha+beta*np.exp(-x))      # 0.043
    # return alpha + beta*(1 / (1+1/x))     # 0.017
    # return alpha + beta*np.arctan(x)      # 0.016

    return alpha + beta * (x / (c + x))  # 0.0088
    # return alpha + beta*np.arctan(c*x)    # 0.0061
    # return alpha + beta*np.log(c*x)         # 0.056


if __name__ == "__main__":
    main()

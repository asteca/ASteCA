"""Figures for the 'foreground' and 'internal' differential reddening models,
A = Av + f*DR (run from this folder). The fits to observed clusters figure
(dr_models_fits.webp) was made separately, from pocoMC fits of UCC members."""

import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, Polygon
from scipy import stats

import asteca

OUT = ".."

C_INT, C_FG = "#d95f02", "#1b6fb5"

# Truncated normal(0.5, 1/6) at +/-3 sigma: f in [0, 1]
F_FG = stats.truncnorm(-3, 3, loc=0.5, scale=1 / 6)


def f_internal_pdf(f):
    """Stars and dust uniformly mixed in a sphere: p(f) = 3/2 (1 - f^2)"""
    return np.where((f >= 0) & (f <= 1), 1.5 * (1 - f**2), 0.0)


def f_moments():
    m_int, s_int = 3 / 8, np.sqrt(0.2 - (3 / 8) ** 2)
    return m_int, s_int, F_FG.mean(), F_FG.std()


def matched_internal(Av_fg, DR_fg):
    """Internal-mode (Av, DR) with the same mean and spread of A as the given
    foreground-mode values"""
    m_int, s_int, m_fg, s_fg = f_moments()
    DR = DR_fg * s_fg / s_int
    return Av_fg + m_fg * DR_fg - m_int * DR, DR


def plot_geometry():
    rng = np.random.default_rng(7)
    fig, axes = plt.subplots(2, 1, figsize=(8, 7.2))

    # Stars uniformly distributed in the (projected) sphere
    v = rng.normal(size=(140, 3))
    v /= np.linalg.norm(v, axis=1)[:, None]
    v *= rng.uniform(size=140)[:, None] ** (1 / 3)
    xs, ys = v[:, 2], v[:, 1]

    for ax, mode in zip(axes, ("internal", "foreground")):
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_xlim(-5.0, 1.6)
        ax.set_ylim(-1.55, 1.55)
        # Observer and line of sight
        ax.annotate("", xy=(-4.95, -1.35), xytext=(-3.75, -1.35),
                    arrowprops=dict(arrowstyle="->", lw=1.2, color="0.3"))
        ax.text(-4.95, -1.25, "to observer", fontsize=9, color="0.3", va="bottom")

        if mode == "internal":
            ax.set_title("'internal' mode: dust mixed with the stars", fontsize=12)
            # Uniform foreground screen
            ax.add_patch(Polygon([(-3.6, -1.2), (-3.3, -1.2), (-3.3, 1.2), (-3.6, 1.2)],
                                 fc="0.55", ec="none", alpha=0.6))
            ax.text(-3.45, 1.27, r"uniform screen: $A_V$", ha="center", fontsize=10)
            ax.add_patch(Circle((0, 0), 1, fc="#c9a26b", ec="#8a6a3b", alpha=0.45))
            ax.scatter(xs, ys, s=6, c="k", zorder=3)
            # Full column through the centre
            ax.plot([-1, 1], [0, 0], ls="--", c="#8a6a3b", lw=1.3)
            ax.text(-0.3, -1.08, r"$DR$: dust column through the centre (dashed)"
                    "\n" r"$\ell_i$: dust in front of star $i$ (orange)",
                    ha="center", va="top", fontsize=9, color="#5b4422")
            # One star and the dust in front of it
            x0, y0 = 0.35, 0.55
            xe = -np.sqrt(1 - y0**2)
            ax.plot([xe, x0], [y0, y0], c=C_INT, lw=3, zorder=4)
            ax.scatter([x0], [y0], s=60, c=C_INT, ec="k", zorder=5)
            ax.annotate(r"$A_i = A_V + f_i\,DR$" "\n" r"$f_i = \ell_i\,/\,2R$",
                        xy=((xe + x0) / 2, y0 + 0.02), xytext=(-3.1, 0.35),
                        fontsize=10, va="center",
                        arrowprops=dict(arrowstyle="->", color=C_INT))
            ax.text(1.05, 0.75, "absorption grows with\nthe depth of each star",
                    fontsize=9, ha="left", color="0.25")
        else:
            ax.set_title("'foreground' mode: patchy screen, dust-free cluster",
                         fontsize=12)
            # Patchy screen: thickness varies across the face of the cluster
            y = np.linspace(-1.2, 1.2, 200)
            t = (np.sin(3.1 * y + 0.4) + 0.6 * np.sin(7.3 * y + 1.3)
                 + 0.3 * np.sin(13 * y))
            t = (t - t.min()) / (t.max() - t.min())
            x_front, w_min, w_rng = -3.8, 0.12, 0.7
            xb = x_front + w_min + w_rng * t
            ax.add_patch(Polygon(np.c_[np.r_[np.full(y.size, x_front), xb[::-1]],
                                       np.r_[y, y[::-1]]],
                                 fc="0.55", ec="none", alpha=0.6))
            ax.plot([x_front + w_min] * 2, [-1.2, 1.2], ls=":", c="0.2", lw=1)
            ax.text(x_front - 0.05, 1.27, r"$A_V$ (minimum)  →  $A_V + DR$",
                    ha="left", fontsize=10)
            ax.add_patch(Circle((0, 0), 1, fc="white", ec="0.4"))
            ax.scatter(xs, ys, s=6, c="k", zorder=3)
            # Two stars at the same height and different depth
            y0 = ys[np.argmin(np.abs(ys - 0.55))]
            j = np.argmin(np.abs(y - y0))
            for x0 in (-0.55, 0.55):
                ax.scatter([x0], [y0], s=60, c=C_FG, ec="k", zorder=5)
            ax.plot([xb[j], 0.55], [y0, y0], c=C_FG, lw=1.2, ls="--")
            ax.plot([x_front, xb[j]], [y0, y0], c=C_FG, lw=3)
            ax.annotate(r"$A_i = A_V + f_i\,DR$" "\n"
                        r"$f_i$: local screen thickness",
                        xy=(x_front + 0.2, y0), xytext=(-3.05, -0.8), fontsize=9.5,
                        va="center", arrowprops=dict(arrowstyle="->", color=C_FG))
            ax.text(1.05, 0.75, "same absorption for stars\non the same line of sight",
                    fontsize=9, ha="left", color="0.25")

    fig.tight_layout()
    fig.savefig(f"{OUT}/dr_models_geometry.webp", dpi=200, bbox_inches="tight")


def plot_distr():
    m_int, s_int, m_fg, s_fg = f_moments()
    f = np.linspace(-0.02, 1.02, 600)
    rng = np.random.default_rng(3)
    u = rng.uniform(size=200_000)
    f_int_s = 2 * np.sin(np.arcsin(u) / 3)
    f_fg_s = F_FG.rvs(200_000, random_state=rng)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))

    ax = axes[0]
    ax.hist(f_int_s, 60, range=(0, 1), density=True, color=C_INT, alpha=0.25)
    ax.hist(f_fg_s, 60, range=(0, 1), density=True, color=C_FG, alpha=0.25)
    ax.plot(f, f_internal_pdf(f), c=C_INT, lw=2,
            label=rf"internal: $\frac{{3}}{{2}}(1-f^2)$, mean {m_int:.3f}, sd {s_int:.3f}")
    ax.plot(f, F_FG.pdf(f), c=C_FG, lw=2,
            label=rf"foreground: $\mathcal{{N}}(0.5, 1/6)$, mean {m_fg:.3f}, sd {s_fg:.3f}")
    ax.set_xlabel(r"$f$")
    ax.set_ylabel("density")
    ax.set_title(r"Distribution of $f$ (both bounded in $[0, 1]$)")
    ax.legend(fontsize=8.5, loc="upper right")
    ax.set_ylim(0, 3.3)

    def a_pdf(a, Av, DR, mode):
        x = (a - Av) / DR
        p = f_internal_pdf(x) if mode == "internal" else F_FG.pdf(x)
        return p / DR

    a = np.linspace(0, 2.0, 800)
    Av, DR = 0.5, 1.0
    ax = axes[1]
    for mode, c, m, s in (("internal", C_INT, m_int, s_int),
                          ("foreground", C_FG, m_fg, s_fg)):
        ax.plot(a, a_pdf(a, Av, DR, mode), c=c, lw=2,
                label=rf"{mode}: $\langle A\rangle$={Av + m * DR:.2f}, "
                      rf"sd={s * DR:.2f}")
        ax.axvline(Av + m * DR, c=c, ls=":", lw=1.2)
    ax.axvspan(Av, Av + DR, color="0.9", zorder=0)
    ax.set_title(rf"Same parameters: $A_V$={Av}, $DR$={DR}")
    ax.set_xlabel(r"$A$")
    ax.legend(fontsize=8.5, loc="upper right")
    ax.set_ylim(0, 3.3)

    Av_i, DR_i = matched_internal(Av, DR)
    ax = axes[2]
    ax.plot(a, a_pdf(a, Av_i, DR_i, "internal"), c=C_INT, lw=2,
            label=rf"internal: $A_V$={Av_i:.2f}, $DR$={DR_i:.2f}")
    ax.plot(a, a_pdf(a, Av, DR, "foreground"), c=C_FG, lw=2,
            label=rf"foreground: $A_V$={Av}, $DR$={DR}")
    ax.axvline(Av + m_fg * DR, c="0.3", ls=":", lw=1.2)
    ax.set_title(rf"Same mean ({Av + m_fg * DR:.2f}) and sd ({s_fg * DR:.2f}) of $A$")
    ax.set_xlabel(r"$A$")
    ax.legend(fontsize=8.5, loc="upper right")
    ax.set_ylim(0, 3.3)

    for ax in axes:
        ax.grid(alpha=0.25, lw=0.5)
    fig.tight_layout()
    fig.savefig(f"{OUT}/dr_models_distr.webp", dpi=200, bbox_inches="tight")


def plot_cmd():
    isochs = asteca.Isochrones(
        model="parsec", isochs_path="../isochrones/parsec", mag="Gmag",
        color=("G_BPmag", "G_RPmag"), verbose=0,
    )
    p0 = {"loga": 9.0, "met": 0.0152, "alpha": 0.0, "beta": 0.0, "dm": 9.0}
    Av, DR = 0.5, 1.0
    Av_i, DR_i = matched_internal(Av, DR)
    m_int, _, m_fg, _ = f_moments()
    panels = (
        ("internal", "internal", Av, DR, C_INT, m_int),
        ("foreground", "foreground", Av, DR, C_FG, m_fg),
        ("internal (same mean and sd as foreground)", "internal", Av_i, DR_i, C_INT,
         m_int),
    )

    fig, axes = plt.subplots(1, 3, figsize=(15, 5.6), sharey=True)
    for ax, (title, method, av, dr, c, m) in zip(axes, panels):
        sc = asteca.Synthetic(
            isochs, ext_law="GAIADR3", DR_distribution=method, seed=42, verbose=0
        )
        mag, col = sc.generate(dict(p0, Av=av, DR=dr), N_stars=5000)[:2]
        ax.scatter(col, mag, s=6, c=c, alpha=0.45, lw=0)
        for a, ls, lab in ((av, "-", r"$A_V$"), (av + m * dr, "--", r"$\langle A\rangle$"),
                           (av + dr, ":", r"$A_V+DR$")):
            iso = sc.get_isochrone(dict(p0, Av=a, DR=0.0))
            ax.plot(iso[1], iso[0], c="k", ls=ls, lw=1.2, label=f"{lab} = {a:.2f}")
        ax.set_title(f"{title}\n$A_V$={av:.2f}, $DR$={dr:.2f}", fontsize=11)
        ax.set_xlabel("BP-RP")
        ax.legend(fontsize=9, loc="upper right", title="isochrone at", title_fontsize=9)
        ax.grid(alpha=0.25, lw=0.5)
        ax.set_xlim(0.4, 2.3)
    axes[0].set_ylim(18.6, 9.4)
    axes[0].set_ylabel("G")
    fig.tight_layout()
    fig.savefig(f"{OUT}/dr_models_cmd.webp", dpi=200, bbox_inches="tight")


if __name__ == "__main__":
    todo = sys.argv[1:] or ["geometry", "distr", "cmd"]
    for name in todo:
        {"geometry": plot_geometry, "distr": plot_distr, "cmd": plot_cmd}[name]()
        print(name, "done")

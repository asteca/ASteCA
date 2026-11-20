import matplotlib.pyplot as plt
import numpy as np

import asteca
from asteca.modules.synthetic import stellar_stats as ssf

# Color shift applied to the turn-off point, same value used by
# ssf.get_bss_probabilities()
COL_OFFSET = -0.05


def main():
    # Isochrones parameters
    isochs = asteca.Isochrones(
        model="parsec",
        isochs_path="../isochrones/parsec/gaia",
        mag="Gmag",
        color=("G_BPmag", "G_RPmag"),
        verbose=3,
    )

    # Synthetic clusters parameters
    synthcl = asteca.Synthetic(
        isochs,
        ext_law="GAIADR3",
        seed=42519,
        verbose=3,
    )

    model = {
        "loga": 8.0,
        "dm": 8.1,
        "Av": 0.42,
    }

    # Isochrone cut at the maximum mass of the synthetic cluster
    isoch_mag, isoch_col = synthcl.get_isochrone(model)

    # Turn-off point, and the point that defines the BSS region
    # (the isochrone points are ordered by mass: their index is used as the mass)
    track = np.array([isoch_mag, isoch_col, np.arange(len(isoch_mag))])
    to = ssf.turn_off_points([track], [np.inf], 0)
    to_col, to_mag = to["to_col"][0], to["to_mag"][0]
    bss_col = to_col + COL_OFFSET

    mplot(isoch_mag, isoch_col, to_col, to_mag, bss_col)


def mplot(isoch_mag, isoch_col, to_col, to_mag, bss_col):
    _fig, ax = plt.subplots(figsize=(4, 4))

    # Isochrone
    ax.plot(isoch_col, isoch_mag, color="black", zorder=3)

    ax.scatter(to_col, to_mag, color="red", s=10, label="TO point", zorder=4)
    ax.scatter(bss_col, to_mag, color="green", s=10, label="BSS point", zorder=4)

    # BSS region: bluer and brighter than the BSS point
    xmin, ymin = -5, 0
    ax.plot([bss_col, bss_col], [to_mag, ymin], color="k", lw=1, ls=":")
    ax.plot([bss_col, xmin], [to_mag, to_mag], color="k", lw=1, ls=":")
    ax.fill(
        [xmin, bss_col, bss_col, xmin],
        [ymin, ymin, to_mag, to_mag],
        color="lightblue",
        alpha=1,
        zorder=1,
        label="BSS region",
    )

    ax.legend(fontsize=8)
    ax.tick_params(axis="both", which="major", labelsize=8)
    ax.set_xlabel("Color", fontsize=8)
    ax.set_ylabel("Magnitude", fontsize=8)
    ax.set_xlim(-0.5, 1.2)
    ax.set_ylim(15, 4)

    # plt.show()
    plt.savefig("../bss_region.webp", dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    main()

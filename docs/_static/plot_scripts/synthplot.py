import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import to_rgba

import asteca


def main():
    # Same style as the plots in the tutorials
    plt.style.use("../asteca.mplstyle")

    df = pd.read_csv("../cluster.csv")

    # Isochrones parameters
    isochs = asteca.Isochrones(
        model="parsec",
        isochs_path="../isochrones/parsec",
        mag="Gmag",
        color=("G_BPmag", "G_RPmag"),
        N_points=5000,
        verbose=3,
    )

    # Synthetic clusters parameters
    synthcl = asteca.Synthetic(isochs, ext_law="GAIADR3", seed=142)

    # Generate a `cluster` object
    my_cluster = asteca.Cluster(
        mag=df["Gmag"],
        e_mag=df["e_Gmag"],
        color=df["BP-RP"],
        e_color=df["e_BP-RP"],
        verbose=2,
    )

    # synthcl.calibrate(my_cluster)

    # Generate synthetic cluster.
    synth_clust = synthcl.generate({}, N_stars=2000)
    isoch_arr = synthcl.get_isochrone({})

    binar_idx = ~np.isnan(synth_clust[-1])
    y_synth = synth_clust[0]
    x_synth = synth_clust[1]

    fig, ax = plt.subplots(figsize=(6, 6))
    # Isochrone below the stars (the single systems lie on it, as no uncertainties
    # are added)
    ax.plot(isoch_arr[1], isoch_arr[0], c="k", lw=1, zorder=0, label="Isochrone")

    # Single and binary synthetic systems: semi-transparent points with a thin black
    # edge, as in the tutorials
    for i, (msk, label) in enumerate(
        ((~binar_idx, "Single systems"), (binar_idx, "Binary systems"))
    ):
        ax.scatter(
            x_synth[msk],
            y_synth[msk],
            s=25,
            fc=to_rgba(f"C{i}", 0.25),
            ec="k",
            lw=0.15,
            label=f"{label}, N={msk.sum()}",
        )

    ax.set_xlabel("$G_{BP}-G_{RP}$")
    ax.set_ylabel("$G$")
    ax.legend(loc="upper right")
    ax.invert_yaxis()

    plt.savefig("../synthplot.webp", dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    main()

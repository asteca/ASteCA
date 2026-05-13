import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import asteca


def main():
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

    plt.figure(figsize=(5, 5))

    # Single synthetic systems
    plt.scatter(
        x_synth[~binar_idx],
        y_synth[~binar_idx],
        c="#519ddb",
        alpha=0.35,
        label="Single systems",
    )
    # Binary synthetic systems
    plt.scatter(
        x_synth[binar_idx],
        y_synth[binar_idx],
        c="#F34C4C",
        alpha=0.35,
        label="Binary systems",
    )
    plt.plot(isoch_arr[1], isoch_arr[0], c="k", label="Isochrone")

    plt.xlabel("color")
    plt.ylabel("mag")
    plt.legend()

    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig("../synthplot.webp", dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    main()

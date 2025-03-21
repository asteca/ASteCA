import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import asteca


def main():
    df = pd.read_csv("cluster.csv")

    # Isochrones parameters
    isochs = asteca.isochrones(
        model="parsec",
        isochs_path="parsec",
        magnitude="Gmag",
        color=("G_BPmag", "G_RPmag"),
        magnitude_effl=6390.7,
        color_effl=(5182.58, 7825.08),
        verbose=3,
    )

    # Synthetic clusters parameters
    synthcl = asteca.synthetic(isochs, seed=42)

    # Generate a `cluster` object
    my_cluster = asteca.cluster(
        magnitude=df["Gmag"],
        e_mag=df["e_Gmag"],
        color=df["BP-RP"],
        e_color=df["e_BP-RP"],
        verbose=2,
    )

    synthcl.calibrate(my_cluster)

    fit_params = {
        "loga": 8.0,
        "Av": 0.1,
        "dm": 8,
        "alpha": 0.0,
        "beta": 1.0,
        "Rv": 3.1,
        "met": 0.0152,
        "DR": 0.0,
    }

    # Generate synthetic cluster.
    synth_clust = synthcl.generate(fit_params)
    binar_idx = ~np.isnan(synth_clust[-1])

    y_synth = synth_clust[0]
    x_synth = synth_clust[1]

    plot1(x_synth, y_synth, binar_idx)
    plot2(my_cluster, x_synth, y_synth, binar_idx)
    isoch_arr = synthcl.get_isochrone(fit_params)
    plot3(my_cluster, x_synth, y_synth, binar_idx, isoch_arr)


def plot1(x_synth, y_synth, binar_idx):
    """ """
    plt.figure(figsize=(7, 5))

    # Single synthetic systems
    plt.scatter(
        x_synth[~binar_idx],
        y_synth[~binar_idx],
        marker="^",
        c="#519ddb",
        alpha=0.5,
        label=f"Synthetic single, N={len(x_synth[~binar_idx])}",
    )
    # Binary synthetic systems
    plt.scatter(
        x_synth[binar_idx],
        y_synth[binar_idx],
        marker="v",
        c="#F34C4C",
        alpha=0.5,
        label=f"Synthetic binary, N={len(x_synth[binar_idx])}",
    )

    plt.xlabel("color")
    plt.ylabel("mag")
    plt.legend()

    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig("synthplot.webp", dpi=300, bbox_inches="tight")


def plot2(my_cluster, x_synth, y_synth, binar_idx):
    plt.figure(figsize=(7, 5))

    plt.scatter(
        my_cluster.colors[0],
        my_cluster.mag,
        c="green",
        alpha=0.5,
        label=f"Observed, N={len(my_cluster.mag)}",
    )

    # Single synthetic systems
    plt.scatter(
        x_synth[~binar_idx],
        y_synth[~binar_idx],
        marker="^",
        c="#519ddb",
        alpha=0.5,
        label=f"Synthetic single, N={len(x_synth[~binar_idx])}",
    )
    # Binary synthetic systems
    plt.scatter(
        x_synth[binar_idx],
        y_synth[binar_idx],
        marker="v",
        c="#F34C4C",
        alpha=0.5,
        label=f"Synthetic binary, N={len(x_synth[binar_idx])}",
    )

    plt.xlabel("color")
    plt.ylabel("mag")
    plt.legend()

    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig("obs_synthplot.webp", dpi=300, bbox_inches="tight")


def plot3(my_cluster, x_synth, y_synth, binar_idx, isoch_arr):
    plt.figure(figsize=(7, 5))

    plt.scatter(
        my_cluster.colors[0],
        my_cluster.mag,
        c="green",
        alpha=0.5,
        label=f"Observed, N={len(my_cluster.mag)}",
    )

    # Single synthetic systems
    plt.scatter(
        x_synth[~binar_idx],
        y_synth[~binar_idx],
        marker="^",
        c="#519ddb",
        alpha=0.5,
        label=f"Synthetic single, N={len(x_synth[~binar_idx])}",
    )
    # Binary synthetic systems
    plt.scatter(
        x_synth[binar_idx],
        y_synth[binar_idx],
        marker="v",
        c="#F34C4C",
        alpha=0.5,
        label=f"Synthetic binary, N={len(x_synth[binar_idx])}",
    )

    plt.xlabel("color")
    plt.ylabel("mag")
    plt.legend()

    plt.plot(isoch_arr[1], isoch_arr[0], c="k")

    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig("obs_synthplot_isoch.webp", dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    main()

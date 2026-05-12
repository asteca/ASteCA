import matplotlib.pyplot as plt
import pandas as pd

import asteca


def main():
    df = pd.read_csv("../cluster.csv")

    # Isochrones parameters
    isochs = asteca.Isochrones(
        model="parsec",
        isochs_path="../parsec",
        mag="Gmag",
        color=("G_BPmag", "G_RPmag"),
        magnitude_effl=6390.7,
        color_effl=(5182.58, 7825.08),
        verbose=3,
    )
    # Generate a `cluster` object
    my_cluster = asteca.Cluster(
        mag=df["Gmag"],
        e_mag=df["e_Gmag"],
        color=df["BP-RP"],
        e_color=df["e_BP-RP"],
        verbose=2,
    )

    # # Synthetic clusters parameters
    # synthcl = asteca.Synthetic(isochs, seed=42)
    # synthcl.calibrate(my_cluster)
    # fit_params = {
    #     "loga": 9,
    #     "alpha": 0.0,
    #     "beta": 0.0,
    #     "DR": 0.0,
    # }
    # synth_clust = synthcl.generate(fit_params)
    # y_synth0 = synth_clust[0]
    # x_synth0 = synth_clust[1]

    synthcl = asteca.Synthetic(isochs, seed=42)
    synthcl.calibrate(my_cluster)
    fit_params = {
        "loga": 9,
        "alpha": 0.0,
        "beta": 0.0,
        "Av": 0.5,
        "DR": 1.0,
    }
    synth_clust = synthcl.generate(fit_params)
    isoch_arr1 = synthcl.get_isochrone(fit_params)
    y_synth1 = synth_clust[0]
    x_synth1 = synth_clust[1]

    synthcl = asteca.Synthetic(isochs, DR_distribution="normal", seed=42)
    synthcl.calibrate(my_cluster)
    fit_params = {
        "loga": 9,
        "alpha": 0.0,
        "beta": 0.0,
        "Av": 0.5,
        "DR": 1.0,
    }
    synth_clust = synthcl.generate(fit_params)
    isoch_arr2 = synthcl.get_isochrone(fit_params)
    y_synth2 = synth_clust[0]
    x_synth2 = synth_clust[1]

    fig, axes = plt.subplots(1, 2, figsize=(10, 5), sharey=True)

    # plt.suptitle(r"$A_{v,\mathrm{DR}} = \max\left(0,\; X\right)$")

    # axes[0].set_title("No DR", fontsize=12)
    # axes[0].scatter(x_synth0, y_synth0, alpha=0.25)

    # axes[0].set_title(r"Uniform: $X \sim \mathcal{U}(A_v-DR, A_v+DR)$", fontsize=12)
    axes[0].set_title("Uniform")
    axes[0].scatter(x_synth1, y_synth1, alpha=0.25)
    axes[0].plot(isoch_arr1[1], isoch_arr1[0], lw=2, zorder=5, c="r")
    axes[0].set_xlabel("color")
    axes[0].set_ylabel("magnitude")

    # axes[1].set_title(r"Normal: $X \sim \mathcal{N}(A_v, DR)$", fontsize=12)
    axes[1].set_title("Normal")
    axes[1].scatter(
        x_synth2,
        y_synth2,
        alpha=0.25,
    )
    axes[1].plot(isoch_arr2[1], isoch_arr2[0], lw=2, zorder=5, c="r")
    axes[1].set_xlabel("color")

    # plt.legend()

    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig("../dr_plot.webp", dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    main()

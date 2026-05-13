import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Import local version
import asteca

print(f"ASteCA version: {asteca.__version__}")


def main(seed=142):
    """
    Generate random synthetic clusters from the three supported isochrone services
    """
    path = "../cluster.csv"
    obs_df = pd.read_csv(path)
    my_cluster = asteca.Cluster(
        ra=obs_df["RA_ICRS"],
        dec=obs_df["DE_ICRS"],
        mag=obs_df["Gmag"],
        e_mag=obs_df["e_Gmag"],
        color=obs_df["BP-RP"],
        e_color=obs_df["e_BP-RP"],
    )

    synthcl_dict = {}
    for i, isoch_model in enumerate(("parsec", "basti", "mist")):
        print("\nLoading:", isoch_model)
        # Isochrones parameters
        if i == 0:
            isochs = asteca.Isochrones(
                model="parsec",
                isochs_path="../isochrones/parsec",
                mag="Gmag",
                color=("G_BPmag", "G_RPmag"),
                verbose=3,
            )
        elif i == 1:
            isochs = asteca.Isochrones(
                model="Basti",
                isochs_path="../isochrones/basti",
                mag="G",
                color=("G_BP", "G_RP"),
                verbose=3,
            )
        elif i == 2:
            isochs = asteca.Isochrones(
                model="MIST",
                isochs_path="../isochrones/mist",
                mag="Gaia_G_EDR3",
                color=("Gaia_BP_EDR3", "Gaia_RP_EDR3"),
                verbose=3,
            )

        synthcl = asteca.Synthetic(isochs, ext_law="GAIADR3", verbose=3, seed=seed)
        synthcl.calibrate(my_cluster)
        synthcl_dict[isoch_model] = synthcl

    xmin, xmax = np.inf, -np.inf
    ymin, ymax = -np.inf, np.inf

    fig, axes = plt.subplots(1, 3, figsize=(10, 5), sharey=True)

    for i, (isoch_model, synthcl) in enumerate(synthcl_dict.items()):
        ax = axes[i]
        dm, Av = 0, 0

        if isoch_model == "basti":
            met = 0.01258
        elif isoch_model == "mist":
            met = 0.0127212
        else:
            met = 0.0152
        for loga in (7.5, 8.5, 9.5):
            model = {
                "met": met,
                "loga": loga,
                "dm": dm,
                "Av": Av,
            }
            isoch_arr = synthcl.get_isochrone(model, full_track=True)
            xmin, xmax = min(xmin, min(isoch_arr[1])), max(xmax, max(isoch_arr[1]))
            ymin, ymax = max(ymin, max(isoch_arr[0])), min(ymax, min(isoch_arr[0]))
            # plot_synthetic(synthcl, ax, model, isoch_arr, title=isoch_model)
            ax.plot(isoch_arr[1], isoch_arr[0], lw=2, label=f"log(age)={loga:.1f}")

            ax.set_title(str(isoch_model).upper())
            ax.set_xlabel("color")
            if i == 1:
                ax.legend()  # loc="upper left")

    xmin, xmax = xmin - xmax * 0.1, xmax + xmax * 0.02
    ymin, ymax = ymin + ymin * 0.02, ymax - (ymin - ymax) * 0.05
    ax = axes[0]
    ax.set_ylabel("magnitude")
    # ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax = axes[1]
    # ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax = axes[2]
    # ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    fig.tight_layout()
    plt.savefig("../isoch_models.webp", dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    main()

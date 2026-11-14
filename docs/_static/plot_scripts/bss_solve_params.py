import matplotlib.pyplot as plt
import pandas as pd

import asteca


def main():
    df_bss = pd.read_csv("../cluster_bss.csv")
    # Generate a `cluster` object
    my_cluster = asteca.Cluster(
        mag=df_bss["Gmag"],
        e_mag=df_bss["e_Gmag"],
        color=df_bss["BP-RP"],
        e_color=df_bss["e_BP-RP"],
    )

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
        verbose=3,
    )

    synthcl.calibrate(my_cluster)

    priors = {
        "loga": (isochs.amin, isochs.amax),
        "dm": (7.5, 8.5),
        "Av": (0.0, 5.0),
    }
    sampler = asteca.Sampler(synthcl, priors, temperature="auto")
    sampler.run()
    summ = sampler.summary()
    print(f"Best distance: {sampler.best_dist:.3f}")
    print(f"Estimated parameters: {summ['median']}")
    model_means, model_std = summ["median"], summ["std"]

    props = synthcl.cluster_properties(model_means, model_std)
    bss_probs = props["stars"]["bss_prob"]
    N_bss = (bss_probs > 0.5).sum()
    print(f"Estimated number of BSS: {N_bss}")

    make_plots(my_cluster, synthcl, model_means, bss_probs)


def make_plots(my_cluster, synthcl, model_means, bss_probs):

    isoch_arr = synthcl.get_isochrone(model_means)

    plt.scatter(
        my_cluster.color,
        my_cluster.mag,
        # color="blue",
        # s=10,
    )
    bss_msk = bss_probs > 0.5
    plt.scatter(
        my_cluster.color[bss_msk],
        my_cluster.mag[bss_msk],
        color="red",
        # s=10,
    )
    plt.plot(isoch_arr[1], isoch_arr[0], color="black")
    plt.gca().invert_yaxis()
    plt.show()


if __name__ == "__main__":
    main()

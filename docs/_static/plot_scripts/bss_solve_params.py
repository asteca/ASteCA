import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import differential_evolution

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

    likelihood = asteca.Likelihood(my_cluster)

    loga_min, loga_max = isochs.amin, isochs.amax
    dm_min, dm_max = 7.5, 8.5

    loga, dm, Av, final_dist = first_pass_DE(
        synthcl, likelihood, loga_min, loga_max, dm_min, dm_max
    )
    print(f"Final distance: {final_dist:.3f}")
    print(f"Estimated parameters: loga={loga:.3f}, dm={dm:.3f}, Av={Av:.3f}")
    model_means = {"loga": loga, "dm": dm, "Av": Av}
    model_std = {"loga": 0.1, "dm": 0.1, "Av": 0.05}

    synthcl.get_models(model_means, model_std)
    stellar_data = synthcl.stellar_parameters()
    bss_probs = stellar_data["bss_prob"]
    N_bss = (bss_probs > 0.5).sum()
    print(f"Estimated number of BSS: {N_bss}")

    make_plots(my_cluster, synthcl, model_means, bss_probs)


def first_pass_DE(synthcl, likelihood, loga_min, loga_max, dm_min, dm_max):
    """ """
    bounds = (
        (loga_min, loga_max),  # loga
        (dm_min, dm_max),  # dn
        (0.0, 5.0),  # Av
    )
    result = differential_evolution(
        objective,
        bounds,
        popsize=50,
        recombination=0.8,
        args=(synthcl, likelihood),
        polish=True,
        seed=42,
    )

    loga, dm, Av = result.x
    return loga, dm, Av, result.fun


def objective(params, synthcl, likelihood):
    loga, dm, Av = params
    fit_params = {
        "loga": loga,
        "dm": dm,
        "Av": Av,
    }
    synth_data = synthcl.generate(fit_params)
    lkl = likelihood.get(synth_data)
    return lkl


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

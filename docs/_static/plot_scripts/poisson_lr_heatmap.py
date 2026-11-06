import asteca
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm

"""
Distance d returned by Likelihood.get() over a (loga, dm) grid, for the observed
cluster used in the parameter fitting tutorials. The 'met' and 'Av' parameters are
fixed to the values estimated in the pocoMC tutorial.
"""

# pocoMC tutorial estimates
fit = {"met": 0.014, "loga": 7.962, "dm": 8.181, "Av": 0.482}


def main():
    isochs = asteca.Isochrones(
        model="parsec",
        isochs_path="../isochrones/parsec/",
        mag="Gmag",
        color=("G_BPmag", "G_RPmag"),
        magnitude_effl=6390.7,
        color_effl=(5182.58, 7825.08),
        verbose=0,
    )
    synthcl = asteca.Synthetic(isochs, seed=457304, verbose=0)
    obs_df = pd.read_csv("../cluster.csv")
    my_cluster = asteca.Cluster(
        ra=obs_df["RA_ICRS"],
        dec=obs_df["DE_ICRS"],
        mag=obs_df["Gmag"],
        e_mag=obs_df["e_Gmag"],
        color=obs_df["BP-RP"],
        e_color=obs_df["e_BP-RP"],
        verbose=0,
    )
    synthcl.calibrate(my_cluster)
    likelihood = asteca.Likelihood(my_cluster)

    # Full range of the priors used in the tutorials, and a zoom around the estimates
    grids = (
        (np.linspace(7.0, 9.5, 101), np.linspace(8.0, 10.5, 101)),
        (
            np.linspace(fit["loga"] - 0.3, fit["loga"] + 0.3, 61),
            np.linspace(fit["dm"] - 0.3, fit["dm"] + 0.3, 61),
        ),
    )

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
    for ax, (loga_arr, dm_arr) in zip(axes, grids):
        d = np.empty((len(dm_arr), len(loga_arr)))
        for j, dm in enumerate(dm_arr):
            for i, loga in enumerate(loga_arr):
                params = dict(fit, loga=loga, dm=dm)
                d[j, i] = likelihood.get(synthcl.generate(params))

        extent = (
            loga_arr[0],
            loga_arr[-1],
            dm_arr[0],
            dm_arr[-1],
        )
        im = ax.imshow(
            d,
            origin="lower",
            extent=extent,
            aspect="auto",
            cmap="viridis",
            norm=LogNorm(vmin=max(d.min(), 1e-2), vmax=d.max()),
            interpolation="nearest",
        )
        cb = fig.colorbar(im, ax=ax)
        cb.set_label("d")
        if d.max() > 1:
            # Boundary between models with and without overlap in the Hess diagram
            ax.contour(
                loga_arr, dm_arr, d, levels=[1.0], colors="w", linewidths=1
            )
        ax.scatter(fit["loga"], fit["dm"], marker="x", c="r", s=60, zorder=5)
        ax.set_xlabel("loga")
        ax.set_ylabel("dm")
        print(f"d range: [{d.min():.3f}, {d.max():.3f}], N(d>1)={(d > 1).sum()}")

    axes[0].set_title("Prior range")
    axes[1].set_title("Zoom around the estimates")
    fig.tight_layout()
    plt.savefig("../lkl_heatmap.webp", dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    main()

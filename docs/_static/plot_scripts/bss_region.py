import matplotlib.pyplot as plt
from asteca.modules import stellar_stats_funcs as ssf

import asteca


def main():
    # Isochrones parameters
    isochs = asteca.Isochrones(
        model="parsec",
        isochs_path="../../../../isochrones/gaia_solar_range",
        mag="Gmag",
        color=("G_BPmag", "G_RPmag"),
        verbose=3,
    )

    # Synthetic clusters parameters
    synthcl = asteca.Synthetic(
        isochs,
        ext_law="GAIADR3",
        max_mass=40000,
        seed=42519,
        verbose=3,
    )

    model = {
        "loga": 8.0,
        "dm": 8.1,
        "Av": 0.42,
    }

    isoch_arr = synthcl.get_isochrone(model)
    cluster_mag, cluster_color = isoch_arr[1], isoch_arr[0]

    # --- Find TO region ---
    to_col_1, to_col_2, to_mag_1, to_mag_2 = ssf.to_point_find(
        isoch_arr[1],
        isoch_arr[0],
        cluster_mag,
        cluster_color,
        isochs_model="parsec",
        mag_offset_1=0.5,
        col_offset_1=-0.1,
        mag_offset_2=-0.5,
        col_offset_2=-0.05,
    )

    mplot(isoch_arr, to_col_1, to_col_2, to_mag_1, to_mag_2)


def mplot(isoch_arr, to_col_1, to_col_2, to_mag_1, to_mag_2):

    fig, ax = plt.subplots(figsize=(4, 4))

    # Isochrone
    ax.plot(isoch_arr[1], isoch_arr[0], color="black")

    ax.scatter(to_col_1, to_mag_1, color="green", s=10, label="TO Point 1", zorder=4)
    ax.scatter(to_col_2, to_mag_2, color="red", s=10, label="TO Point 2", zorder=4)

    ax.plot([to_col_2, to_col_2], [to_mag_2 - 5, to_mag_2], color="k", lw=1, ls=":")
    ax.plot([to_col_2, to_col_1], [to_mag_2, to_mag_2], color="k", lw=1, ls=":")
    ax.plot([to_col_1, to_col_1], [to_mag_1, to_mag_2], color="k", lw=1, ls=":")
    ax.plot([to_col_1, to_col_1 - 5], [to_mag_1, to_mag_1], color="k", lw=1, ls=":")

    # Define polygon vertices (following the drawn boundary)
    x_poly = [
        to_col_2,     # top-left
        to_col_1,     # top-right
        to_col_1,     # down to (to_mag_1)
        -5,           # extend left to x = -5
        -5,           # down to y = 0
        to_col_2      # back toward left boundary
    ]

    y_poly = [
        to_mag_2,     # top
        to_mag_2,
        to_mag_1,
        to_mag_1,
        0,
        to_mag_2 - 5
    ]
    # Fill region
    ax.fill(x_poly, y_poly, color="lightblue", alpha=1, zorder=1, label="BSS Region")

    ax.legend(fontsize=8)
    # smaller x and y ticks
    ax.tick_params(axis="both", which="major", labelsize=8)
    ax.set_xlabel("Color", fontsize=8)
    ax.set_ylabel("Magnitude", fontsize=8)
    ax.set_xlim(-.5, 1.2)
    ax.set_ylim(15, 4)

    # plt.show()
    plt.savefig("../bss_region.webp", dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    main()

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import asteca

"""
Compare the CMDs of the fastMP members (P>=0.5) obtained with and without
'mag_correct=True', for the UCC (https://ucc.ar) cluster Kronberger 31.

The Gaia DR3 frame around the cluster (G<20) is read from
'../kronberger31_field.csv'. The center values and number of members are fixed to
the values listed in the UCC catalogue (version 261022).
"""

seed = 457304

# UCC values for Kronberger 31
radec_c = (295.01889, 26.25776)
pms_c = (-2.3499, -4.4781)
plx_c = 0.0594
N_cluster = 125


def main():
    # Same style as the plots in the tutorials
    plt.style.use("../asteca.mplstyle")

    df = pd.read_csv("../kronberger31_field.csv")
    my_field = asteca.Cluster(
        ra=df["RA_ICRS"],
        dec=df["DE_ICRS"],
        mag=df["Gmag"],
        pmra=df["pmRA"],
        pmde=df["pmDE"],
        plx=df["Plx"],
        e_mag=df["e_Gmag"],
        e_pmra=df["e_pmRA"],
        e_pmde=df["e_pmDE"],
        e_plx=df["e_Plx"],
        verbose=0,
    )
    my_field.radec_c = radec_c
    my_field.pms_c = pms_c
    my_field.plx_c = plx_c
    my_field.N_cluster = N_cluster

    my_field.membership.fastmp(seed=seed)
    probs_0 = my_field.probs.copy()
    my_field.membership.fastmp(seed=seed, mag_correct=True)
    probs_1 = my_field.probs.copy()

    # The color is not used by fastMP, only to draw the CMDs
    mag, color = np.asarray(df["Gmag"]), np.asarray(df["BP-RP"])

    fig, axes = plt.subplots(1, 2, figsize=(11, 5), sharex=True, sharey=True)
    titles = ("No correction", "mag_correct=True")
    for ax, probs, title in zip(axes, (probs_0, probs_1), titles):
        ax.set_box_aspect(1)
        msk = probs >= 0.5
        ax.scatter(color, mag, c="grey", s=4, alpha=0.15, lw=0)
        im = ax.scatter(
            color[msk],
            mag[msk],
            c=probs[msk],
            s=25,
            ec="k",
            lw=0.15,
            vmin=0.5,
            vmax=1,
        )
        ax.set_title(f"{title} ({msk.sum()} stars with P>=0.5)")
        ax.set_xlabel("$G_{BP}-G_{RP}$")
    axes[0].set_ylabel("$G$")
    axes[0].set_xlim(0, 3.5)
    axes[0].set_ylim(20.2, 13)
    # Colorbar with the same height as the (square) panels
    for ax in axes:
        ax.apply_aspect()
    pos = axes[1].get_position()
    cax = fig.add_axes((pos.x1 + 0.015, pos.y0, 0.015, pos.height))
    fig.colorbar(im, cax=cax, label="Membership probability")
    fig.suptitle("Kronberger 31", y=pos.y1 + 0.13)

    plt.savefig("../mag_correct.webp", dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    main()

import asteca
import asteca.modules.fastmp as fmp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

"""
Compare the fastMP probabilities obtained with and without 'mag_correct=True',
for the field used in the membership notebook.
"""

seed = 457304

# Store the reference selection (no error resampling) used by the magnitude correction
_ref = {}
_mag_dist_correct = fmp.mag_dist_correct


def _capture(probs, ref_msk, ref_dist, mag_clean, *args, **kwargs):
    _ref["msk"], _ref["mag"] = ref_msk.copy(), mag_clean.copy()
    return _mag_dist_correct(probs, ref_msk, ref_dist, mag_clean, *args, **kwargs)


fmp.mag_dist_correct = _capture


def main():
    df = pd.read_csv("../field.csv")
    my_field = asteca.Cluster(
        ra=df["RA_ICRS"],
        dec=df["DE_ICRS"],
        mag=df["Gmag"],
        color=df["BP-RP"],
        pmra=df["pmRA"],
        pmde=df["pmDE"],
        plx=df["Plx"],
        e_mag=df["e_Gmag"],
        e_color=df["e_BP-RP"],
        e_pmra=df["e_pmRA"],
        e_pmde=df["e_pmDE"],
        e_plx=df["e_Plx"],
        verbose=0,
    )
    my_field.get_center()
    my_field.get_nmembers()

    my_field.membership.fastmp(seed=seed)
    probs_0 = my_field.probs.copy()
    my_field.membership.fastmp(seed=seed, mag_correct=True)
    probs_1 = my_field.probs.copy()

    mag = np.asarray(df["Gmag"])
    msk_0, msk_1 = probs_0 >= 0.5, probs_1 >= 0.5
    promoted, demoted = msk_1 & ~msk_0, msk_0 & ~msk_1
    mag_ref = _ref["mag"][_ref["msk"]]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))

    # Left: probabilities versus magnitude
    show = probs_0 > 0.05
    ax1.scatter(
        mag[show], probs_0[show], s=12, c="grey", alpha=0.5, lw=0, label="No correction"
    )
    ax1.scatter(
        mag[promoted],
        probs_1[promoted],
        marker="^",
        s=40,
        c="C0",
        ec="k",
        lw=0.5,
        label=f"Promoted ({promoted.sum()})",
    )
    ax1.scatter(
        mag[demoted],
        probs_1[demoted],
        marker="v",
        s=40,
        c="C3",
        ec="k",
        lw=0.5,
        label=f"Demoted ({demoted.sum()})",
    )
    for i in np.flatnonzero(promoted | demoted):
        ax1.plot([mag[i], mag[i]], [probs_0[i], probs_1[i]], c="k", lw=0.5, alpha=0.5)
    ax1.axhline(0.5, ls=":", c="k", lw=1)
    ax1.set_xlabel("G")
    ax1.set_ylabel("P")
    ax1.legend(loc="lower left")

    # Right: magnitude distribution of the P>=0.5 members
    edges = np.arange(np.floor(np.nanmin(mag_ref)), np.ceil(np.nanmax(mag)) + 0.5, 0.5)
    ax2.hist(
        mag_ref,
        edges,
        histtype="stepfilled",
        color="grey",
        alpha=0.3,
        label="Reference (no resampling)",
    )
    ax2.hist(
        mag[msk_0],
        edges,
        histtype="step",
        lw=2,
        color="C1",
        label=f"P>=0.5, no correction ({msk_0.sum()})",
    )
    ax2.hist(
        mag[msk_1],
        edges,
        histtype="step",
        lw=2,
        ls="--",
        color="C0",
        label=f"P>=0.5, mag_correct ({msk_1.sum()})",
    )
    ax2.set_xlim(edges[0], np.nanmax(mag_ref) + 1)
    ax2.set_ylim(0, ax2.get_ylim()[1] * 1.3)
    ax2.set_xlabel("G")
    ax2.set_ylabel("N")
    ax2.legend(loc="upper left")

    fig.tight_layout()
    plt.savefig("../mag_correct.webp", dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    main()
